#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* Priority levels: 0 - highest (audio), 1 - next (video), 2 - lowest (data) */
#define PRIORITY_LEVELS 3

/* Event types */
#define EVENT_ARRIVAL 0
#define EVENT_DEPARTURE 1

enum PacketType { AUDIO = 0, VIDEO, DATA, REFERENCE };


typedef struct Packet {
    int type;             // Packet type: AUDIO, VIDEO, DATA or REFERENCE
    double arrival_time;
    double service_time;
    int is_reference;
    double creation_time;
    struct Packet *next;
} Packet;

/* FIFO queue structure */
typedef struct Queue {
    Packet *front;
    Packet *rear;
    int count;
} Queue;


typedef struct Node {
    Queue queues[PRIORITY_LEVELS];
    int capacity;       // Total queue capacity (K)
    int current_count;
    int busy;
} Node;

typedef struct Event {
    int type;         // Event type: arrival or departure
    double time;      // Time at which the event occurs
    int node_id;
    Packet *packet;
    struct Event *next;
} Event;

/* Global simulation time */
double sim_time = 0.0;

/* Event list: sorted in ascending order by event time */
Event *event_list = NULL;

/* Global Allocated Variables for Statistics */
static double **sum_delay;
static int    **count_served;
static int    **count_arrivals;
static int    **count_drops;
static double **area_queue;

/* Statistics for the reference flow */
static int    count_ref_arrivals = 0;
static int    count_ref_drops = 0;
static double sum_end_to_end_ref_delay = 0;
static int    count_ref_served = 0;

/* Last event time used for time-average queue length calculation */
static double last_event_time = 0.0;

void initQueue(Queue *q) {
    q->front = q->rear = NULL;
    q->count = 0;
}

void enqueue(Queue *q, Packet *p) {
    p->next = NULL;
    if(q->rear == NULL) {
        q->front = q->rear = p;
    } else {
        q->rear->next = p;
        q->rear = p;
    }
    q->count++;
}

Packet* dequeue(Queue *q) {
    if(q->front == NULL)
        return NULL;
    Packet *p = q->front;
    q->front = q->front->next;
    if(q->front == NULL)
        q->rear = NULL;
    q->count--;
    return p;
}

/* Insert a new event into the event list in order of time */
void schedule_event(Event *new_event) {
    if(event_list == NULL || new_event->time < event_list->time) {
        new_event->next = event_list;
        event_list = new_event;
    } else {
        Event *temp = event_list;
        while(temp->next != NULL && temp->next->time < new_event->time) {
            temp = temp->next;
        }
        new_event->next = temp->next;
        temp->next = new_event;
    }
}

/* Remove and return the earliest event */
Event* get_next_event() {
    if(event_list == NULL)
        return NULL;
    Event *e = event_list;
    event_list = event_list->next;
    return e;
}

/* Generate an exponentially distributed random number with the given mean */
double exp_random(double mean) {
    double u = (rand() + 1.0) / (RAND_MAX + 1.0);
    return -mean * log(u);
}

/* Compute the service time based on link capacity C (bps) and packet size (bytes) */
double compute_service_time(int packet_size, double C) {
    return (packet_size * 8.0) / C;
}

/* Before processing each event, update the accumulated area under the queue length curve for each node based on the time interval since the last event */
void update_area_queues(double time_diff, Node *nodes, int M) {
    for(int i = 0; i < M; i++){
        for(int p = 0; p < PRIORITY_LEVELS; p++){
            area_queue[i][p] += nodes[i].queues[p].count * time_diff;
        }
    }
}

void process_arrival(Event *e, Node *nodes, int M, double C) {
    int node_id = e->node_id;
    Node *node = &nodes[node_id];
    Packet *pkt = e->packet;

    int queue_index;
    if(pkt->type == AUDIO || pkt->type == REFERENCE) {
        queue_index = 0;  // Premium
    } else if(pkt->type == VIDEO) {
        queue_index = 1;  // Assured
    } else {
        queue_index = 2;  // Best-Effort
    }

    count_arrivals[node_id][queue_index]++;

    /* Check if the node's queue is full */
    if(node->current_count >= node->capacity) {
        // Update drop statistics
        count_drops[node_id][queue_index]++;
        if(pkt->is_reference) {
            count_ref_drops++;
        }
        // Packet is dropped
        free(pkt);
    } else {
        enqueue(&node->queues[queue_index], pkt);
        node->current_count++;

        /* If the node is currently idle, schedule a departure event */
        if(!node->busy) {
            Packet *next_pkt = NULL;
            for(int i = 0; i < PRIORITY_LEVELS; i++) {
                if(node->queues[i].count > 0) {
                    next_pkt = node->queues[i].front;
                    break;
                }
            }
            if(next_pkt) {
                Event *de = (Event*)malloc(sizeof(Event));
                de->type = EVENT_DEPARTURE;
                de->time = sim_time + next_pkt->service_time;
                de->node_id = node_id;
                de->packet = next_pkt;
                de->next = NULL;
                schedule_event(de);
                node->busy = 1;
            }
        }
    }
}

void process_departure(Event *e, Node *nodes, int M, double C) {
    int node_id = e->node_id;
    Node *node = &nodes[node_id];

    /* Find the queue where the packet is located (assumed to be at the front) */
    int queue_index = -1;
    for(int i = 0; i < PRIORITY_LEVELS; i++){
        if(node->queues[i].front == e->packet) {
            queue_index = i;
            break;
        }
    }
    if(queue_index >= 0) {
        Packet *pkt = dequeue(&node->queues[queue_index]);
        node->current_count--;

        /* Calculate the time the packet spent in the node */
        double delay = sim_time - pkt->arrival_time;
        /* Accumulate delay and update served packet count */
        sum_delay[node_id][queue_index] += delay;
        count_served[node_id][queue_index]++;

        /* If it is a reference flow packet and it has not reached the final node, schedule an arrival event at the next node immediately */
        if(pkt->is_reference && node_id < M - 1) {
            Event *new_arrival = (Event*)malloc(sizeof(Event));
            new_arrival->type = EVENT_ARRIVAL;
            new_arrival->time = sim_time; // Propagation delay can be added here if needed
            new_arrival->node_id = node_id + 1;

            /* Set arrival_time for the next node to the current time */
            pkt->arrival_time = sim_time;
            new_arrival->packet = pkt;
            new_arrival->next = NULL;
            schedule_event(new_arrival);
        } else if(pkt->is_reference && node_id == M - 1) {
            /* For a reference flow packet reaching the final node, update end-to-end delay statistics */
            double end_to_end = sim_time - pkt->creation_time;
            sum_end_to_end_ref_delay += end_to_end;
            count_ref_served++;
            free(pkt);
        } else {
            /* For non-reference flow packets or those already forwarded, free the packet */
            free(pkt);
        }
    }

    /* Check if there are more packets waiting in the node; if so, schedule the next departure event */
    Packet *next_pkt = NULL;
    for(int i = 0; i < PRIORITY_LEVELS; i++){
        if(node->queues[i].count > 0) {
            next_pkt = node->queues[i].front;
            break;
        }
    }
    if(next_pkt) {
        Event *de = (Event*)malloc(sizeof(Event));
        de->type = EVENT_DEPARTURE;
        de->time = sim_time + next_pkt->service_time;
        de->node_id = node_id;
        de->packet = next_pkt;
        de->next = NULL;
        schedule_event(de);
        node->busy = 1;
    } else {
        node->busy = 0;
    }
}

/*
  Command-line arguments:
    argv[1]: Number of nodes (M)
    argv[2]: Queue capacity per node (K)
    argv[3]: Number of audio sources
    argv[4]: Number of video sources
    argv[5]: Number of data sources
    argv[6]: Reference flow type (0: audio, 1: video, 2: data)
*/
int main(int argc, char *argv[]) {
    if(argc < 7) {
        printf("Usage: %s M K N_audio N_video N_data ref_type\n", argv[0]);
        return 1;
    }
    int M = atoi(argv[1]);
    int K = atoi(argv[2]);
    int N_audio = atoi(argv[3]);
    int N_video = atoi(argv[4]);
    int N_data = atoi(argv[5]);
    int ref_type = atoi(argv[6]); // Reference flow type: 0 - audio, 1 - video, 2 - data
    double C = 10e6; // Link capacity: 10 Mbps

    /* Allocate and initialize M nodes */
    Node *nodes = (Node*)malloc(M * sizeof(Node));
    for(int i = 0; i < M; i++){
        nodes[i].capacity = K;
        nodes[i].current_count = 0;
        nodes[i].busy = 0;
        for(int j = 0; j < PRIORITY_LEVELS; j++){
            initQueue(&nodes[i].queues[j]);
        }
    }
    srand(time(NULL));

    /* Initialize statistics variables */
    sum_delay    = (double**)malloc(M * sizeof(double*));
    count_served = (int**)malloc(M * sizeof(int*));
    count_arrivals = (int**)malloc(M * sizeof(int*));
    count_drops  = (int**)malloc(M * sizeof(int*));
    area_queue   = (double**)malloc(M * sizeof(double*));
    for(int i = 0; i < M; i++){
        sum_delay[i]    = (double*)calloc(PRIORITY_LEVELS, sizeof(double));
        count_served[i] = (int*)calloc(PRIORITY_LEVELS, sizeof(int));
        count_arrivals[i] = (int*)calloc(PRIORITY_LEVELS, sizeof(int));
        count_drops[i]  = (int*)calloc(PRIORITY_LEVELS, sizeof(int));
        area_queue[i]   = (double*)calloc(PRIORITY_LEVELS, sizeof(double));
    }

    /* Compute effective arrival rates for each packet type */
    double lambda_audio = (64e3 / (120 * 8.0)) * (0.36 / (0.36 + 0.64));
    double lambda_video = (384e3 / (1000 * 8.0)) * (0.33 / (0.33 + 0.73));
    double lambda_data  = (256e3 / (583 * 8.0)) * (0.35 / (0.35 + 0.65));

    /* Schedule initial arrival events for background traffic for each node */
    for(int i = 0; i < M; i++){
        /* Audio source */
        for(int s = 0; s < N_audio; s++){
            double t = exp_random(1.0 / lambda_audio);
            Packet *p = (Packet*)malloc(sizeof(Packet));
            p->type = AUDIO;
            p->arrival_time = t;
            p->service_time = compute_service_time(120, C);
            p->is_reference = 0;
            p->creation_time = 0.0;
            p->next = NULL;
            Event *ev = (Event*)malloc(sizeof(Event));
            ev->type = EVENT_ARRIVAL;
            ev->time = t;
            ev->node_id = i;
            ev->packet = p;
            ev->next = NULL;
            schedule_event(ev);
        }
        /* Video source */
        for(int s = 0; s < N_video; s++){
            double t = exp_random(1.0 / lambda_video);
            Packet *p = (Packet*)malloc(sizeof(Packet));
            p->type = VIDEO;
            p->arrival_time = t;
            p->service_time = compute_service_time(1000, C);
            p->is_reference = 0;
            p->creation_time = 0.0;
            p->next = NULL;
            Event *ev = (Event*)malloc(sizeof(Event));
            ev->type = EVENT_ARRIVAL;
            ev->time = t;
            ev->node_id = i;
            ev->packet = p;
            ev->next = NULL;
            schedule_event(ev);
        }
        /* Data source */
        for(int s = 0; s < N_data; s++){
            double t = exp_random(1.0 / lambda_data);
            Packet *p = (Packet*)malloc(sizeof(Packet));
            p->type = DATA;
            p->arrival_time = t;
            p->service_time = compute_service_time(583, C);
            p->is_reference = 0;
            p->creation_time = 0.0;
            p->next = NULL;
            Event *ev = (Event*)malloc(sizeof(Event));
            ev->type = EVENT_ARRIVAL;
            ev->time = t;
            ev->node_id = i;
            ev->packet = p;
            ev->next = NULL;
            schedule_event(ev);
        }
    }

    /* Schedule the initial arrival event for the reference flow (starting at node 0) */
    double t_ref = exp_random(1.0 / (ref_type == 0 ? lambda_audio :
                                     (ref_type == 1 ? lambda_video : lambda_data)));
    Packet *p_ref = (Packet*)malloc(sizeof(Packet));
    p_ref->type = ref_type;  // Reference flow uses the corresponding type
    p_ref->arrival_time = t_ref;
    p_ref->creation_time = t_ref;  // For end-to-end delay calculation
    int packet_size = (ref_type == 0) ? 120 : (ref_type == 1 ? 1000 : 583);
    p_ref->service_time = compute_service_time(packet_size, C);
    p_ref->is_reference = 1;
    p_ref->next = NULL;

    count_ref_arrivals++;

    Event *ev_ref = (Event*)malloc(sizeof(Event));
    ev_ref->type = EVENT_ARRIVAL;
    ev_ref->time = t_ref;
    ev_ref->node_id = 0;
    ev_ref->packet = p_ref;
    ev_ref->next = NULL;
    schedule_event(ev_ref);


    while(event_list != NULL) {
        Event *e = get_next_event();
        // Before processing each event, update the time-integral of the queue lengths
        double time_diff = e->time - last_event_time;
        if(time_diff > 0) {
            update_area_queues(time_diff, nodes, M);
            last_event_time = e->time;
        }

        sim_time = e->time; // Current simulation time
        if(e->type == EVENT_ARRIVAL) {
            process_arrival(e, nodes, M, C);
        } else if(e->type == EVENT_DEPARTURE) {
            process_departure(e, nodes, M, C);
        }
        free(e);
    }

    printf("Simulation ended at time %.3f\n\n", sim_time);

    /* -- Results -- */
    /* (a), (b), (c): Per node/queue statistics */
    for(int i = 0; i < M; i++){
        for(int p = 0; p < PRIORITY_LEVELS; p++){
            double avg_delay = 0.0;
            if(count_served[i][p] > 0) {
                avg_delay = sum_delay[i][p] / count_served[i][p];
            }
            double blocking_ratio = 0.0;
            int total_arrivals = count_arrivals[i][p];
            if(total_arrivals > 0) {
                blocking_ratio = (double)count_drops[i][p] / (double)total_arrivals;
            }
            double avg_queue_length = 0.0;
            if(sim_time > 0) {
                avg_queue_length = area_queue[i][p] / sim_time;
            }
            printf("Node %d, Priority %d: avg_delay=%.6f s, blocking_ratio=%.4f, avg_queue_len=%.4f\n",
                   i, p, avg_delay, blocking_ratio, avg_queue_length);
        }
    }

    /* (d) and (e): Reference flow statistics */
    double ref_avg_delay = 0.0;
    if(count_ref_served > 0) {
        ref_avg_delay = sum_end_to_end_ref_delay / count_ref_served;
    }
    double ref_block_ratio = 0.0;
    if(count_ref_arrivals > 0) {
        ref_block_ratio = (double)count_ref_drops / (double)count_ref_arrivals;
    }
    printf("\nReference flow: end-to-end avg_delay=%.6f s, block_ratio=%.4f\n",
           ref_avg_delay, ref_block_ratio);

    free(nodes);
    for(int i = 0; i < M; i++){
        free(sum_delay[i]);
        free(count_served[i]);
        free(count_arrivals[i]);
        free(count_drops[i]);
        free(area_queue[i]);
    }
    free(sum_delay);
    free(count_served);
    free(count_arrivals);
    free(count_drops);
    free(area_queue);

    return 0;
}