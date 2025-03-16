# Packet Switching Simulation
This program simulates a network of nodes using a strict priority queue scheduling discipline. It models different types of traffic (audio, video, data) along with a reference flow. The simulation gathers statistics such as average packet delay, packet drop (blocking) ratios, and average queue lengths per node and per priority level, as well as end-to-end statistics for the reference flow.

## Requirements
	• C/C++ compiler (e.g., GCC or G++)

## Compilation

To compile the program, open a terminal in the directory containing the source file (e.g., main.cpp) and run:

```
gcc main.cpp -o simulation
```

## Usage

The program is executed from the command line with the following parameters:
```
./simulation M K N_audio N_video N_data ref_type
```
### Parameters
	• M: Number of nodes in the network
	• K: Queue capacity per node
	• N_audio: Number of audio sources
	• N_video: Number of video sources
	• N_data: Number of data sources
	• ref_type: Reference flow type: 0 for audio, 1 for video, 2 for data

### Example
To simulate a network with 5 nodes, each with a queue capacity of 100, 2 audio sources, 2 video sources, 2 data sources, and a reference flow set to audio, run:
```
./simulation 5 100 2 2 2 0
```
### Output
 ##### After running, the program will process all scheduled events and then output statistics, including:
	• Per Node/Queue Statistics:
	    • Average packet delay at each node and for each priority.
        • Packet blocking ratio (i.e., the ratio of dropped packets to total arrivals) for each node and priority.
 	    • Average queue length (time-average) for each node and priority.

 	• Reference Flow Statistics:
 	    • End-to-end average delay.
 	    • Overall blocking ratio for the reference flow.