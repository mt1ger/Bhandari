// #define DEBUG_get_predecessor_list
// #define DEBUG_generate_dijkstra_routing_table
// #define DEBUG_get_shortest_path

#include <iostream>
#include <vector>
#include <stack>
#include "RoutingTable.h"

using namespace std;

vector< vector <int> > RoutingTable::single_src_routing_table (int src) {
	stack<int> pathStack;
	vector<int> shortestPath;
	vector< vector<int> > singleSrcRoutingTable;

	for (int j = 0; j < network->NumofNodes; j++) {
		if (j == src) {
			pathStack.push (-1);	
		}
		else
			pathStack.push (j);
		int temp = j;
		while (temp != src) {
			pathStack.push (predecessors[src].at (temp));	
			temp = predecessors[src].at (temp);
		}

		while(!pathStack.empty())
		{
			shortestPath.push_back(pathStack.top());
			pathStack.pop();		
		}	

		singleSrcRoutingTable.push_back (shortestPath);
		shortestPath.clear ();
	}

	return singleSrcRoutingTable;
}

void RoutingTable::get_predecessor_list () {
	vector<int> hpredecessors;
	Topology topology (network);
	Dijkstra dijkstra (network);

	topology.read_topology ();
	dijkstra.ajacent_nodes (dijkstra.AjacentNodes, network->NodesWeight);
	
	for (int i = 0; i < network->NumofNodes; i++) {
		dijkstra.shortest_path (i, -1, hpredecessors, network->NodesWeight);

		predecessors.push_back (hpredecessors);
		hpredecessors.clear ();
	}
	
	#ifdef DEBUG_get_predecessor_list
	cout << endl;
	cout << "Predecessor List: " << endl;
	for (int i = 0; i < network->NumofNodes; i++) {
		for (int j = 0; j < network->NumofNodes; j++) {
			cout << predecessors[i][j] + 1 << ' ';
		}
		cout << endl;
	}
	#endif
}


void RoutingTable::generate_dijkstra_routing_table () {

	get_predecessor_list ();

	for (int i = 0; i < network->NumofNodes; i++) {
		network->DRoutingTable.push_back (single_src_routing_table (i));
	}
	
	#ifdef  DEBUG_generate_dijkstra_routing_table	
	cout << endl;
	for (int i = 0; i < network->NumofNodes; i++) {
		cout << endl;
		cout << "Start to print table " << i << endl;
		for (int j = 0; j < network->NumofNodes; j++) {
			for (int k = 0; k < network->DRoutingTable[i][j].size (); k++) {
				cout << ' ' << network->DRoutingTable[i][j][k] + 1;	
			}
			cout << endl;
		}	
	}
	#endif
}

void RoutingTable::generate_routing_table () {
	vector< vector< vector<int> > > YroutingTable;
	Bhandari bhandari (network);
	
	generate_dijkstra_routing_table ();

	for (int src = 0; src < network->NumofNodes; src++) {
		for (int dest = 0; dest < network->NumofNodes; dest++) {
			cout << "***** ";
			cout << "SRC is: " << src + 1 << " and DEST is: " << dest + 1 << endl;
			if (src == dest) {
				vector<int> T (1, -1);
				vector< vector<int> > Z;
				Z.push_back (T);
				YroutingTable.push_back (Z);
				continue;
			}
			YroutingTable.push_back (bhandari.eliminate_common_links (src, dest, 10));
		}
		network->routingTable.push_back (YroutingTable);
		YroutingTable.clear ();
	}
}


// vector<int> RoutingTable::get_shortest_path (int src, int dest) {
// 	vector<int> shortestPath;
//
// 	
// 	for (int i = 0; i < network->DRoutingTable[src][dest].size (); i++) {
// 		shortestPath.push_back (network->DRoutingTable[src][dest][i]);	
// 	}
//
// 	#ifdef DEBUG_get_shortest_path
// 	for (int i = 0; i < network->DRoutingTable[src][dest].size (); i++) {
// 		cout << shortestPath[i] + 1 << ' '; 	
// 	}
// 	cout << endl;
// 	#endif
// 	
// 	return shortestPath;
// }

