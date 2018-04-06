
#include <pthread.h>
#include <time.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string.h>
#include "Network.h"
#include "RoutingTable.h"


using namespace std;


int main (int argc, char *argv[]) {
	string Path;
	vector<int> pred;

	Network * network;
	Network net;
	network = &net;

	if (argc != 2) {
		cout << "Please input arguments in the following order: " << endl;	
		cout << "\tThe file for network topology" << endl; 
		cout << endl;
		exit (0);
	}
	strcpy (network->FileName, argv[1]);
	
	network->init ();

	return 1;
}
