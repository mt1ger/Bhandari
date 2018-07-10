/**************************************************
 * Best-Fit  
 **************************************************/
#include <iostream>
#include <string>
#include <climits>
#include <cfloat>
#include <list>
#include <cmath>
 
#include "ModulationFormats.h"
#include "ResourceAssignment.h"


void merge (vector<int> * ar, int low, int mid, int high)
{
	int n, m, i, j, k;
	n = mid - low + 1;
	m = high - mid;
	vector<int> R, L;

	for (i = 0; i < n; i++)
	{
		L.push_back (ar->at (low + i));
	}
	for (j = 0; j < m; j++)
	{
		R.push_back (ar->at (mid + j + 1));
	}
	for (i = 0, j = 0, k = low; k <= high; k++) //
	{
		if (i < n && j < m)
		{
			if (L[i] <= R[j])
			{
				ar->at (k) = L[i];
				i++;
			}
			else
			{
				ar->at (k) = R[j];
				j++;
			}
		}
		else if (i < n && j >= m)
		{
			ar->at (k) = L[i];
			i++;
		}
		else if (i >= n && j < m)
		{
			ar->at (k) = R[j];
			j++;
		}
	}
}


void merge_sort(vector<int> * ar, int low, int high)
{
	int mid;
	if (low < high)
	{
		mid = ((high + low) / 2);
		merge_sort (ar, low, mid);
		merge_sort (ar, mid+1, high);
		merge (ar, low, mid, high);
	}
}

double logrithm (double x, double base) {
	return log (x) / log (base);
}

void ResourceAssignment::check_availability_source (unsigned int predecessor, unsigned int successor, CircuitRequest * circuitRequest) {
	vector<int> HAvailableSpecSlots;

	AvailableSpecSlots.clear ();
	for (int c = 0; c < network->NumofCores; c++) {
		for (int i = 0; i < NUMOFSPECTRALSLOTS; i++) {
			if (network->SpectralSlots[predecessor][successor][c][i] == false) {
				HAvailableSpecSlots.push_back (c);
				HAvailableSpecSlots.push_back (i);
				AvailableSpecSlots.push_back (HAvailableSpecSlots);
				HAvailableSpecSlots.clear ();	
			}
		}
	}
}


void ResourceAssignment::check_availability_link (vector<int> * CircuitRoute) {

	list< vector<int> >::iterator i;

	for (int r = 2; r < CircuitRoute->size (); r++) {
		for (i = AvailableSpecSlots.begin (); i != AvailableSpecSlots.end (); i++) {
			if (network->SpectralSlots[CircuitRoute->at (r - 1)][CircuitRoute->at (r)][i->at (0)][i->at (1)] == true) {
				i = AvailableSpecSlots.erase (i);
				i--;
				// AvailableSpecSlots.erase (i); // This line will only work on Mac, will cause Seg Fault (dump core) on ubuntu.
			}
		}
	}

	#ifdef DEBUG_print_AvailableSpecSlots
	cout << "PRINT Available Spectral Slots:" << endl;
	for (i = AvailableSpecSlots.begin (); i != AvailableSpecSlots.end (); i++) {
		for (int j = 0; j < i->size (); j++) {
			cout << i->at (j) << ' ';
		}
		cout << "    ";
	}
	cout << endl;
	#endif
}


double ResourceAssignment::path_compactness (list< vector<int> > &PotentialSec) {
	int TotalOSS = 0; // OSS = Occupied  Spec Slots
	int TotalASS = 0; // ASS = Available Spec Slots
	double PC = 0; // PC = Path Compactness
	vector<int> HPotentialSections;
	list< vector<int> >::iterator i;

	/* Use the Available Spectral Slots to form Potential Sections */
	for (i = AvailableSpecSlots.begin (); i != AvailableSpecSlots.end (); i++) {
		if (HPotentialSections.empty ()) {
			if (i != (--AvailableSpecSlots.end ())) {
				HPotentialSections.push_back (i->at (0));
				HPotentialSections.push_back (i->at (1));
				HPotentialSections.push_back (i->at (1));
				PotentialSec.push_back (HPotentialSections);
			}
		}
		else {
			if ((i->at (1) == (HPotentialSections.at (2) + 1)) && (i->at (0) == HPotentialSections.at (0))) {
				HPotentialSections.at (2) = i->at(1);
				*(--PotentialSec.end ()) = HPotentialSections;
			}
			else {
				if (HPotentialSections.at (2) != HPotentialSections.at (1)) {
					HPotentialSections.clear ();
					if (i != (--AvailableSpecSlots.end ())) {
						HPotentialSections.push_back (i->at (0));
						HPotentialSections.push_back (i->at (1));
						HPotentialSections.push_back (i->at (1));
						PotentialSec.push_back (HPotentialSections);
					}
				}
				else {
					HPotentialSections.clear ();
					PotentialSec.erase (--PotentialSec.end ());
					if (i != (--AvailableSpecSlots.end ())) {
						HPotentialSections.push_back (i->at (0));
						HPotentialSections.push_back (i->at (1));
						HPotentialSections.push_back (i->at (1));
						PotentialSec.push_back (HPotentialSections);
					}
				}
			}
		}
	}

	for (i = PotentialSec.begin (); i != PotentialSec.end (); i++) {
		TotalASS += i->at (2) - i->at (1) + 1;
	}
	// TotalOSS = network->NumofCores * NUMOFSPECTRALSLOTS - TotalASS;

	if (TotalOSS == 0) PC = PC_MAX; 
	else if (TotalOSS == NUMOFSPECTRALSLOTS) PC = 0;
	else { 
		// PC = ((double) network->NumofCores * NUMOFSPECTRALSLOTS / TotalOSS) * ((double)TotalASS / PotentialSec.size ());
		PC = ((double) TotalASS / PotentialSec.size ());
	}

	#ifdef DEBUG_print_PotentialSections_and_PC
	cout << "PRINT Potential Sections: " << endl;
	for (i = PotentialSec.begin (); i != PotentialSec.end (); i++) {
		for (int j = 0; j < i->size (); j++) {
			cout << i->at (j) << ' ';
		}
		cout << "    ";
	}
	cout << endl;
	cout << "Num of Core is: " << network->NumofCores << "    Num of TotalOSS is: " << TotalOSS << "    TotalASS is: " << TotalASS << "    Num of AS is: " << PotentialSec.size () << endl; 
	cout << "PC is: " << PC << endl;
	#endif

	return PC;
}


void ResourceAssignment::handle_requests (CircuitRequest * circuitRequest) {
	RoutingTable routingTable (network);	
	ModulationFormats modulationFormats (circuitRequest, network);

	vector< vector<int> > PathList;
	vector<int> CircuitRoute;
	bool AvailableFlag = true;
	vector< vector<int> > AssignedSpectralSection;
	vector<int> HAssignedSpectralSection;
	unsigned int TempNumofTransponders = 0;
	string MF = "BPSK";
	unsigned int mfTimes = 0, NumofOccupiedSpectralSlots = 0;
	vector< vector<int> > mfTimesList, NoOSSList; 
	list< vector<double> > WoPList;
	vector<int> HmfTimesList, HNoOSSList;
	vector<double> HWoPList;
	vector<int> PathSortedIndex;
	vector<int> AddrTranslation;
	list<double> OnlyWeight; // Only Weight stored in the vector without indication of which path it belongs
	vector<string> MFList;

	PotentialSections.clear ();
	PathList = network->routingTable[circuitRequest->Src][circuitRequest->Dest];

	#ifdef DISPLAY_Available_Path
	cout << "The Available " << PathList.size () << " Paths are:" << endl;
	for (int i = 0; i < PathList.size (); i++) {
		for (int j = 0; j < PathList[i].size () - 1; j++) {
			cout << PathList[i][j] << "-->";
		}
		cout << PathList[i][PathList[i].size () - 1] << endl;
	}
	#endif

	for (int i = 0; i < PathList.size (); i++) {
		int NoOSS;
		double WoP;
		double PC;
		bool AvFlag = true;
		list< vector<int> > PSections;

		#ifdef DEBUG_path_checking
		cout << endl << "*** BEGIN: Path checking ***" << endl;
		cout << "Path under checking is: " << endl;
		for (int j = 0; j < PathList[i].size () - 1; j++) {
			cout << PathList[i][j] << "->"; 
		}
		cout << PathList[i][PathList[i].size () - 1] << endl;
		#endif

		NoOSS = modulationFormats.mf_chosen (PathList[i], &circuitRequest->OccupiedSpectralSlots, &circuitRequest->DataSize, &MF, &mfTimes);
		MFList.push_back (MF);

		#ifdef DEBUG_path_checking
		cout << "Bits Per Signal is: " << mfTimes << " " << MF << endl;
		#endif

		// Calculate possible SpectralSlotSections on the link between every two nodes
		check_availability_source (PathList[i][0], PathList[i][1], circuitRequest);
		check_availability_link (&PathList[i]);

		if (NoOSS <= AvailableSpecSlots.size ()) 
		{

			PC = path_compactness (PSections);
			list< vector<int> >::iterator Pointer;
			int AvailSize = 0;
			for (Pointer = PSections.begin (); Pointer != PSections.end (); Pointer++) {
				AvailSize += Pointer->at (2) - Pointer->at (1) + 1; 
			}

			if (NoOSS > AvailSize - PSections.size ()) {
				AvFlag = false;
			}


			int NoH = PathList[i].size () - 1;

			#ifdef DEBUG_path_checking
			cout << "Bits Per Signal is: " << mfTimes << "	Num of Hops is: " << NoH << endl;
			#endif

			if (AvFlag == true)
			{
				AddrTranslation.push_back (i);
				PotentialSections.push_back (PSections);
				WoP = network->a * ((double)(mfTimes - LmfTimes) / (HmfTimes - LmfTimes)) + network->b * ((double)1 / NoH - (double)1 / network->MaxNoH) / ((double)1 / network->MinNoH - (double)1 / network->MaxNoH) + network->c * (PC - network->NumofCores * PC_MIN) / (network->NumofCores * (PC_MAX - PC_MIN)); 

				#ifdef DEBUG_path_checking
				cout << "Weight of BpS is: " << ((double) (mfTimes - LmfTimes) / (HmfTimes - LmfTimes)) << "    Weights of Hops is: " << ((double)1 / NoH - (double)1 / network->MaxNoH) / ((double)1 / network->MinNoH - (double)1 /network->MaxNoH) << "    Weight of PC is: " << (double)(PC - PC_MIN) / (network->NumofCores * (PC_MAX - PC_MIN)) << endl;
				cout << "WoP is: " << WoP << endl;
				#endif

				HmfTimesList.push_back (i);
				HmfTimesList.push_back (mfTimes);
				mfTimesList.push_back (HmfTimesList);
				HmfTimesList.clear ();

				HWoPList.push_back (i);
				HWoPList.push_back (WoP);
				OnlyWeight.push_back (WoP);
				WoPList.push_back (HWoPList);
				HWoPList.clear ();
			}
		}

		#ifdef DEBUG_path_checking
		cout << "*** END: Path checking ***" << endl;
		#endif
	}

	// Sort WoPList
	OnlyWeight.sort ();
	list<double>::iterator IOW;
	list< vector<double> >::iterator IWL;
	for (IOW = OnlyWeight.begin (); IOW != OnlyWeight.end (); IOW++) {
		for (IWL = WoPList.begin (); IWL != WoPList.end (); IWL++) {
			if (IWL->at (1) == *IOW) {
				PathSortedIndex.push_back (IWL->at (0));
				IOW = OnlyWeight.erase (IOW);
				IOW--;
				WoPList.erase (IWL);
				IWL--;
				break;
			}
		}
	}

	#ifdef DISPLAY_path_order
	cout << "Sorted Path Index are:" << endl;
	for (auto i : PathSortedIndex) {
		cout << i << ' ';
	}
	cout << endl;
	#endif

	if (PathSortedIndex.size () == 0) {
		AvailableFlag = false;
	}
	
	// Put the selected PathList[i] into CircuitRoute
	for (int i = PathSortedIndex.size () - 1; i >= 0 ; i--) {
		CircuitRoute = PathList[PathSortedIndex[i]];

		#ifdef DISPLAY_selected_path
		cout << "Trying to Allocate request on Path: " << endl;
		for (int i = 0; i < CircuitRoute.size () - 1; i++) {
			cout << CircuitRoute[i] << "-->";
		}
		cout << CircuitRoute[CircuitRoute.size () - 1] << endl;
		#endif

		int SizeSum = 0;
		int temp = 0;
		vector<int> Size;
		int SizeCnt = 0;
		list< vector<int> > SortedSections;
		list< vector<int> >::iterator IPS;

		int addr;
		MF = MFList[PathSortedIndex[i]];
		for (addr = 0; addr < AddrTranslation.size (); addr++) {
			if (PathSortedIndex[i] == AddrTranslation[addr]) break;
		}

		for (IPS = PotentialSections[addr].begin (); IPS != PotentialSections[addr].end (); IPS++) {
			Size.push_back (IPS->at (2) - IPS->at (1) + 1); 
			SizeCnt++;
		}
		merge_sort (&Size, 0, SizeCnt - 1);
		for (int s = SizeCnt - 1; s >= 0; s--) {
			for (IPS = PotentialSections[addr].begin (); IPS != PotentialSections[addr].end (); IPS++) {
				if ((IPS->at (2) - IPS->at (1) + 1) == Size[s]) {
					SortedSections.push_back (*IPS);
					PotentialSections[addr].erase (IPS);
					break;
				}
			}
		}

		#ifdef DEBUG_print_SortedSections
		cout << "PRINT Sorted Sections:" << endl;
		for (IPS = SortedSections.begin (); IPS != SortedSections.end (); IPS++) {
			for (int j = 0; j < IPS->size (); j++) {
				cout << IPS->at (j) << ' ';
			}
			cout << "    ";
		}
		cout << endl;
		#endif 

		#ifdef DEBUG_print_resource_state_on_the_path
		cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
		cout << "PRINT resources BEFORE Allcation: " << endl;
		for (int i = 1; i < CircuitRoute.size (); i++) {
			cout << "On link " << CircuitRoute[i - 1] << " to " << CircuitRoute[i] << endl;
			for (int c = 0; c < network->NumofCores; c++) {
				cout << "On Core " << c << endl;
				for (int j = 0; j < NUMOFSPECTRALSLOTS; j++) {
					cout <<  network->SpectralSlots[CircuitRoute[i - 1]][CircuitRoute[i]][c][j] << ' ';
				}
				cout << endl;
			}
		}
		#endif

		// Allocation 
		bool BiggerFlag = false;
		bool SmallerFlag = false;
		int OSS;  
		for (int j = 0; j < mfTimesList.size (); j++) {
			if (PathSortedIndex[i] == mfTimesList[j][0]) {
				mfTimes = mfTimesList[j][1];
				OSS = ceil ((double) circuitRequest->OccupiedSpectralSlots / mfTimes);
			}
		}
		
		list< vector<int> >::iterator ISS = SortedSections.begin ();
		list< vector<int> >::iterator IT;
		while (OSS > 0) {
			/* Allocation */
			if (ISS == SortedSections.end ()) {
				if (BiggerFlag == true) {
					ISS--;
					SmallerFlag = true;
					continue;
				}
			}
			else if (ISS->at (2) - ISS->at (1) + 1 == OSS + GB) {
				HAssignedSpectralSection.push_back (ISS->at (0));
				HAssignedSpectralSection.push_back (ISS->at (1));
				HAssignedSpectralSection.push_back (ISS->at (2));
				AssignedSpectralSection.push_back (HAssignedSpectralSection);
				TempNumofTransponders++;
				break;
			}
			else if (ISS->at (2) - ISS->at (1) + 1 > OSS + GB) {
				BiggerFlag = true;
				if (SmallerFlag == true) {
					HAssignedSpectralSection.push_back (ISS->at (0));
					HAssignedSpectralSection.push_back (ISS->at (1));
					HAssignedSpectralSection.push_back (ISS->at (1) + OSS + GB - 1); 
					AssignedSpectralSection.push_back (HAssignedSpectralSection);
					TempNumofTransponders++;
					break;
				}
				ISS++;
			}
			else if (ISS->at (2) - ISS->at (1) + 1 < OSS + GB) {
				if (BiggerFlag == true) {
					ISS--;
					SmallerFlag = true;
					continue;
				}
				HAssignedSpectralSection.push_back (ISS->at (0));
				HAssignedSpectralSection.push_back (ISS->at (1));
				HAssignedSpectralSection.push_back (ISS->at (2));
				AssignedSpectralSection.push_back (HAssignedSpectralSection);
				OSS -= ISS->at (2) - ISS->at (1) + 1 - GB;
				HAssignedSpectralSection.clear ();
				ISS = SortedSections.erase (ISS);
				TempNumofTransponders++;
			}
		}
		if (AvailableFlag == true) break;
	}
	
	if (AssignedSpectralSection.size () > network->SectionNumLimitation) {
		AvailableFlag = false;
	}

	
	if (AvailableFlag == false) {
		network->NumofDoneRequests++;

		#ifdef DEBUG_collect_EventID_of_blocked_requests
		network->BlockedRequests.push_back (circuitRequest->EventID);
		#endif

		#ifdef PRINT_allocation_block_release
		cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
		cout << "Request " << circuitRequest->EventID << " is blocked" << endl;
		cout << "Modulation Format: " << MF << endl;
		cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
		#endif

			network->NumofFailedRequests++;
	}
	else if (AvailableFlag == true) {
		int temp = AssignedSpectralSection[0][0];
		int CoreCnter = 1;
		for (int i = 0; i < AssignedSpectralSection.size (); i++) {
			for (int p = 1; p < CircuitRoute.size (); p++) {
				for (int j = AssignedSpectralSection[i].at (1); j <= AssignedSpectralSection[i].at (2); j++) {
					network->SpectralSlots[CircuitRoute[p -1]][CircuitRoute[p]][AssignedSpectralSection[i].at (0)][j] = true;
			}
			}
			if (AssignedSpectralSection[i][0] != temp) {
				temp = AssignedSpectralSection[i][0];
				CoreCnter++;
			}
		}
		for (int i = 0; i < AssignedSpectralSection.size (); i++) {
			NumofOccupiedSpectralSlots += AssignedSpectralSection[i][2] - AssignedSpectralSection[i][1] + 1; 
		}

		#ifdef PRINT_allocation_block_release
		cout << "------------------------------------------------------------" << endl;
		cout << "Request ID: " << circuitRequest->EventID << "\tStart: " << circuitRequest->EventTime << "\tEnd: " << circuitRequest->StartTime + circuitRequest->Duration << endl;
		cout << "Source: " << circuitRequest->Src << "  Destination: " << circuitRequest->Dest << "  ModulationFormats: " << MF << endl;
		cout << "Path: ";
		for(unsigned int t = 0; t < CircuitRoute.size() - 1; t++)
			cout << CircuitRoute.at(t) << " --> ";
		cout << CircuitRoute.at (CircuitRoute.size() - 1) << endl;

		for (int i = 0; i < AssignedSpectralSection.size (); i++) {
			cout << "Core: " << AssignedSpectralSection[i][0] << "  Spectral Section: " << AssignedSpectralSection[i][1] << " to " << AssignedSpectralSection[i][2] << endl; 
		}

		cout << "# of Transponders Used: " << TempNumofTransponders << endl;
		cout << "# of Core Used: " << CoreCnter << endl;
		cout << "------------------------------------------------------------" << endl;
		#endif

		CircuitRelease * circuitRelease;
		circuitRelease = new CircuitRelease (circuitRequest->EventID, CircuitRoute, AssignedSpectralSection, circuitRequest->StartTime + circuitRequest->Duration, TempNumofTransponders);
		eventQueue->queue_insert (circuitRelease);

		network->NumofAllocatedRequests++;
		network->NumofSections = TempNumofTransponders;
		// network->TotalHoldingTime += circuitRequest->Duration;
		network->NumofTransponders = network->NumofTransponders + TempNumofTransponders;
		network->TotalTranspondersUsed += TempNumofTransponders;
		network->TotalCoresUsed += CoreCnter;
		network->TotalGBUsed += TempNumofTransponders;
		network->TotalDataSize += circuitRequest->DataSize;
		network->TotalSSUsed += NumofOccupiedSpectralSlots * mfTimes;
		network->TotalSSOccupied += (NumofOccupiedSpectralSlots + TempNumofTransponders) * mfTimes;
		
		#ifdef DISPLAY_metrics
		cout << "*** METRICS ***" << endl;
		cout << "    Bits Per Signal: " << mfTimes << endl;
		cout << "    NumofSections: " << network->NumofSections << endl; 
		cout << "    NumofTransponders: " << network->NumofTransponders << "  This time: " << TempNumofTransponders << endl; 
		cout << "    TotalTranspondersUsed: " << network->TotalTranspondersUsed << "  This time: " << TempNumofTransponders << endl; 
		cout << "    TotalGBUsed: " << network->TotalGBUsed << "  This time: " << TempNumofTransponders << endl; 
		cout << "    TotalDataSize: " << network->TotalDataSize << "  This time: " << circuitRequest->DataSize << endl; 
		cout << "    TotalSSUsed: " << network->TotalSSUsed << "  This time: " << NumofOccupiedSpectralSlots * mfTimes << endl; 
		cout << "    TotalSSOccupied: " << network->TotalSSOccupied << "  This time: " << (NumofOccupiedSpectralSlots + TempNumofTransponders) * mfTimes << endl; 
		#endif
	}

	#ifdef DEBUG_print_resource_state_on_the_path
	cout << "PRINT resources AFTER Allocation:" << endl;
	for (int i = 1; i < CircuitRoute.size (); i++) {
		cout << "On link " << CircuitRoute[i - 1] << " to " << CircuitRoute[i] << endl;
		for (int c = 0; c < network->NumofCores; c++) {
			cout << "On Core " << c << endl;
			for (int j = 0; j < NUMOFSPECTRALSLOTS; j++) {
				cout <<  network->SpectralSlots[CircuitRoute[i - 1]][CircuitRoute[i]][c][j] << ' ';
			}
			cout << endl;
		}
	}
	cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	cout << endl << endl;
	#endif
}


void ResourceAssignment::handle_releases (CircuitRelease * circuitRelease) {
	#ifdef DEBUG_print_resource_state_on_the_path
	cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	cout << "PRINT resources BEFORE Release:" << endl;
	for (int i = 1; i < circuitRelease->CircuitRoute.size (); i++) {
		cout << "On link " << circuitRelease->CircuitRoute[i-1] << " to " << circuitRelease->CircuitRoute[i] << endl;
		for (int c = 0; c < network->NumofCores; c++) {
			cout << "On Core " << c << endl;
			for (int j = 0; j < NUMOFSPECTRALSLOTS; j++) {
				cout <<  network->SpectralSlots[circuitRelease->CircuitRoute[i - 1]][circuitRelease->CircuitRoute[i]][c][j] << ' ';
			}
			cout << endl;
		}
	}
	#endif

	for (int i = 1; i < circuitRelease->CircuitRoute.size (); i++) {
		for (int j = 0; j < circuitRelease->OccupiedSpectralSection.size (); j++) {
			for (int k = circuitRelease->OccupiedSpectralSection[j][1]; k <= circuitRelease->OccupiedSpectralSection[j][2]; k++) {
				network->SpectralSlots[circuitRelease->CircuitRoute[i - 1]][circuitRelease->CircuitRoute[i]][circuitRelease->OccupiedSpectralSection[j][0]][k] = false;	
			}
		}
	}

	network->NumofDoneRequests++;
	network->NumofTransponders = network->NumofTransponders - circuitRelease->TranspondersUsed;

	#ifdef PRINT_allocation_block_release
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "Released Event on Path:" << endl;
	for (int i = 0; i < circuitRelease->CircuitRoute.size () - 1; i++) {
		cout << circuitRelease->CircuitRoute[i] << "-->";
	}
	cout << circuitRelease->CircuitRoute[circuitRelease->CircuitRoute.size () - 1] << endl;
	cout << "Release Event: " << circuitRelease->EventID << "\tTime: " << circuitRelease->EventTime << endl;
	for (int i = 0; i < circuitRelease->OccupiedSpectralSection.size (); i++) {
		cout << "Core: " << circuitRelease->OccupiedSpectralSection[i][0] << "  Spectral Section: " << circuitRelease->OccupiedSpectralSection[i][1] << " to " << circuitRelease->OccupiedSpectralSection[i][2] << endl;
	}
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	#endif

	#ifdef DEBUG_print_resource_state_on_the_path
	cout << "PRINT resources AFTER Release:" << endl;
	for (int i = 1; i < circuitRelease->CircuitRoute.size (); i++) {
		cout << "On link " << circuitRelease->CircuitRoute[i-1] << " to " << circuitRelease->CircuitRoute[i] << endl;
		for (int c = 0; c < network->NumofCores; c++) {
			cout << "On Core " << c << endl;
			for (int j = 0; j < NUMOFSPECTRALSLOTS; j++) {
				cout <<  network->SpectralSlots[circuitRelease->CircuitRoute[i - 1]][circuitRelease->CircuitRoute[i]][c][j] << ' ';
			}
			cout << endl;
		}
	}
	cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
	cout << endl << endl;
	#endif
}

