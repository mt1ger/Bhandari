CC=g++-7 
LDFLAGS=-pthread
CFLAGS=-c -Wall
EXEC=Sim

SRCS=Network.cpp\
	 Dijkstra.cpp\
	 Bhandari.cpp\
	 RoutingTable.cpp\
	 Topology.cpp\
	 main.cpp\
	 ModulationFormats.cpp\
	 RandomVariable.cpp\
	 ResourceAssignment.cpp\
	 TrafficGenerator.cpp\
	 Event.cpp\
	 EventQueue.cpp\


OBJS=$(SRCS:.cpp=.o)
CLEANFILES=*.o $(EXEC)

all: $(OBJS) $(EXEC)
$(EXEC): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $@

Network.o: Network.cpp Network.h EventQueue.h RoutingTable.h TrafficGenerator.h Event.h ResourceAssignment.h
	$(CC) $(CFLAGS) Network.cpp

Topology.o: Topology.cpp Topology.h Network.h
	$(CC) $(CFLAGS) Topology.cpp 

Dijkstra.o: Dijkstra.cpp Dijkstra.h Network.h Topology.h
	$(CC) $(CFLAGS) Dijkstra.cpp

Bhandari.o: Bhandari.cpp Bhandari.h Network.h Dijkstra.h Topology.h
	$(CC) $(CFLAGS) Bhandari.cpp

RoutingTable.o: RoutingTable.cpp RoutingTable.h Network.h Dijkstra.h
	$(CC) $(CFLAGS) RoutingTable.cpp

ResourceAssignment.o: ResourceAssignment.cpp ResourceAssignment.h Event.h EventQueue.h RoutingTable.h ModulationFormats.h
	$(CC) $(CFLAGS) ResourceAssignment.cpp

TrafficGenerator.o: TrafficGenerator.cpp TrafficGenerator.h Network.h EventQueue.h RandomVariable.h
	$(CC) $(CFLAGS) TrafficGenerator.cpp

Event.o: Event.cpp Event.h
	$(CC) $(CFLAGS) Event.cpp

EventQueue.o: EventQueue.cpp EventQueue.h Event.h
	$(CC) $(CFLAGS) EventQueue.cpp

RandomVariable.o: RandomVariable.cpp RandomVariable.h 
	$(CC) $(CFLAGS) RandomVariable.cpp

ModulationFormats.o: ModulationFormats.cpp ModulationFormats.h Event.h Network.h
	$(CC) $(CFLAGS) ModulationFormats.cpp

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

.PHONY: clean
clean:
	rm $(CLEANFILES)


