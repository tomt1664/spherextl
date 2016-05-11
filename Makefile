# compiler flags:
CFLAGS  = 

exe = spherextl

CC = g++

OBJ = main.o graph_beads.o

$(exe): $(OBJ)
	$(CC) $(CFLAGS) -o $(exe) $(OBJ) 

main.o : main.cpp
	$(CC) $(CFLAGS) -c $*.cpp
graph_beads.o : graph_beads.cpp
	$(CC) $(CFLAGS) -c $*.cpp

clean:
	$(RM) $(OBJ) $(exe)
