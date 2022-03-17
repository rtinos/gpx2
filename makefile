CC = g++
FLAGS1 = -Wall
FLAGS2 = -g -DREPORT

CPPFLAGS = $(FLAGS1) $(FLAGS2)

objects = file_man.o global.o gpx.o util_functions.o

main: main.cpp $(objects)
	$(CC) $(CPPFLAGS) main.cpp $(objects) -o gpx2
