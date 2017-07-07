
CC = g++
CFLAGS =  -O3

all: SVHapSimulator

main.o: main.cpp 
	$(CC) $(CFLAGS) -c main.cpp

vcf_parse.o: vcf_parse.cpp
	$(CC) $(CFLAGS) -c vcf_parse.cpp
	
file_io.o: file_io.cpp
	$(CC) $(CFLAGS) -c file_io.cpp

public_func.o: public_func.cpp
	$(CC) $(CFLAGS) -c public_func.cpp

SVsim: main.o vcf_parse.o file_io.o public_func.o
	$(CC) $(CFLAGS) -o SVsim main.o vcf_parse.o file_io.o public_func.o
