CC=gcc
CFLAGS=-I/usr/include/gsl
LDIR=/usr/lib

LIBS=-lm -lgsl -lgslcblas

Default:
	$(CC) -o main.exe main.c PBC.c LJ.c -lm -lgsl -lgslcblas


debug:
	$(CC) -Wall -g -c -lm
	$(CC) -Wall -g -o main-debug.exe main.c PBC.c LJ.c -lm -lgsl -lgslcblas

profile: 
	$(CC) -pg -no-pie -fno-builtin -c main.c PBC.c LJ.c -lm -fno-stack-protector
	$(CC) -pg -no-pie -fno-builtin -o main-gprof.exe main.c PBC.o LJ.o -lm -fno-stack-protector
	
#save the gmon.out data(which you get after running the .exe) to a text file
profileText:
	$ gprof main-gprof.exe gmon.out > analyze.txt
	
plot:
	$ xmgrace -free -batch plot.bfile -nosafe

#CLEAN
clean:
	rm -f *.o *.exe 
#*.exe *.txt *.out *.dump *.dat
