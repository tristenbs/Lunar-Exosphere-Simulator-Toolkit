
CFLAGS  = -g -pthread -Wunused -Warray-bounds $(DEFINITIONS) 

all :  L4tools.o

L4tools.o : L4tools.c  L4vector.c L4matpac.c L4ranlib.h L4ranlib.c ../L4common.h
	cc -c L4tools.c $(CFLAGS)


