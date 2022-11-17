
CFLAGS  = -g -pthread -Wunused -Warray-bounds $(DEFINITIONS)

all : L4surface.o


L4surface.o :   ../L4common.h  L4surface.c  \
                  L4temperature.c L4ephemeris.c L4-lola.c
	cc -c L4surface.c $(CFLAGS) 



