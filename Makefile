CC       = g++
CFLAGS   = -std=c++14  -O3 -Wall -Winline -Wshadow
RM       = rm -f


.PHONY: all clean

all: lbm cmdparser.o

clean:
	$(RM) lbm cmdparser.o

lbm: cmdparser.o lbm.cpp 
	$(CC) $(CFLAGS) -o lbm cmdparser.o lbm.cpp

cmdparser.o: cmdparser.cpp
	$(CC) $(CFLAGS) -c $<