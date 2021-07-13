CC       = g++
CFLAGS   = -std=c++14  -O3 -Wall -Winline -Wshadow
RM       = rm -f


.PHONY: all clean

all: lbm

clean:
	$(RM) lbm 

lbm: lbm.cpp 
	$(CC) $(CFLAGS) -o lbm lbm.cpp
