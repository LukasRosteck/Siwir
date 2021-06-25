CC       = g++
CFLAGS   = -std=c++14 -Wall -Wextra -Winline -Wshadow -Werror -O3 -DNDEBUG -pedantic  -D_XOPEN_SOURCE=700 
RM       = rm -f


.PHONY: all clean

all: mgsolve

clean:
	$(RM) mgsolve

mgsolve: multigrid.cpp
	$(CC) $(CFLAGS) -o mgsolve $<