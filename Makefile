CXX = /usr/bin/c++

CFLAGS = -std=c++11 -Wall

OBJS = delaunay.o


delaunay: ${OBJS}
	${CXX} ${CFLAGS} -o $@ $^

%.o: %.cpp
	${CXX} ${CFLAGS} -c -o $@ $^

clean:
	rm -f delaunay *.o

PHONY: clean
