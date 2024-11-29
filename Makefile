all: main.o matrix.o gaussAll.o functions.o
	g++ main.o matrix.o gaussAll.o functions.o

main.o: main.cpp matrix.h gaussAll.h functions.h
	g++ -Wall -std=c++20 -O3 -c main.cpp

matrix.o: matrix.cpp matrix.h
	g++ -Wall -std=c++20 -O3 -c matrix.cpp

gaussAll.o: gaussAll.cpp gaussAll.h
	g++ -Wall -std=c++20 -O3 -c gaussAll.cpp

functions.o: functions.cpp functions.h
	g++ -Wall -std=c++20 -O3 -c functions.cpp

clean:
	rm -f *.o *.out
