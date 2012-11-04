all:
	g++ -c matrix.cpp -o matrix.o -02
	g++ matrix.o -g -o matrix_ss
	#if parallel only
	#g++ -MMD -MP  -DCACHE_LINE_SIZE=64 -DPROMETHEUS_TLS_SUPPORT -O3 -march=core2  -std=c++0x   -I../include -c matrix.cpp -o matrix.o
	#g++ matrix.o -Werror -g -std=c99 -lstdc++  -L../lib  -lpthread -o matrix_ss -lrt

clean:
	rm -f *.o *.d matrix_ss 
