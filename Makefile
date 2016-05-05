
gldemo: gldemo.o
	clang++ gldemo.o -lGL -lGLU -lglut -o gldemo

gldemo.o: gmath.hh gldemo.cc
	clang++ -c -std=c++11 gldemo.cc -o gldemo.o

clean:
	rm -f gldemo.o gldemo

test: gldemo
	./gldemo
