cflags= "-std=c++11"
includes= -I "/Users/jiahaozhang/Documents/homework/mse561/mse561assignmentthree/head"
objdir=../src
vpath=$(objdir)
obejects=main.o atom.o
cxx=g++
atom.o: atom.cpp ndarrays.h atom.h parameters.h
	$(cxx) atom.cpp -c $(cflags) $(includes) ndarrays.h atom.h parameters.h
main.o: main.cpp ndarrays.h atom.h parameters.h
	$(cxx) main.cpp -c $(cflags) $(includes) ndarrays.h atom.h parameters.h
three: $(objects)
	$(cxx) $(objects) -o three
