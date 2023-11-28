cuda: main.cu
	nvcc -arch=sm_61 main.cu -o cuda
diago: scalapack.cpp
	g++ scalapack.cpp -o diago -lscalapack-openmpi
clean:
	rm -rf cuda diago