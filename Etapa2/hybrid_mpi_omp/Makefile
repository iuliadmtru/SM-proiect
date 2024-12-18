CPP=mpic++
CPPFLAGS=-O3 -g -fno-omit-frame-pointer -Wall -fopenmp -I include
LDLIBS=-ljpeg

build: gaussian_blur

gaussian_blur: src/gaussian_blur.cpp
	$(CPP) $^ -o gaussian_blur $(CPPFLAGS) $(LDLIBS)

clean:
	rm -rf gaussian_blur time.log

.PHONY: clean
