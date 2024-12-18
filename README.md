# Gaussian blur

## Team

- Maria Sfîrăială
- Gabriela Grosu
- Iulia Dumitru

## Description

The Gaussian blur feature is obtained by blurring (smoothing) an image using
a Gaussian function to reduce the noise level. It can be considered as a nonuniform
low-pass filter that preserves low spatial frequency and reduces image noise
and negligible details in an image.  It is typically achieved by convolving an
image with a Gaussian kernel.

In this project we started from a
[serial implementation for Gaussian blur](https://github.com/anishmrao/AA_2018/blob/master/gaussian_blur.cpp)
that we paralelized using OpenMP, MPI and a hybrid approach with both OpenMP and MPI.

### OpenMP

### MPI

The image is divided into overlapping slices. Each slice is sent to an MPI process
and processed individually. They are then combined back together to form the
blurred image.

This is accomplished by scattering (`MPI_Scatterv`) the original image among
the processes, providing information regarding the size of the slice (different
for each process) and the offset from which to get the slice from (also varying).
The result is then gathered (`MPI_Gatherv`) into the final (blurred) image.

The slices need to be overlapping because the blurring function (`gaussian_blur`)
needs information from the cells immediately outside the slice. Hence, the slice
is originally bigger than the final blurred slice. This is handled by calculating
`edge_rows` and resizing the slices' sizes and displacement accordingly, both
when scattering and when gathering.

### Pthreads

### Hybrid: MPI + OpenMP

