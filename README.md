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

The Gaussian Blur algorithm was parallelized using OpenMP to reduce execution time. 
The three most time-consuming functions, `box_blur_h()`, `box_blur_t()`, and `copy_image()`, 
were parallelized using OpenMP directives. 

Parallelization is achieved by dividing the workload into equal chunks, 
calculated as the total size divided by the number of threads (NUM_THREADS). 
The scheduling is static, ensuring that each thread processes a fixed number of rows, 
columns, or pixels, providing a uniform workload distribution.
This approach leverages available cores to speed up image processing.

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

Similarly to the OpenMP implementation, we got a speedrun by improving the `box_blur_h()`, `box_blur_t()`, and `copy_image()` functions using the simple `for` parallelization method.
We divide the loop between threads based on individual indeces with a formula dependant on the thread id, the total number of threads and the number of loops:

```C
int start = args->thread_id * (double)args->h / args->num_threads;
int end = MIN((args->thread_id + 1) * (double)args->h / args->num_threads, args->h);C
```

As such, no other synchronization primitive is needed, because each thread has exclusive access on a distinct portion of the image.

### Hybrid: MPI + OpenMP

The hybrid MPI and OpenMP implementation combines message-passing and multithreading
to optimize the Gaussian Blur algorithm. 

The image is divided into overlapping slices, which are distributed to MPI processes 
using MPI_Scatterv. Each MPI process processes its slice in parallel using OpenMP, 
leveraging multiple threads to execute the most time-consuming functions 
(box_blur_h, box_blur_t, and copy_image). 
Within each process, the workload is divided into chunks among the threads, ensuring 
efficient utilization of computational resources. After processing, the results are
gathered back into the final blurred image using MPI_Gatherv. 

This hybrid approach maximizes parallelism by distributing work across both nodes and threads.


# Bonus

### Hybrid: MPI + PThreads

This hybrid version combines the use of MPI for distributing image slices across different
processes and pthreads for parallelizing the blur operation on each process. 

The Gaussian blur algorithm is applied in two stages: the first stage applies the blur 
horizontally using pthreads, and the second stage applies it vertically.
After the main process (ROOT) divides the image into slices and sends them to the worker
processes, each process performs the blur on its respective slice. 

The results are then gathered and combined using MPI. This hybrid approach enables 
efficient use of resources in a distributed environment, improving performance for large images.
