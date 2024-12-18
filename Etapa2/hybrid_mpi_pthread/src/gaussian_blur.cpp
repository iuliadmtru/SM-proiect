#include "image.hpp"
#include <mpi.h>

#define ROOT 0

int NUM_THREADS;
#define MIN(a,b) (((a)<(b))?(a):(b))

// Utility function to display the image
void display(Image *source)
{
	int size = source->height * source->width * source->components;
	cout << "Height : " << source->height << endl;
	cout << "Width : " << source->width << endl;
	cout << "Components : " << source->components << endl;
	cout << "Size : " << size << endl;
	for (int i = 0; i < size; ++i) {
		cout << (int)source->image[i] << " ";
		if ((i + 1) % source->components == 0)
			cout << endl;
	}

	cout << endl << "---------------------------------------------" << endl;
}



void set_pixel(Image *source, int index, int r, int g, int b)
{
	if (r < 0)
		r = 0;
	if (r > 255)
		r = 255;

	if (g < 0)
		g = 0;
	if (g > 255)
		g = 255;

	if (b < 0)
		b = 0;
	if (b > 255)
		b = 255;

	int index_rgb = index * source->components;
	source->image[index_rgb] = r;
	source->image[index_rgb + 1] = g;
	source->image[index_rgb + 2] = b;
}

int *boxes_for_gauss(double sigma, int n)
{
	double w_ideal = sqrt((12 * sigma * sigma / n) + 1);
	int wl = (int)floor(w_ideal);
	if (wl % 2 == 0)
		wl--;

	int wu = wl + 2;

	double m_ideal =
	    (12 * sigma * sigma - n * wl * wl - 4 * n * wl - 3 * n) /
	    (-4 * wl - 4);
	int m = round(m_ideal);

	int *sizes = new int[n]();
	for (int i = 0; i < n; i++) {
		sizes[i] = i < m ? wl : wu;
	}
	return sizes;
}

struct args_box_blur {
	double iarr;

	Image *source;
	Image *target;
	int w, h, radius;
	int num_threads;
	int thread_id;
};

void *box_blur_h_tfunc(void *arg) {
	args_box_blur *args = (args_box_blur *)arg;

	int start = args->thread_id * (double)args->h / args->num_threads;
    int end = MIN((args->thread_id + 1) * (double)args->h / args->num_threads, args->h);

	for (int i = start; i < end; i++) {
		int ti = i * args->w;
		int li = ti;
		int ri = ti + args->radius;
		int fv = ti * args->source->components;
		int lv = (ti + args->w - 1) * args->source->components;

		unsigned currennt_r = args->source->image[fv] * (args->radius + 1);
		unsigned currennt_g = args->source->image[fv + 1] * (args->radius + 1);
		unsigned currennt_b = args->source->image[fv + 2] * (args->radius + 1);

		int pixel;
		for (int j = 0; j < args->radius; j++) {
			pixel = (ti + j) * args->source->components;
			currennt_r += args->source->image[pixel];
			currennt_g += args->source->image[pixel + 1];
			currennt_b += args->source->image[pixel + 2];
		}

		for (int j = 0; j <= args->radius; j++) {
			pixel = ri * args->source->components;
			ri++;
			currennt_r += (args->source->image[pixel] - args->source->image[fv]);
			currennt_g += (args->source->image[pixel + 1] - args->source->image[fv + 1]);
			currennt_b += (args->source->image[pixel + 2] - args->source->image[fv + 2]);

			set_pixel(args->target, ti++, currennt_r * args->iarr,
				  currennt_g * args->iarr, currennt_b * args->iarr);
		}

		for (int j = args->radius + 1; j < args->w - args->radius; j++) {
			int first_pixel = ri * args->source->components;
			int second_pixel = li * args->source->components;
			ri++;
			li++;

			currennt_r += (args->source->image[first_pixel] -
				       args->source->image[second_pixel]);
			currennt_g += (args->source->image[first_pixel + 1] -
				       args->source->image[second_pixel + 1]);
			currennt_b += (args->source->image[first_pixel + 2] -
				       args->source->image[second_pixel + 2]);

			set_pixel(args->target, ti++, currennt_r * args->iarr,
				  currennt_g * args->iarr, currennt_b * args->iarr);
		}

		for (int j = args->w - args->radius; j < args->w; j++) {
			pixel = li * args->source->components;
			li++;

			currennt_r += (args->source->image[lv] - args->source->image[pixel]);
			currennt_g += (args->source->image[lv + 1] - args->source->image[pixel + 1]);
			currennt_b += (args->source->image[lv + 2] - args->source->image[pixel + 2]);

			set_pixel(args->target, ti++, currennt_r * args->iarr,
				  currennt_g * args->iarr, currennt_b * args->iarr);
		}
	}

	return NULL;
}

void box_blur_h(Image *source, Image *target, int w, int h, int radius)
{
	double iarr = (double)1 / (radius + radius + 1);

    pthread_t threads[NUM_THREADS];
    int rc;
    void *status;
    args_box_blur arg[NUM_THREADS];

    for (int i = 0; i < NUM_THREADS; ++i) {
        arg[i].iarr = iarr;
        arg[i].source = source;
        arg[i].target = target;
        arg[i].w = w;
		arg[i].h = h;
		arg[i].radius = radius;
		arg[i].num_threads = NUM_THREADS;
		arg[i].thread_id = i;
    }

    for (int i = 0; i < NUM_THREADS; ++i) {
        rc = pthread_create(&threads[i], NULL, box_blur_h_tfunc, &arg[i]);

		if (rc) {
			printf("Couldn't create thread %d\n", i);
			exit(-1);
		}
	}

	for (int i = 0; i < NUM_THREADS; ++i) {
		rc = pthread_join(threads[i], &status);

		if (rc) {
			printf("Couldn't join thread %d\n", i);
			exit(-1);
		}
	}
}

void *box_blur_t_tfunc(void *arg)
{
	args_box_blur *args = (args_box_blur *)arg;

	int start = args->thread_id * (double)args->w / args->num_threads;
    int end = MIN((args->thread_id + 1) * (double)args->w / args->num_threads, args->w);

	for (int i = start; i < end; i++) {
		int ti = i;
		int li = ti;
		int ri = ti + args->radius * args->w;

		int fv = ti * args->source->components;
		int lv = (ti + args->w * (args->h - 1)) * args->source->components;

		unsigned currennt_r = args->source->image[fv] * (args->radius + 1);
		unsigned currennt_g = args->source->image[fv + 1] * (args->radius + 1);
		unsigned currennt_b = args->source->image[fv + 2] * (args->radius + 1);

		int pixel;
		
		for (int j = 0; j < args->radius; j++) {
			pixel = (ti + j * args->w) * args->source->components;
			currennt_r += args->source->image[pixel];
			currennt_g += args->source->image[pixel + 1];
			currennt_b += args->source->image[pixel + 2];
		}

		for (int j = 0; j <= args->radius; j++) {
			pixel = ri * args->source->components;
			currennt_r += (args->source->image[pixel] - args->source->image[fv]);
			currennt_g += (args->source->image[pixel + 1] - args->source->image[fv + 1]);
			currennt_b += (args->source->image[pixel + 2] - args->source->image[fv + 2]);

			set_pixel(args->target, ti, currennt_r * args->iarr,
				  currennt_g * args->iarr, currennt_b * args->iarr);

			ri += args->w;
			ti += args->w;
		}

		for (int j = args->radius + 1; j < args->h - args->radius; j++) {
			int first_pixel = ri * args->source->components;
			int second_pixel = li * args->source->components;

			currennt_r += (args->source->image[first_pixel] -
				       args->source->image[second_pixel]);
			currennt_g += (args->source->image[first_pixel + 1] -
				       args->source->image[second_pixel + 1]);
			currennt_b += (args->source->image[first_pixel + 2] -
				       args->source->image[second_pixel + 2]);

			set_pixel(args->target, ti, currennt_r * args->iarr,
				  currennt_g * args->iarr, currennt_b * args->iarr);

			li += args->w;
			ri += args->w;
			ti += args->w;
		}

		for (int j = args->h - args->radius; j < args->h; j++) {
			int pixel = li * args->source->components;

			currennt_r += (args->source->image[lv] - args->source->image[pixel]);
			currennt_g += (args->source->image[lv + 1] - args->source->image[pixel + 1]);
			currennt_b += (args->source->image[lv + 2] - args->source->image[pixel + 2]);

			set_pixel(args->target, ti, currennt_r * args->iarr,
				  currennt_g * args->iarr, currennt_b * args->iarr);

			li += args->w;
			ti += args->w;
		}
	}

	return NULL;
}


void box_blur_t(Image *source, Image *target, int w, int h, int radius)
{
	double iarr = (double)1 / (radius + radius + 1);
    pthread_t threads[NUM_THREADS];
    int rc;
    void *status;
    args_box_blur arg[NUM_THREADS];

    for (int i = 0; i < NUM_THREADS; ++i) {
        arg[i].iarr = iarr;
        arg[i].source = source;
        arg[i].target = target;
        arg[i].w = w;
		arg[i].h = h;
		arg[i].radius = radius;
		arg[i].num_threads = NUM_THREADS;
		arg[i].thread_id = i;
    }

    for (int i = 0; i < NUM_THREADS; ++i) {
        rc = pthread_create(&threads[i], NULL, box_blur_t_tfunc, &arg[i]);

		if (rc) {
			printf("Couldn't create thread %d\n", i);
			exit(-1);
		}
	}

	for (int i = 0; i < NUM_THREADS; ++i) {
		rc = pthread_join(threads[i], &status);

		if (rc) {
			printf("Couldn't join thread %d\n", i);
			exit(-1);
		}
	}
}

struct args_copy_image {
	Image *source;
	Image *target;
	int thread_id, num_threads;
};

void *copy_image_tfunc(void *arg) {
	args_copy_image *args = (args_copy_image *)arg;
	int size = args->source->height * args->source->width * args->source->components;

	int start = args->thread_id * (double)size / args->num_threads;
    int end = MIN((args->thread_id + 1) * (double)size / args->num_threads, size);

	for (int i = start; i < end; i++) {
		args->target->image[i] = args->source->image[i];
	}

	return NULL;
}

void copy_image(Image *source, Image *target)
{
	pthread_t threads[NUM_THREADS];
    int rc;
    void *status;
    args_copy_image arg[NUM_THREADS];

    for (int i = 0; i < NUM_THREADS; ++i) {
        arg[i].source = source;
        arg[i].target = target;
		arg[i].num_threads = NUM_THREADS;
		arg[i].thread_id = i;
    }

    for (int i = 0; i < NUM_THREADS; ++i) {
        rc = pthread_create(&threads[i], NULL, copy_image_tfunc, &arg[i]);

		if (rc) {
			printf("Couldn't create thread %d\n", i);
			exit(-1);
		}
	}

	for (int i = 0; i < NUM_THREADS; ++i) {
		rc = pthread_join(threads[i], &status);

		if (rc) {
			printf("Couldn't join thread %d\n", i);
			exit(-1);
		}
	}
}

void box_blur(Image *source, Image *target, int w, int h, int radius)
{
	copy_image(source, target);
	box_blur_h(target, source, w, h, radius);
	box_blur_t(source, target, w, h, radius);
}

/**
    Performs a guassian blur operation in O(n) time.

    @param source the source Image
    @param radius the intensity of the blur

*/
Image *gaussian_blur(Image *source, double radius)
{
    int width = source->width;
    int height = source->height;

    Image *tar = new Image(source->width, source->height, source->components);

    int *bxs = boxes_for_gauss(radius, 3);
    box_blur(source, tar, width, height, (bxs[0] - 1) / 2);
    box_blur(tar, source, width, height, (bxs[1] - 1) / 2);
    box_blur(source, tar, width, height, (bxs[2] - 1) / 2);

    delete[] bxs;

    return tar;
}

void mpi_run(int argc, char *argv[])
{
    int num_edge_rows, edge_offset;
    int row_size, image_size;
    int *send_counts, *displs;

    Image *image, *slice, *blurred_slice;

    int blur_intensity;
    int width, height, components;

    // Initialize MPI
    int rank, num_tasks;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);

    blur_intensity = atoi(argv[3]);

    // Calculate the number of edge rows sent
    if (num_tasks > 1) {
        num_edge_rows = blur_intensity * 3;  // Magic number 3 (`n` from `boxes_for_gauss`)
    } else {
        num_edge_rows = 0;
        edge_offset = 0;
    }

    send_counts = new int[num_tasks]();
    displs = new int[num_tasks]();

    if (rank == ROOT) {
        image = new Image(argv[1]);
        width = image->width;
        height = image->height;
        components = image->components;

        if (blur_intensity > height / 2) {
            fprintf(stderr, "Blur intensity must not surpass half of image height\n");
            MPI_Finalize();
            exit(0);
        }

        int num_rows;
        row_size = width * components;
        // Calculate for each processes the image slice size and offset
        for (int i = 0; i < num_tasks; i++) {
            num_rows =
                i < height % num_tasks ? height / num_tasks + 1 : height / num_tasks;
            num_rows += num_edge_rows;

            send_counts[i] = num_rows * row_size;

            // Overlap slices to eliminate edge rows after blurring
            displs[i] =
                i == 0 ? 0 : displs[i - 1] + send_counts[i - 1] - (num_edge_rows * row_size);
            // The first process after ROOT (which has displacement 0) and the last process
            // should have decreased displacements
            if (num_tasks > 1 && (i == ROOT + 1 || i == num_tasks - 1)) {
                if (num_tasks == 2) {
                    displs[i] -= num_edge_rows * row_size;
                } else {
                    displs[i] -= (num_edge_rows / 2) * row_size;
                }
            }
        }
    }

    MPI_Bcast(&width, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&components, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    slice = new Image(width, height, components);

    row_size = width * components;
    image_size = height * row_size;


    MPI_Scatterv(image->image, send_counts, displs, MPI_BYTE,
                 slice->image, image_size, MPI_BYTE,
                 ROOT, MPI_COMM_WORLD);


    blurred_slice = gaussian_blur(slice, blur_intensity);

    // Ignore edges
    if (rank == ROOT) {
        // Recalculate displacement and slice size for each process
        int offset = (num_edge_rows / 2) * row_size;
        for (int i = 1; i < num_tasks; i++) {
            send_counts[i] -= offset;
            displs[i] += offset;
        }
    } else {
        edge_offset = (num_edge_rows / 2) * row_size;
        image_size -= edge_offset;
    }

    MPI_Gatherv(blurred_slice->image + edge_offset, image_size,
            MPI_BYTE, image->image, send_counts, displs,
            MPI_BYTE, ROOT, MPI_COMM_WORLD);


    if (rank == ROOT) {
        image->write_jpeg(argv[2]);

        delete image;
    }

    delete slice;
    delete blurred_slice;
    delete[] displs;
    delete[] send_counts;

    MPI_Finalize();
}

int main(int argc, char *argv[])
{
	if (argc < 5) {
		printf("./gaussian_blur img_in img_out blur_intensity num_threads\n");
		return 0;
	}

	Image *image = new Image(argv[1]);
	int blur_intensity = atoi(argv[3]);

	if (blur_intensity > (image->height / 2)) {
		fprintf(stderr, "Blur intensity must not surpass half of image height\n");
		exit(0);
	}

	NUM_THREADS = strtol(argv[4], NULL, 10);

	mpi_run(argc, argv);
	return 0;
}
