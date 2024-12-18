#include "image.hpp"

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
    @param savename the filename to save the result into
    @return void

*/
void gaussian_blur(Image *source, double radius, char *savename)
{
	int width = source->width;
	int height = source->height;

	Image *tar = new Image(source->width, source->height, source->components);

	int *bxs = boxes_for_gauss(radius, 3);
	box_blur(source, tar, width, height, (bxs[0] - 1) / 2);
	box_blur(tar, source, width, height, (bxs[1] - 1) / 2);
	box_blur(source, tar, width, height, (bxs[2] - 1) / 2);

    tar->write_jpeg(savename);

    delete tar;

    delete[] bxs;
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

	gaussian_blur(image, blur_intensity, argv[2]);
	return 0;
}
