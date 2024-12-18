#include "image.hpp"
#include <omp.h>

int NUM_THREADS;

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

void box_blur_h(Image *source, Image *target, int w, int h, int radius)
{
	double iarr = (double)1 / (radius + radius + 1);

    int chunk_size = h / NUM_THREADS;
    #pragma omp parallel for schedule(static, chunk_size)
	for (int i = 0; i < h; i++) {
		int ti = i * w;
		int li = ti;
		int ri = ti + radius;
		int fv = ti * source->components;
		int lv = (ti + w - 1) * source->components;

		unsigned currennt_r = source->image[fv] * (radius + 1);
		unsigned currennt_g = source->image[fv + 1] * (radius + 1);
		unsigned currennt_b = source->image[fv + 2] * (radius + 1);

		int pixel;
		for (int j = 0; j < radius; j++) {
			pixel = (ti + j) * source->components;
			currennt_r += source->image[pixel];
			currennt_g += source->image[pixel + 1];
			currennt_b += source->image[pixel + 2];
		}

		for (int j = 0; j <= radius; j++) {
			pixel = ri * source->components;
			ri++;
			currennt_r += (source->image[pixel] - source->image[fv]);
			currennt_g += (source->image[pixel + 1] - source->image[fv + 1]);
			currennt_b += (source->image[pixel + 2] - source->image[fv + 2]);

			set_pixel(target, ti++, currennt_r * iarr,
				  currennt_g * iarr, currennt_b * iarr);
		}

		for (int j = radius + 1; j < w - radius; j++) {
			int first_pixel = ri * source->components;
			int second_pixel = li * source->components;
			ri++;
			li++;

			currennt_r += (source->image[first_pixel] -
				       source->image[second_pixel]);
			currennt_g += (source->image[first_pixel + 1] -
				       source->image[second_pixel + 1]);
			currennt_b += (source->image[first_pixel + 2] -
				       source->image[second_pixel + 2]);

			set_pixel(target, ti++, currennt_r * iarr,
				  currennt_g * iarr, currennt_b * iarr);
		}

		for (int j = w - radius; j < w; j++) {
			pixel = li * source->components;
			li++;

			currennt_r += (source->image[lv] - source->image[pixel]);
			currennt_g += (source->image[lv + 1] - source->image[pixel + 1]);
			currennt_b += (source->image[lv + 2] - source->image[pixel + 2]);

			set_pixel(target, ti++, currennt_r * iarr,
				  currennt_g * iarr, currennt_b * iarr);
		}
	}
}


void box_blur_t(Image *source, Image *target, int w, int h, int radius)
{

	double iarr = (double)1 / (radius + radius + 1);

    int chunk_size = w / NUM_THREADS;
    #pragma omp parallel for schedule(static, chunk_size)
	for (int i = 0; i < w; i++) {
		int ti = i;
		int li = ti;
		int ri = ti + radius * w;

		int fv = ti * source->components;
		int lv = (ti + w * (h - 1)) * source->components;

		unsigned currennt_r = source->image[fv] * (radius + 1);
		unsigned currennt_g = source->image[fv + 1] * (radius + 1);
		unsigned currennt_b = source->image[fv + 2] * (radius + 1);

		int pixel;
		
		for (int j = 0; j < radius; j++) {
			pixel = (ti + j * w) * source->components;
			currennt_r += source->image[pixel];
			currennt_g += source->image[pixel + 1];
			currennt_b += source->image[pixel + 2];
		}

		for (int j = 0; j <= radius; j++) {
			pixel = ri * source->components;
			currennt_r += (source->image[pixel] - source->image[fv]);
			currennt_g += (source->image[pixel + 1] - source->image[fv + 1]);
			currennt_b += (source->image[pixel + 2] - source->image[fv + 2]);

			set_pixel(target, ti, currennt_r * iarr,
				  currennt_g * iarr, currennt_b * iarr);

			ri += w;
			ti += w;
		}

		for (int j = radius + 1; j < h - radius; j++) {
			int first_pixel = ri * source->components;
			int second_pixel = li * source->components;

			currennt_r += (source->image[first_pixel] -
				       source->image[second_pixel]);
			currennt_g += (source->image[first_pixel + 1] -
				       source->image[second_pixel + 1]);
			currennt_b += (source->image[first_pixel + 2] -
				       source->image[second_pixel + 2]);

			set_pixel(target, ti, currennt_r * iarr,
				  currennt_g * iarr, currennt_b * iarr);

			li += w;
			ri += w;
			ti += w;
		}

		for (int j = h - radius; j < h; j++) {
			int pixel = li * source->components;

			currennt_r += (source->image[lv] - source->image[pixel]);
			currennt_g += (source->image[lv + 1] - source->image[pixel + 1]);
			currennt_b += (source->image[lv + 2] - source->image[pixel + 2]);

			set_pixel(target, ti, currennt_r * iarr,
				  currennt_g * iarr, currennt_b * iarr);

			li += w;
			ti += w;
		}
	}
}

void copy_image(Image *source, Image *target)
{
	int size = source->height * source->width * source->components;

    int chunk_size = size / NUM_THREADS;
    #pragma omp parallel for schedule(static, chunk_size)
	for (int i = 0; i < size; i++) {
		target->image[i] = source->image[i];
		
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
	omp_set_num_threads(NUM_THREADS);

	gaussian_blur(image, blur_intensity, argv[2]);
	return 0;
}
