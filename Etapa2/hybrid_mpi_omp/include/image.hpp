#include <bits/stdc++.h>
#include <jpeglib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <typeinfo>
#include <unistd.h>
#define GRAYSCALE JCS_GRAYSCALE
#define RGB JCS_RGB

using namespace std;
/* High level class Image */

enum save_flags { SAVE_FIRST, DO_NOT_SAVE };

class Image
{

      public:
	int width;
	int height;
	int components;
	unsigned char *image; // contains image pixels in form R,G,B,R,G,B,...
	char *filename;

	Image()
	{
		filename = NULL;
		width = 0;
		height = 0;
		components = 0;
		image = NULL;
	}

	Image(const char *name) // constructor takes filepath as parameter
	{
		int size = strlen(name) + 1;
		filename = new char[size];
		strcpy(filename, name);
		read_jpeg();
	}

	Image(int wid, int hei, int component, unsigned char *img)
	{
		width = wid;
		height = hei;
		components = component;
		filename = NULL;

		int size = width * height * components;
		image = new unsigned char[size];

		int i;
		for (i = 0; i < size; ++i)
			image[i] = img[i];
	}

	Image(int wid, int hei, int component)
	{

		width = wid;
		height = hei;
		components = component;
		filename = NULL;
		image = new unsigned char[width * height * components]();
	}

	Image(const Image &src)
	{

		width = src.width;
		height = src.height;
		components = src.components;
		int len = strlen(src.filename) + 1;
		if (len > 0) {
			filename = new char[len];
			strcpy(filename, src.filename);
		}

		int size = width * height * components;
		image = new unsigned char[size];

		int i;
		for (i = 0; i < size; ++i)
			image[i] = src.image[i];
	}

	unsigned char **read_jpeg() // puts jpeg filename in image
	{
		struct jpeg_decompress_struct cinfo;
		FILE *infile;
		struct jpeg_error_mgr jerr;

		infile = fopen(filename, "rb");
		cinfo.err = jpeg_std_error(&jerr);
		jpeg_create_decompress(&cinfo);

		jpeg_stdio_src(&cinfo, infile);
		jpeg_read_header(&cinfo, TRUE);
		jpeg_start_decompress(&cinfo);

		width = static_cast<int>(cinfo.output_width);
		height = static_cast<int>(cinfo.output_height);
		components = cinfo.output_components;

		image = new unsigned char[height * width * components];

		int row_stride = width * components;
		int iter = 0;
		while (cinfo.output_scanline < cinfo.output_height) {
			(void)jpeg_read_scanlines(&cinfo, &image, 1);
			image += row_stride;
			iter += 1;
		}

		image -= iter * row_stride;

		jpeg_destroy_decompress(&cinfo);
		fclose(infile);
		return 0;
	}

	void write_jpeg(const char *savename, J_COLOR_SPACE colorspace = RGB,
			int quality = 98) // writes image with name savename
					  // with passed quality (1-100)
	{

		struct jpeg_compress_struct cinfo;
		struct jpeg_error_mgr jerr;

		FILE *outfile;		 /* target file */
		JSAMPROW row_pointer[1]; /* pointer to JSAMPLE row[s] */
		int row_stride; /* physical row width in image buffer */

		cinfo.err = jpeg_std_error(&jerr);

		jpeg_create_compress(&cinfo);

		if ((outfile = fopen(savename, "wb")) == NULL) {
			fprintf(stderr, "can't open %s\n", savename);
			exit(1);
		}
		jpeg_stdio_dest(&cinfo, outfile);

		cinfo.image_width =
		    width; /* image width and height, in pixels */
		cinfo.image_height = height;
		cinfo.input_components =
		    components; /* # of color components per pixel */
		cinfo.in_color_space =
		    colorspace; /* colorspace of input image */

		jpeg_set_defaults(&cinfo);

		jpeg_set_quality(&cinfo, quality,
				 TRUE /* limit to baseline-JPEG values */);

		jpeg_start_compress(&cinfo, TRUE);

		row_stride = width * components; /* JSAMPLEs per row in image */

		while (cinfo.next_scanline < cinfo.image_height) {
			row_pointer[0] =
			    &image[cinfo.next_scanline * row_stride];
			(void)jpeg_write_scanlines(&cinfo, row_pointer, 1);
		}

		jpeg_finish_compress(&cinfo);
		fclose(outfile);
		jpeg_destroy_compress(&cinfo);
	}

	void display(
	    const char *savename, int save_flag = DO_NOT_SAVE,
	    J_COLOR_SPACE colorspace =
		RGB) // function to display image after saving at path savename
	{
		if (save_flag == SAVE_FIRST)
			write_jpeg(savename, colorspace); // saves jpeg first
		/* eog <savename> & - displays image on linux machines in
		background so that execution can continue */
		char eog[] = "eog ";
		char *command = new char[10 + strlen(savename) + 1];
		strcpy(command, eog);
		strcat(command, savename);
		char temp[] = " &";
		strcat(command, temp);
		system(command);
	}

	~Image()
	{
		if (image)
			delete[] image;

		delete[] filename;
	}
};
