#include "image.hpp"
#define C char
#define UC unsigned C

using namespace std;
class BGR : public Image
{
      public:
	UC **B;
	UC **G;
	UC **R;

	BGR(Image &source) : Image(source)
	{
		B = (UC **)malloc(sizeof(UC *) * height);
		G = (UC **)malloc(sizeof(UC *) * height);
		R = (UC **)malloc(sizeof(UC *) * height);
		int i, j, count = 0;
		for (i = 0; i < height; i++) {
			R[i] = (UC *)malloc(sizeof(UC) * width);
			G[i] = (UC *)malloc(sizeof(UC) * width);
			B[i] = (UC *)malloc(sizeof(UC) * width);
			for (j = 0; j < width; j++) {
				R[i][j] = image[count];
				G[i][j] = image[count + 1];
				B[i][j] = image[count + 2];
				count += components;
			}
		}
	}

	BGR(int width, int height, int components)
	    : Image(width, height, components)
	{
		B = (UC **)malloc(sizeof(UC *) * height);
		G = (UC **)malloc(sizeof(UC *) * height);
		R = (UC **)malloc(sizeof(UC *) * height);
		int i, j, count = 0;
		for (i = 0; i < height; i++) {
			R[i] = (UC *)malloc(sizeof(UC) * width);
			G[i] = (UC *)malloc(sizeof(UC) * width);
			B[i] = (UC *)malloc(sizeof(UC) * width);
			for (j = 0; j < width; j++) {
				R[i][j] = image[count];
				G[i][j] = image[count + 1];
				B[i][j] = image[count + 2];
				count += components;
			}
		}
	}

	BGR(int width, int height, int components, unsigned char *image)
	    : Image(width, height, components, image)
	{
		B = (UC **)malloc(sizeof(UC *) * height);
		G = (UC **)malloc(sizeof(UC *) * height);
		R = (UC **)malloc(sizeof(UC *) * height);
		int i, j, count = 0;
		for (i = 0; i < height; i++) {
			R[i] = (UC *)malloc(sizeof(UC) * width);
			G[i] = (UC *)malloc(sizeof(UC) * width);
			B[i] = (UC *)malloc(sizeof(UC) * width);
			for (j = 0; j < width; j++) {
				R[i][j] = image[count];
				G[i][j] = image[count + 1];
				B[i][j] = image[count + 2];
				count += components;
			}
		}
	}

	BGR(const BGR &source)
	    : Image(source.width, source.height, source.components,
		    source.image)
	{
		B = (UC **)malloc(sizeof(UC *) * height);
		G = (UC **)malloc(sizeof(UC *) * height);
		R = (UC **)malloc(sizeof(UC *) * height);
		int i, j, count = 0;
		for (i = 0; i < height; i++) {
			R[i] = (UC *)malloc(sizeof(UC) * width);
			G[i] = (UC *)malloc(sizeof(UC) * width);
			B[i] = (UC *)malloc(sizeof(UC) * width);
			for (j = 0; j < width; j++) {
				R[i][j] = image[count];
				G[i][j] = image[count + 1];
				B[i][j] = image[count + 2];
				count += components;
			}
		}
	}

	void save(char *filename)
	{

		int count = 0;
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				image[count++] = R[i][j];
				image[count++] = G[i][j];
				image[count++] = B[i][j];
			}
		}

		write_jpeg(filename);
	}

	Image *save_to_img()
	{
		int count = 0;
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				image[count++] = R[i][j];
				image[count++] = G[i][j];
				image[count++] = B[i][j];
			}
		}

		return new Image(this->width, this->height, this->components, this->image);
	}

	~BGR()
	{
		for (int i = 0; i < height; i++) {
			free(R[i]);
			free(G[i]);
			free(B[i]);
		}

		free(R);
		free(G);
		free(B);
	}
};
