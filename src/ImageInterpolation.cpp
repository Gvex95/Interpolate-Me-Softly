#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include <math.h>
#include "ImageFilter.h"


void sampleAndHold(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */
	//prvo mi trebaju neophodne stvari za yuv i rgb konverziju

	uchar* y_old = new uchar[xSize*ySize]();
	char* u_old = new char[xSize * ySize / 4]();
	char* v_old = new char[xSize * ySize / 4]();

	uchar* y_new = new uchar[newXSize*newYSize]();
	char* u_new = new char[newXSize*newYSize / 4]();
	char* v_new = new char[newXSize*newYSize / 4]();

	// konverzija
	RGBtoYUV420(input, xSize, ySize, y_old, u_old, v_old);
	//faktor skaliranja
	const double s_x = (double)newXSize / xSize;
	const double s_y = (double)newYSize / ySize;

	// za y komponentu
	for (int i = 0; i < newXSize; i++) {
		for (int j = 0; j < newYSize; j++) {
			
			int new_i, new_j;

			//ako sam ispao iz slike
			new_i = (i - 1) / s_x;
			new_j = (j - 1) / s_y;
			
			//nadji nove vrednosti, po formuli sa pdf-a, ako nisa izasao iz granica
			if (new_i < xSize - 1) {
				new_i += 1;
			}
			if (new_j < ySize - 1) {
				new_j += 1;
			}
			
			

			y_new[j * newXSize + i] = y_old[new_j * xSize + new_i];

		}
	}
	//za u i v komponente
	for (int i = 0; i < newXSize/2; i++) {
		for (int j = 0; j < newYSize/2; j++) {

			int new_i, new_j;

			//ako sam ispao iz slike
			new_i = (i - 1) / s_x;
			new_j = (j - 1) / s_y;
			
			
			//nadji nove vrednosti, po formuli sa pdf-a, ako nisa izasao iz granica
			if (new_i < xSize/2 - 1) {
				new_i += 1;
			}
			if (new_j < ySize/2 - 1) {
				new_j += 1;
			}

			u_new[j * newXSize/2 + i] = u_old[new_j * xSize/2 + new_i];
			v_new[j * newXSize/2 + i] = v_old[new_j * xSize/2 + new_i];


		}
	}

	YUV420toRGB(y_new, u_new, v_new, newXSize, newYSize, output);

	delete[] y_old;
	delete[] u_old;
	delete[] v_old;

	delete[] y_new;
	delete[] u_new;
	delete[] v_new;



}

void bilinearInterpolate(const uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */
	/* TO DO */
	//prvo mi trebaju neophodne stvari za yuv i rgb konverziju

	uchar* y_old = new uchar[xSize*ySize]();
	char* u_old = new char[xSize * ySize / 4]();
	char* v_old = new char[xSize * ySize / 4]();

	uchar* y_new = new uchar[newXSize*newYSize]();
	char* u_new = new char[newXSize*newYSize / 4]();
	char* v_new = new char[newXSize*newYSize / 4]();

	// konverzija
	RGBtoYUV420(input, xSize, ySize, y_old, u_old, v_old);
	//faktor skaliranja
	const double s_x = (double)newXSize / xSize;
	const double s_y = (double)newYSize / ySize;
	//a i b za formulu
	
	for (int i = 0; i < newXSize; i++) {
		for (int j = 0; j < newYSize; j++) {

			int new_i, new_j;
			int new_new_i, new_new_j;

			
			double a = (i / s_x) - floor(i / s_x);
			double b = (j / s_y) - floor(j / s_y);
			
			//ako sam ispao iz slike
			new_i = i / s_x;
			new_j = j / s_y;

			new_new_i = new_i;
			new_new_j = new_j;

			//nadji nove vrednosti, po formuli sa pdf-a, ako nisa izasao iz granica
			if (new_i < xSize - 1) {
				new_new_i = new_i + 1;
			}
			if (new_j < ySize - 1) {
				new_new_j = new_j + 1;
			}



			y_new[j * newXSize + i] = 
				(1 - a) * (1 - b) * y_old[new_j * xSize + new_i] +
				(1 - a) * b * y_old[new_new_j * xSize + new_i] +
				a * (1 - b) * y_old[new_j * xSize + new_new_i] +
				a * b * y_old[new_new_j * xSize + new_new_i];


		}
	}

	for (int i = 0; i < newXSize/2; i++) {
		for (int j = 0; j < newYSize/2; j++) {

			int new_i, new_j;
			int new_new_i, new_new_j;

			double a = (i / s_x) - floor(i / s_x);
			double b = (j / s_y) - floor(j / s_y);
			
			//ako sam ispao iz slike
			new_i = i / s_x;
			new_j = j / s_y;

			new_new_i = new_i;
			new_new_j = new_j;

			//nadji nove vrednosti, po formuli sa pdf-a, ako nisa izasao iz granica
			if (new_i < xSize/2 - 1) {
				new_new_i = new_i + 1;
			}
			if (new_j < ySize/2 - 1) {
				new_new_j = new_j + 1;
			}


			u_new[j * newXSize/2 + i] = 
				(1 - a) * (1 - b) * u_old[new_j * xSize/2 + new_i] +
				(1 - a) * b * u_old[new_new_j * xSize/2 + new_i] +
				a * (1 - b) * u_old[new_j * xSize/2 + new_new_i] +
				a * b * u_old[new_new_j * xSize/2 + new_new_i];

			v_new[j * newXSize/2 + i] = 
				(1 - a) * (1 - b) * v_old[new_j * xSize/2 + new_i] +
				(1 - a) * b * v_old[new_new_j * xSize/2 + new_i] +
				a * (1 - b) * v_old[new_j * xSize/2 + new_new_i] +
				a * b * v_old[new_new_j * xSize/2 + new_new_i];
			

		}
	}

	YUV420toRGB(y_new, u_new, v_new, newXSize, newYSize, output);

	delete[] y_old;
	delete[] u_old;
	delete[] v_old;

	delete[] y_new;
	delete[] u_new;
	delete[] v_new;




}

void bicubicInterpolate(uchar input[], int xSize, int ySize, uchar output[], int newXSize, int newYSize)
{
	/* TO DO */
}

void imageRotate(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	uchar* y_old = new uchar[xSize*ySize]();
	char* u_old = new char[xSize * ySize / 4]();
	char* v_old = new char[xSize * ySize / 4]();

	uchar* y_new = new uchar[xSize*ySize]();
	char* u_new = new char[xSize*ySize / 4]();
	char* v_new = new char[xSize*ySize / 4]();

	// konverzija
	RGBtoYUV420(input, xSize, ySize, y_old, u_old, v_old);

	double ugao = 3.14 * angle / 180;

	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {

			int new_i, new_j;
		
			new_i = (int) (i * cos(ugao) - j * sin(ugao) - m * cos(ugao) + n * sin(ugao) + m);
			new_j = (int) (j * cos(ugao) + i * sin(ugao) - m * sin(ugao) - n * cos(ugao) + n);

			
			if (new_i < 0 || new_i >= xSize || new_j < 0 || new_j >= ySize)
				y_new[j * xSize + i] = 0;
			else
				y_new[j * xSize + i] = y_old[new_j * xSize + new_i];
		}
	}

	for (int i = 0; i < xSize/2; i++) {
		for (int j = 0; j < ySize/2; j++) {

			int new_i, new_j;

			new_i = (int) (i * cos(ugao) - j * sin(ugao) - m/2 * cos(ugao) + n/2 * sin(ugao) + m/2);
			new_j = (int) (j * cos(ugao) + i * sin(ugao) - m/2 * sin(ugao) - n/2 * cos(ugao) + n/2);


			if (new_i < 0 || new_i >= xSize/2 || new_j < 0 || new_j >= ySize/2)
			{ 
				u_new[j * xSize/2 + i] = 0;
				v_new[j * xSize/2 + i] = 0;
			}
			else
			{
				u_new[j * xSize/2 + i] = u_old[new_j * xSize/2 + new_i];
				v_new[j * xSize/2 + i] = v_old[new_j * xSize/2 + new_i];
			}
		}
	}

	YUV420toRGB(y_new, u_new, v_new, xSize, ySize, output);

	delete[] y_old;
	delete[] u_old;
	delete[] v_old;

	delete[] y_new;
	delete[] u_new;
	delete[] v_new;


}

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m, int n, double angle)
{
	/* TO DO */
	uchar* y_old = new uchar[xSize*ySize]();
	char* u_old = new char[xSize * ySize / 4]();
	char* v_old = new char[xSize * ySize / 4]();

	uchar* y_new = new uchar[xSize*ySize]();
	char* u_new = new char[xSize*ySize / 4]();
	char* v_new = new char[xSize*ySize / 4]();

	// konverzija
	RGBtoYUV420(input, xSize, ySize, y_old, u_old, v_old);

	uchar *y_extended = new uchar[(xSize + 2) * (ySize + 2)];

	extendBorders(y_old, xSize, ySize, y_extended, 1);

	double ugao = 3.14 * angle / 180;

	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {

			int new_i, new_j;
			int new_new_i, new_new_j;

			new_i = (int)(i * cos(ugao) - j * sin(ugao) - m * cos(ugao) + n * sin(ugao) + m);
			new_j = (int)(j * cos(ugao) + i * sin(ugao) - m * sin(ugao) - n * cos(ugao) + n);

			double a = i * cos(ugao) - j * sin(ugao) - m * cos(ugao) + n * sin(ugao) + m - new_i;
			double b = j * cos(ugao) + i * sin(ugao) - m * sin(ugao) - n * cos(ugao) + n - new_j;

			new_new_i = new_i + 1;
			new_new_j = new_j + 1;
			
			
			if (new_i < 0 || new_i >= xSize || new_j < 0 || new_j >= ySize)
				y_new[j * xSize + i] = 0;
			else
			{
				y_new[j * xSize + i] =
					(1 - a) * (1 - b) * y_extended[new_j * (xSize + 2) + new_i] +
					(1 - a) * b * y_extended[new_new_j * (xSize + 2) + new_i] +
					a * (1 - b) * y_extended[new_j * (xSize + 2) + new_new_i] +
					a * b * y_extended[new_new_j * (xSize + 2) + new_new_i];
			}
		}
	}

	
	
	char *u_extended = new char[(xSize + 8) * (ySize + 8) / 4]();
	char *v_extended = new char[(xSize + 8) * (ySize + 8) / 4]();

	extendBorders_char(u_old, xSize / 2, ySize / 2, u_extended, 1);
	extendBorders_char(v_old, xSize / 2, ySize / 2, v_extended, 1);
	
	for (int i = 0; i < xSize / 2; i++) {
		for (int j = 0; j < ySize / 2; j++) {

			int new_i, new_j;
			int new_new_i, new_new_j;

			new_i = (int)(i * cos(ugao) - j * sin(ugao) - m / 2 * cos(ugao) + n / 2 * sin(ugao) + m / 2);
			new_j = (int)(j * cos(ugao) + i * sin(ugao) - m / 2 * sin(ugao) - n / 2 * cos(ugao) + n / 2);

			double a = i * cos(ugao) - j * sin(ugao) - m/2 * cos(ugao) + n/2 * sin(ugao) + m/2 - new_i;
			double b = j * cos(ugao) + i * sin(ugao) - m/2 * sin(ugao) - n/2 * cos(ugao) + n/2 - new_j;

			new_new_i = new_i + 1;
			new_new_j = new_j + 1;


			if (new_i < 0 || new_i >= xSize / 2 || new_j < 0 || new_j >= ySize / 2)
			{
				u_new[j * xSize / 2 + i] = 0;
				v_new[j * xSize / 2 + i] = 0;
			}
			else
			{
				u_new[j * xSize/2 + i] =
					(1 - a) * (1 - b) * u_extended[new_j * (xSize/2 + 2) + new_i] +
					(1 - a) * b * u_extended[new_new_j * (xSize/2 + 2) + new_i] +
					a * (1 - b) * u_extended[new_j * (xSize/2 + 2) + new_new_i] +
					a * b * u_extended[new_new_j * (xSize/2 + 2) + new_new_i];


				v_new[j * xSize/2 + i] =
					(1 - a) * (1 - b) * v_extended[new_j * (xSize/2 + 2) + new_i] +
					(1 - a) * b * v_extended[new_new_j * (xSize/2 + 2) + new_i] +
					a * (1 - b) * v_extended[new_j * (xSize/2 + 2) + new_new_i] +
					a * b * v_extended[new_new_j * (xSize/2 + 2) + new_new_i];
			}
		}
	}

	YUV420toRGB(y_new, u_new, v_new, xSize, ySize, output);

	delete[] y_old;
	delete[] u_old;
	delete[] v_old;

	delete[] y_new;
	delete[] u_new;
	delete[] v_new;

}