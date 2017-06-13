#include "ImageInterpolation.h"
#include "ColorSpaces.h"
#include <math.h>
#include "ImageFilter.h"
#include <iostream>


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
			
			int p, q;

			//ako sam ispao iz slike
			p = (i - 1) / s_x;
			q = (j - 1) / s_y;
			
			//nadji nove vrednosti, po formuli sa pdf-a, ako nisa izasao iz granica
			if (p < xSize - 1) {
				p += 1;
			}
			if (q < ySize - 1) {
				q += 1;
			}
			
			

			y_new[j * newXSize + i] = y_old[q * xSize + p];

		}
	}
	//za u i v komponente
	for (int i = 0; i < newXSize/2; i++) {
		for (int j = 0; j < newYSize/2; j++) {

			int p, q;

			//ako sam ispao iz slike
			p = (i - 1) / s_x;
			q = (j - 1) / s_y;
			
			
			//nadji nove vrednosti, po formuli sa pdf-a, ako nisa izasao iz granica
			if (p < xSize/2 - 1) {
				p += 1;
			}
			if (q < ySize/2 - 1) {
				q += 1;
			}

			u_new[j * newXSize/2 + i] = u_old[q * xSize/2 + p];
			v_new[j * newXSize/2 + i] = v_old[q * xSize/2 + p];


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

			//m je po y n je po x
			int m, n;
			int m_prim, n_prim;

			
			double a = (i / s_x) - floor(i / s_x);
			double b = (j / s_y) - floor(j / s_y);
			
			//ako sam ispao iz slike
			n = i / s_x;
			m = j / s_y;

			n_prim = n;
			m_prim = m;

			//nadji nove vrednosti, po formuli sa pdf-a, ako nisa izasao iz granica
			if (n < xSize - 1) {
				n_prim = n + 1;
			}
			if (m < ySize - 1) {
				m_prim = m + 1;
			}



			y_new[j * newXSize + i] = 
				(1 - a) * (1 - b) * y_old[m * xSize + n] +
				(1 - a) * b * y_old[m_prim * xSize + n] +
				a * (1 - b) * y_old[m * xSize + n_prim] +
				a * b * y_old[m_prim * xSize + n_prim];


		}
	}

	for (int i = 0; i < newXSize/2; i++) {
		for (int j = 0; j < newYSize/2; j++) {

			//m je po y n je po x
			int m, n;
			int m_prim, n_prim;


			double a = (i / s_x) - floor(i / s_x);
			double b = (j / s_y) - floor(j / s_y);

			//ako sam ispao iz slike
			n = i / s_x;
			m = j / s_y;

			n_prim = n;
			m_prim = m;

			//nadji nove vrednosti, po formuli sa pdf-a, ako nisa izasao iz granica
			if (n < xSize/2 - 1) {
				n_prim = n + 1;
			}
			if (m < ySize/2 - 1) {
				m_prim = m + 1;
			}



			u_new[j * newXSize/2 + i] =
				(1 - a) * (1 - b) * u_old[m * xSize/2 + n] +
				(1 - a) * b * u_old[m_prim * xSize/2 + n] +
				a * (1 - b) * u_old[m * xSize/2 + n_prim] +
				a * b * u_old[m_prim * xSize/2 + n_prim];

			v_new[j * newXSize/2 + i] =
				(1 - a) * (1 - b) * v_old[m * xSize/2 + n] +
				(1 - a) * b * v_old[m_prim * xSize/2 + n] +
				a * (1 - b) * v_old[m * xSize/2 + n_prim] +
				a * b * v_old[m_prim * xSize/2 + n_prim];
			

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

void imageRotateBilinear(const uchar input[], int xSize, int ySize, uchar output[], int m_rot, int n_rot, double angle)
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

	double ugao = 3.14159265359 * angle / 180;

	for (int i = 0; i < xSize; i++) {
		for (int j = 0; j < ySize; j++) {

			int n;
			int m;
			int n1;
			int m1;

			n = (i * cos(ugao) - j * sin(ugao) - m_rot * cos(ugao) + n_rot * sin(ugao) + m_rot);
			m = (j * cos(ugao) + i * sin(ugao) - m_rot * sin(ugao) - n_rot * cos(ugao) + n_rot);
	
			double a = n - floor(n);
			double b = m - floor(m);

			n1 = n;
			m1 = m;

			if (m < ySize - 1)
				m1++;

			if (n < xSize - 1)
				n1++;

			
			if (n < 0 || n >= xSize || m < 0 || m >= ySize)
				y_new[j * xSize + i] = 0;
			
			else
			{
				y_new[j * xSize + i] =
					(1 - a) * (1 - b) * y_old[m * xSize + n] +
					(1 - a) * b * y_old[m1 * xSize + n] +
					a * (1 - b) * y_old[m * xSize + n1] +
					a * b * y_old[m1  * xSize + n1];

			}
			
			
		}
	}

	
	
	
	
	for (int i = 0; i < xSize/2; i++) {
		for (int j = 0; j < ySize/2; j++) {
			int n;
			int m;
			int n1;
			int m1;

			n = (i * cos(ugao) - j * sin(ugao) - m_rot/2 * cos(ugao) + n_rot/2 * sin(ugao) + m_rot/2);
			m = (j * cos(ugao) + i * sin(ugao) - m_rot/2 * sin(ugao) - n_rot/2 * cos(ugao) + n_rot/2);

			double a = n - floor(n);
			double b = m - floor(m);

			n1 = n;
			m1 = m;

			if (m < ySize/2 - 1)
				m1++;

			if (n < xSize/2 - 1)
				n1++;


			if (n < 0 || n >= xSize / 2 || m < 0 || m >= ySize / 2)
			{
				u_new[j * xSize / 2 + i] = 0;
				v_new[j * xSize / 2 + i] = 0;
			}
			else
			{
				u_new[j * xSize/2 + i] =
					(1 - a) * (1 - b) * u_old[m * xSize/2 + n] +
					(1 - a) * b * u_old[m1 * xSize/2 + n] +
					a * (1 - b) * u_old[m * xSize/2 + n1] +
					a * b * u_old[m1  * xSize/2 + n1];

				v_new[j * xSize / 2 + i] =
					(1 - a) * (1 - b) * v_old[m * xSize / 2 + n] +
					(1 - a) * b * v_old[m1 * xSize / 2 + n] +
					a * (1 - b) * v_old[m * xSize / 2 + n1] +
					a * b * v_old[m1  * xSize / 2 + n1];

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