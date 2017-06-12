
#include "ImageProcessing.h"
#include "ImageInterpolation.h"

#include <QDebug>

int newSize(int x, int y)
{
	return (x + y - 1) & ~(y - 1);
}

void imageProcessingFun(const QString& progName, QImage* const outImgs, const QImage* const inImgs, const QVector<double>& params) 
{
	int X_SIZE = inImgs->width();
	int Y_SIZE = inImgs->height();
	int NEW_X_SIZE, NEW_Y_SIZE;

	/* NOTE: Calculate output image resolution and construct output image object */

	if(progName == "Sample and hold") 
	{	
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */
		NEW_X_SIZE = newSize(X_SIZE * params[1], 4);
		NEW_Y_SIZE = newSize(Y_SIZE * params[0], 4);

		/* TO DO: Calculate output image resolution and construct output image object */
		new (outImgs) QImage(NEW_X_SIZE, NEW_Y_SIZE, inImgs->format());

		/* TO DO: Perform Sample and hold interpolation  */
		sampleAndHold(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), NEW_X_SIZE, NEW_Y_SIZE);


	}
	else if (progName == "Bilinear") 
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		NEW_X_SIZE = newSize(X_SIZE * params[1], 4);
		NEW_Y_SIZE = newSize(Y_SIZE * params[0], 4);
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */

		/* TO DO: Calculate output image resolution and construct output image object */
		new (outImgs) QImage(NEW_X_SIZE, NEW_Y_SIZE, inImgs->format());
		/* TO DO: Perform Bilinear interpolation  */
		bilinearInterpolate(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), NEW_X_SIZE, NEW_Y_SIZE);
	}
	else if (progName == "Bicubic")
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Vertical scale factor is params[0] */
		/* Horizontal scale factor is params[1] */
		/* TO DO: Calculate output image resolution and construct output image object */
		/* TO DO: Perform Bicubic interpolation  */
	}
	else if(progName == "Rotation") 
	{	
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Rotation angle in degrees is params[0]*/
		/* Center of rotation coordinates are (XSIZE/2, YSIZE/2) */
		int m = X_SIZE / 2;
		int n = Y_SIZE / 2;
		int ugao = params[0];

		/* TO DO: Construct output image object */
		new (outImgs) QImage(X_SIZE, Y_SIZE, inImgs->format());

		/* TO DO: Perform image rotation */
		imageRotate(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), m, n, ugao);

	}
	else if (progName == "Rotation Bilinear") 
	{
		/* Input image data in RGB format can be obtained with inImgs->bits() */
		/* Rotation angle in degrees is params[0]*/
		/* Center of rotation coordinates are (XSIZE/2, YSIZE/2) */
		int m = X_SIZE / 2;
		int n = Y_SIZE / 2;
		int ugao = params[0];
		/* TO DO: Construct output image object */
		new (outImgs) QImage(X_SIZE, Y_SIZE, inImgs->format());
		/* TO DO: Perform image rotation with bilinear interpolation */
		imageRotateBilinear(inImgs->bits(), X_SIZE, Y_SIZE, outImgs->bits(), m, n, ugao);
	}

}

