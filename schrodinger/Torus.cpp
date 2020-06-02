#include "Torus.h"

Torus::Torus( int resx, int resy, int resz, float size )
{
	gridx = resx;
	gridy = resy;
	gridz = resz;

	dx = size;
	dy = size;
	dz = size;

	sizex = size * resx;
	sizey = size * resy;
	sizez = size * resz;
}

void Torus::derivativeOfFunction()
{
}

void Torus::derivativeOfOneForm()
{
}

void Torus::derivativeOfTwoForm()
{
}

void Torus::Div()
{
}

void Torus::Sharp()
{
}

void Torus::staggerSharp()
{
}

void Torus::poissonSolve()
{
}
