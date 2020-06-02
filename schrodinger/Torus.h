#pragma once

#include "ComplexWrapper.h"
#include "core.h"

class Torus
{
public:
	//size [ -1.0, 1.0 ]
	int gridx;
	int gridy;
	int gridz;

	float sizex;
	float sizey;
	float sizez;
	 
	float dx;
	float dy;
	float dz;
	Torus() {}
	Torus( int resx, int resy, int resz, float size = 0.04 );
	void derivativeOfFunction();
	void derivativeOfOneForm();
	void derivativeOfTwoForm();
	void Div();
	void Sharp();
	void staggerSharp();
	void poissonSolve();
};

