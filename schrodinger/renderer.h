#pragma once
#include "core.h"

class Renderer
{
public:
	float *_density;
	Renderer(float* volumeData);
	~Renderer();
	void Render();
};
