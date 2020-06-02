#include "renderer.h"

Renderer::Renderer(float* density)
{
	glCullFace(GL_FRONT);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	_density = density;
}

Renderer::~Renderer() {}

void Renderer::Render(void)
{
	GLdouble mvMatrix[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);	
	glBegin(GL_POINTS);
	GRID_LOOP_BIGIN
	{//density > 0 then draw
		if (!FLOAT_EQUAL(_density[_IndexTool(i,j,k)], 0))
		{
			glVertex3f(((float)i - RES_X / 2) / RES_X * 4, ((float)j - RES_Y / 2) / RES_X * 4, ((float)k - RES_Z / 2) / RES_X * 4);
		}
	} GRID_LOOP_END
	glEnd();
}
