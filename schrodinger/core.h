#pragma once

#define M_PI 3.141592653
#define INF 0x3f3f3f3f

#define HBAR  0.1
#define PARTICLE_NUM 2000
#define DT  1.0f/24		
#define eps = 1e-8
#define RES_X 128
#define RES_Y 32
#define RES_Z 32
#define BUFFER_SIZE ((RES_X)*(RES_Y)*(RES_Z))
#define _IndexTool(x,y,z) ( ( ( x ) * ( RES_Y ) * ( RES_Z ) ) + ( ( y ) * ( RES_Z ) ) + ( z ) )

#define GRID_LOOP_BIGIN for (int i=0; i<(RES_X); i++) {\
	                        for (int j=0; j<(RES_Y); j++) {\
		                        for (int k=0; k<(RES_Z); k++) {

#define GRID_LOOP_END }}}

#define FLOAT_EQUAL(a, b) ( ( fabs(a-b) < 0.00001f ) ? true : false )

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctype.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#define GLEW_STATIC
#include <GL/glew.h>
#define GLFW_INCLUDE_GLU
#include <glfw3.h>

#include <GL/gl.h>
#include <GL/glu.h>

#include "event_listener.h"
