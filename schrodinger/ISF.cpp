#include "core.h"
#include "ISF.h"

ISF::ISF()
{
	_isLeftKeyPressed = false;
	Reset();
}

ISF::~ISF()
{
	delete[] _density;
	delete[] _velX;
	delete[] _velY;
	delete[] _velZ;
	delete[] Phase;
	delete[] PoissonParameters;
}



void ISF::InitSmoke()
{
	memset(_density, 0, BUFFER_SIZE * sizeof(float));
	const int centerY = RES_Y / 2;
	const int centerZ = RES_Z / 2;
	float dens = 1;

	for (int i = 0; i < PARTICLE_NUM; i++)
	{
		double rt = (double)(i) / PARTICLE_NUM;
		rt *= 2 * M_PI;

		double ix = 1;
		double jy = centerY + 0.7 * RES_Y / 2 * cos(rt);
		double kz = centerZ + 0.7 * RES_Z / 2 * sin(rt);

		particleX[i] = ix * _Torus.dx;
		particleY[i] = jy * _Torus.dy;
		particleZ[i] = kz * _Torus.dz;

		particleGridX[i] = (int)ix;
		particleGridY[i] = (int)jy;
		particleGridZ[i] = (int)kz;

		this->_density[_IndexTool((int)ix, (int)jy, (int)kz)] += dens;
		this->_velX[_IndexTool((int)ix, (int)jy, (int)kz)] = 0.6;
		particleDie[i] = 0;
	}
	for (int i = 0; i < BUFFER_SIZE; i++)
		_velX[i] = 0.6;
}

void ISF::BuildISF()
{
	_density = new float[BUFFER_SIZE];
	for (int i = 0; i < BUFFER_SIZE; i++)
		_density[i] = .0;
	_velX = new float[BUFFER_SIZE];
	for (int i = 0; i < BUFFER_SIZE; i++)
		_velX[i] = .0;
	_velY = new float[BUFFER_SIZE];
	for (int i = 0; i < BUFFER_SIZE; i++)
		_velY[i] = .0;
	_velZ = new float[BUFFER_SIZE];
	for (int i = 0; i < BUFFER_SIZE; i++)
		_velZ[i] = .0;
	Phase = new float[BUFFER_SIZE];
	for (int i = 0; i < BUFFER_SIZE; i++)
		Phase[i] = .0;
	PoissonParameters = new float[BUFFER_SIZE];
	for (int i = 0; i < BUFFER_SIZE; i++)
		PoissonParameters[i] = .0;
	_renderer = new Renderer(_density);
	Omega = 0;


	_Torus = Torus(RES_X, RES_Y, RES_Z);
	SchrodingerMask = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * BUFFER_SIZE);
	psi1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * BUFFER_SIZE);
	psi2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * BUFFER_SIZE);
	_DIV = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * BUFFER_SIZE);
	PoissonSolutions = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * BUFFER_SIZE);

	jetVel[0] = 1.0;
	jetVel[1] = 0.0;
	jetVel[2] = 0.0;

	for (int i = 0; i < 3; i++)
	{
		kvec[i] = jetVel[i] / HBAR;
		Omega += jetVel[i] * jetVel[i];
	}
	Omega /= 2.0 * HBAR;

	double nx = _Torus.gridx, ny = _Torus.gridy, nz = _Torus.gridz;
	double fac = -4.0 * M_PI * M_PI * HBAR;

	GRID_LOOP_BIGIN{
		double kx = (i - nx / 2) / _Torus.sizex;
		double ky = (j - ny / 2) / _Torus.sizey;
		double kz = (k - nz / 2) / _Torus.sizez;
		double lambda = fac * (kx * kx + ky * ky + kz * kz);
		fftw_complex tmp;
		tmp[0] = 0;
		tmp[1] = lambda * DT / 2;
		iexp(SchrodingerMask[_IndexTool(i, j, k)], tmp);


		psi1[_IndexTool(i, j, k)][0] = 1.0;
		psi1[_IndexTool(i, j, k)][1] = 0.0;
		psi2[_IndexTool(i, j, k)][0] = 0.01;
		psi2[_IndexTool(i, j, k)][1] = 0.0;

		double psi_norm = sqrt(
			  psi1[_IndexTool(i, j, k)][0] * psi1[_IndexTool(i, j, k)][0]
			+ psi1[_IndexTool(i, j, k)][1] * psi1[_IndexTool(i, j, k)][1]
			+ psi2[_IndexTool(i, j, k)][0] * psi2[_IndexTool(i, j, k)][0]
			+ psi2[_IndexTool(i, j, k)][1] * psi2[_IndexTool(i, j, k)][1]);

		idiv(psi1[_IndexTool(i, j, k)], psi_norm);
		idiv(psi2[_IndexTool(i, j, k)], psi_norm);

		Phase[_IndexTool(i, j, k)] = kvec[0] * i * _Torus.dx + kvec[1] * j * _Torus.dy + kvec[2] * k * _Torus.dz;
		//std::cout << Phase[_IndexTool(i, j, k)] << std::endl;
	}GRID_LOOP_END
}


void ISF::constrainVel(float totaltime)
{
	GRID_LOOP_BIGIN{
		float radius = (j - (RES_Y) / 2.0) * (j - (RES_Y) / 2.0) + (k - (RES_Z) / 2.0) * (k - (RES_Z) / 2.0);
		if (i <= RES_X/8.0 && radius <= 0.49 * RES_Y / 2.0 * RES_Y / 2.0)
		{
			double amp1 = sqrt(psi1[_IndexTool(i, j, k)][0] * psi1[_IndexTool(i, j, k)][0] + psi1[_IndexTool(i, j, k)][1] * psi1[_IndexTool(i, j, k)][1]);
			double amp2 = sqrt(psi2[_IndexTool(i, j, k)][0] * psi2[_IndexTool(i, j, k)][0] + psi2[_IndexTool(i, j, k)][1] * psi2[_IndexTool(i, j, k)][1]);
			fftw_complex tmp;
			tmp[0] = 0.0;
			tmp[1] = Phase[_IndexTool(i, j, k)] - Omega * totaltime;

			iexp(psi1[_IndexTool(i, j, k)], tmp);
			imul(psi1[_IndexTool(i, j, k)], amp1);

			iexp(psi2[_IndexTool(i, j, k)], tmp);
			imul(psi2[_IndexTool(i, j, k)], amp2);
		}
	}GRID_LOOP_END

	pressureProject();
}

void ISF::pressureProject()
{
	velocityOneForm(1.0);
	divergence();
	poissonSolve();
	GaugeTransform();
	//psiNormalize();
}

void ISF::psiNormalize()
{
	GRID_LOOP_BIGIN{
	double psi_norm = sqrt(
		  psi1[_IndexTool(i, j, k)][0] * psi1[_IndexTool(i, j, k)][0]
		+ psi1[_IndexTool(i, j, k)][1] * psi1[_IndexTool(i, j, k)][1]
		+ psi2[_IndexTool(i, j, k)][0] * psi2[_IndexTool(i, j, k)][0]
		+ psi2[_IndexTool(i, j, k)][1] * psi2[_IndexTool(i, j, k)][1]);

	idiv(psi1[_IndexTool(i, j, k)], psi_norm);
	idiv(psi2[_IndexTool(i, j, k)], psi_norm);

	}GRID_LOOP_END
}

void ISF::velocityOneForm(float localhbar)
{
	GRID_LOOP_BIGIN{
		int ixp = (i + 1) % _Torus.gridx;
		int iyp = (j + 1) % _Torus.gridy;
		int izp = (k + 1) % _Torus.gridz;

		fftw_complex vx,vy,vz;
		fftw_complex tmpres1, tmpres2;
		fftw_complex tmpl1, tmpl2;
		tmpl1[0] = psi1[_IndexTool(i, j, k)][0];
		tmpl1[1] = -psi1[_IndexTool(i, j, k)][1];

		tmpl2[0] = psi2[_IndexTool(i, j, k)][0];
		tmpl2[1] = -psi2[_IndexTool(i, j, k)][1];

		imul(tmpres1, tmpl1, psi1[_IndexTool(ixp, j, k)]);
		imul(tmpres2, tmpl2, psi2[_IndexTool(ixp, j, k)]);

		vx[0] = tmpres1[0] + tmpres2[0];
		vx[1] = tmpres1[1] + tmpres2[1];

		imul(tmpres1, tmpl1, psi1[_IndexTool(i, iyp, k)]);
		imul(tmpres2, tmpl2, psi2[_IndexTool(i, iyp, k)]);

		vy[0] = tmpres1[0] + tmpres2[0];
		vy[1] = tmpres1[1] + tmpres2[1];

		imul(tmpres1, tmpl1, psi1[_IndexTool(i, j, izp)]);
		imul(tmpres2, tmpl2, psi2[_IndexTool(i, j, izp)]);

		vz[0] = tmpres1[0] + tmpres2[0];
		vz[1] = tmpres1[1] + tmpres2[1];

		_velX[_IndexTool(i, j, k)] = iangle(vx) * localhbar;
		_velY[_IndexTool(i, j, k)] = iangle(vy) * localhbar;
		_velZ[_IndexTool(i, j, k)] = iangle(vz) * localhbar;
	}GRID_LOOP_END
}

void ISF::divergence()
{
	GRID_LOOP_BIGIN{
		double dx2 = _Torus.dx * _Torus.dx;
		double dy2 = _Torus.dy * _Torus.dy;
		double dz2 = _Torus.dz * _Torus.dz;

		int ixm = (i - 1 + _Torus.gridx) % _Torus.gridx;
		int iym = (j - 1 + _Torus.gridy) % _Torus.gridy;
		int izm = (k - 1 + _Torus.gridz) % _Torus.gridz;

		_DIV[_IndexTool(i, j, k)][0] = (_velX[_IndexTool(i, j, k)] - _velX[_IndexTool(ixm, j, k)]) / dx2
			+ (_velY[_IndexTool(i, j, k)] - _velY[_IndexTool(i, iym, k)]) / dy2
			+ (_velZ[_IndexTool(i, j, k)] - _velZ[_IndexTool(i, j, izm)]) / dz2;
		//std::cout << _DIV[_IndexTool(i, j, k)][0] << std::endl;
		_DIV[_IndexTool(i, j, k)][1] = 0.0;
	}GRID_LOOP_END
}

void ISF::poissonSolve()
{
	fftw_complex* FFTOUT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) *BUFFER_SIZE);
	fftw_plan p = fftw_plan_dft_3d(RES_X, RES_Y, RES_Z, _DIV, FFTOUT, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);

	GRID_LOOP_BIGIN{
		imul(FFTOUT[_IndexTool(i, j, k)], PoissonParameters[_IndexTool(i, j, k)]);
	}GRID_LOOP_END

	fftw_plan p2 = fftw_plan_dft_3d(RES_X, RES_Y, RES_Z, FFTOUT, PoissonSolutions, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p2);
	fftw_destroy_plan(p2);

	GRID_LOOP_BIGIN{
		PoissonSolutions[_IndexTool(i, j, k)][0] /= BUFFER_SIZE;
		PoissonSolutions[_IndexTool(i, j, k)][1] /= BUFFER_SIZE;
	}GRID_LOOP_END
	if (FFTOUT != NULL) fftw_free(FFTOUT);
}

void ISF::initPoisson()
{
	GRID_LOOP_BIGIN{
		double sx = sin(M_PI * i / _Torus.gridx) / _Torus.dx;
		double sy = sin(M_PI * j / _Torus.gridy) / _Torus.dy;
		double sz = sin(M_PI * k / _Torus.gridz) / _Torus.dz;
		double size = sx * sx + sy * sy + sz * sz;
		PoissonParameters[_IndexTool(i, j, k)] = -0.25 / size;
	}GRID_LOOP_END

		PoissonParameters[0] = 0;
}


void ISF::GaugeTransform()
{
	GRID_LOOP_BIGIN{
		fftw_complex negi = {0.0, -1.0};
		fftw_complex tmp, tmp2;
		imul(tmp, negi, PoissonSolutions[_IndexTool(i, j, k)]);
		iexp(tmp2, tmp);
		imul(psi1[_IndexTool(i, j, k)], psi1[_IndexTool(i, j, k)], tmp2);
		imul(psi2[_IndexTool(i, j, k)], psi2[_IndexTool(i, j, k)], tmp2);
	}GRID_LOOP_END
}

void ISF::Setup()
{
	itrNum = 0;

	BuildISF();
	InitSmoke();

	initPoisson();

	for (int i = 0; i < 10; i++)
	{
		constrainVel(0.0);
	}
}


void ISF::SchrodingerFlowIntegrator()
{
	fftw_complex* FFTOUT1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * BUFFER_SIZE);
	fftw_plan p1 = fftw_plan_dft_3d(RES_X, RES_Y, RES_Z, psi1, FFTOUT1, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p1);
	fftw_destroy_plan(p1);

	fftw_complex* FFTOUT2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * BUFFER_SIZE);
	fftw_plan p2 = fftw_plan_dft_3d(RES_X, RES_Y, RES_Z, psi2, FFTOUT2, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p2);
	fftw_destroy_plan(p2);

	fftShift(FFTOUT1);
	fftShift(FFTOUT2);

	GRID_LOOP_BIGIN{
		imul(FFTOUT1[_IndexTool(i, j, k)], FFTOUT1[_IndexTool(i, j, k)], SchrodingerMask[_IndexTool(i, j, k)]);
	    imul(FFTOUT2[_IndexTool(i, j, k)], FFTOUT2[_IndexTool(i, j, k)], SchrodingerMask[_IndexTool(i, j, k)]);
	}GRID_LOOP_END

	fftShift(FFTOUT1);
	fftShift(FFTOUT2);

	fftw_plan p3 = fftw_plan_dft_3d(RES_X, RES_Y, RES_Z, FFTOUT1, psi1, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p3);
	fftw_destroy_plan(p3);

	fftw_plan p4 = fftw_plan_dft_3d(RES_X, RES_Y, RES_Z, FFTOUT2, psi2, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p4);
	fftw_destroy_plan(p4);

	if (FFTOUT1 != NULL) fftw_free(FFTOUT1);
	if (FFTOUT2 != NULL) fftw_free(FFTOUT2);

	psiNormalize();
}

void ISF::StaggeredAdvect()
{
	for (int i = 0; i < PARTICLE_NUM; i++)
	{
		if (particleDie[i])
			continue;
		double k1x, k1y, k1z;
		double k2x, k2y, k2z;
		double k3x, k3y, k3z;
		double k4x, k4y, k4z;


		double x = particleX[i];
		double y = particleY[i];
		double z = particleZ[i];

		int ix = particleGridX[i];
		int iy = particleGridY[i];
		int iz = particleGridZ[i];
		_density[_IndexTool(ix, iy, iz)] -= 1;

		StaggeredVelocity( x, y, z, k1x, k1y, k1z );
		StaggeredVelocity( x + DT * k1x / 2, y + DT * k1y / 2, z + DT * k1z / 2, k2x, k2y, k2z );
		StaggeredVelocity( x + DT * k2x / 2, y + DT * k2y / 2, z + DT * k2z / 2, k3x, k3y, k3z );
		StaggeredVelocity( x + DT * k3x, y + DT * k3y, z+ DT * k3z, k4x, k4y, k4z );
		particleX[i] += DT / 6 * ( k1x + 2 * k2x + 2 * k3x + k4x );
		particleY[i] += DT / 6 * ( k1y + 2 * k2y + 2 * k3y + k4y );
		particleZ[i] += DT / 6 * ( k1z + 2 * k2z + 2 * k3z + k4z );

		particleX[i] += _Torus.sizex;
		particleX[i] = fmod(particleX[i], (double)_Torus.sizex);
		particleY[i] += _Torus.sizey;
		particleY[i] = fmod(particleY[i], (double)_Torus.sizey);
		particleZ[i] += _Torus.sizez;
		particleZ[i] = fmod(particleZ[i], (double)_Torus.sizez);

		/*if (x[i] > _Torus.sizex || y[i] > _Torus.sizey || z[i] > _Torus.sizez)
		{
			particleDie[i] = true;
			continue;
		}*/


		ix = particleX[i] / _Torus.dx;
		ix += _Torus.gridx;
		ix %= _Torus.gridx;
		iy = particleY[i] / _Torus.dy;
		iy += _Torus.gridy;
		iy %= _Torus.gridy;
		iz = particleZ[i] / _Torus.dz;
		iz += _Torus.gridz;
		iz %= _Torus.gridz;

		particleGridX[i] = ix;
		particleGridY[i] = iy;
		particleGridZ[i] = iz;
		//std::cout << x[i] << std::endl;
		_density[_IndexTool(ix, iy, iz)] += 1;
	}
}

double fma_lerp(double v0, double v1, double t)
{
	return fma(t, v1, fma(-t, v0, v0));
}

void ISF::StaggeredVelocity(double ox, double oy, double oz, double& kx, double& ky, double& kz)
{
	ox = fmod(ox + _Torus.sizex, (double)_Torus.sizex);
	oy = fmod(oy + _Torus.sizey, (double)_Torus.sizey);
	oz = fmod(oz + _Torus.sizez, (double)_Torus.sizez);

	int ix = floor(ox/ _Torus.dx);
	int iy = floor(oy/ _Torus.dy);
	int iz = floor(oz/ _Torus.dz);

	int ixp = (ix + 1 + _Torus.gridx) % _Torus.gridx;
	int iyp = (iy + 1 + _Torus.gridy) % _Torus.gridy;
	int izp = (iz + 1 + _Torus.gridz) % _Torus.gridz;

	int ind0 = _IndexTool(ix, iy, iz);
	int indxp = _IndexTool(ixp, iy, iz);
	int indyp = _IndexTool(ix, iyp, iz);
	int indzp = _IndexTool(ix, iy, izp);
	int indxpyp = _IndexTool(ixp, iyp, iz);
	int indypzp = _IndexTool(ix, iyp, izp);
	int indxpzp = _IndexTool(ixp, iy, izp);

	double wx = ox - ix * _Torus.dx;
	wx /= _Torus.dx;
	double wy = oy - iy * _Torus.dy;
	wy /= _Torus.dy;
	double wz = oz - iz * _Torus.dz;
	wz /= _Torus.dz;

	kx = fma_lerp(fma_lerp(_velX[ind0], _velX[indyp], wy),
		fma_lerp(_velX[indzp], _velX[indypzp], wy), wz);
	ky = fma_lerp(fma_lerp(_velY[ind0], _velY[indzp], wz),
		fma_lerp(_velY[indxp], _velY[indxpzp], wz), wx);
	kz = fma_lerp(fma_lerp(_velZ[ind0], _velZ[indxp], wx),
		fma_lerp(_velZ[indyp], _velZ[indxpyp], wx), wy);
}



void ISF::SimulateStep()
{
	itrNum++;
	float totalTime = itrNum * DT;

	SchrodingerFlowIntegrator();
    pressureProject();

	constrainVel(totalTime);
	velocityOneForm(HBAR / _Torus.dx);
	StaggeredAdvect();
	ToFile(itrNum);
}


void ISF::Show()
{
	_renderer->Render();
}


void ISF::fftShift(fftw_complex* complexArray)
{
	for ( int idx = 0; idx < BUFFER_SIZE / 2; idx++ ) 
	{
		int x = idx / RES_Y / RES_Z;
		int mod = idx % ( RES_Y * RES_Z );
		int y = mod / RES_Z;
		int z = mod % RES_Z;

		fftw_complex temp;
		temp[0] = complexArray[idx][0];
		temp[1] = complexArray[idx][1];

		x += RES_X / 2;
		y += RES_Y / 2;
		z += RES_Z / 2;

		x %= RES_X;
		y %= RES_Y;
		z %= RES_Z;

		complexArray[idx][0] = complexArray[_IndexTool(x, y, z)][0];
		complexArray[idx][1] = complexArray[_IndexTool(x, y, z)][1];
		complexArray[_IndexTool(x, y, z)][0] = temp[0];
		complexArray[_IndexTool(x, y, z)][1] = temp[1];		
	}
}





void ISF::MouseButton(GLFWwindow* window, int button, int action, int mods)
{
	double mouseX, mouseY;
	if (button == GLFW_MOUSE_BUTTON_LEFT) 
	{
		glfwGetCursorPos(window, &mouseX, &mouseY);
		if (action == GLFW_PRESS) 
		{
			_isLeftKeyPressed = true;
			_arcball.StartRotation(mouseX, mouseY);
		}
		else if (action == GLFW_RELEASE) 
		{
			_isLeftKeyPressed = false;
			_arcball.StopRotation();
		}
	}
}

void ISF::MouseMotion(GLFWwindow* window, double nx, double ny)
{
	if (_isLeftKeyPressed) 
	{
		_arcball.UpdateRotation(nx, ny);
	}
}

void ISF::Resize(GLFWwindow* windowHandle, int x, int y)
{
	_arcball.SetWidthHeight(x, y);
}

void ISF::MouseScroll(GLFWwindow* window, double nx, double ny)
{
	_arcball.StartZooming(0, 0);
	_arcball.UpdateZooming(-ny, nx);
	_arcball.StopZooming();
}

void ISF::Reset()
{
	_arcball.SetWidthHeight(_winX, _winY);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0, 0, -5);
	glRotatef(0, 1, 0, 0);
	glRotatef(0, 0, 1, 0);
}

void ISF::RegisterParentWindow(GLFWwindow* windowHandle)
{
	_windowHandle = windowHandle;
	glfwGetWindowSize(_windowHandle, &_winX, &_winY);
}

//export .vol files, adopt from Mitsuba
void ISF::ToFile( int numF )
{
	std::ifstream fin;
	fin.open("output/scene.xml");
	char xmlline[128];
	snprintf(xmlline, sizeof(xmlline), "			\<string name=\"filename\" value=\"density-%04i.vol\"/>", numF);
	int line = 0;
	char tmpData[1024] = { 0 };
	std::string xmlData = "";
	while (fin.getline(tmpData, sizeof(tmpData)))
	{
		if (line == 11)
		{
			xmlData += xmlline;
			xmlData += "\n";
		}
		else
		{
			xmlData += tmpData;
			xmlData += "\n";
		}
		line++;
	}
	fin.close();

	char outname[128];
	snprintf(outname, sizeof(outname), "output/density-%04i.xml", numF);
	std::ofstream fout;
	fout.open(outname);
	fout.flush();
	fout << xmlData;
	fout.close();


	char fname[128], cmd[128];
	snprintf(fname, sizeof(fname), "output/density-%04i.vol", numF);
	//snprintf(cmd, sizeof(fname), "gzip -f output/density-%04i.vol&", numF);
	std::cout << "    + Saving \"" << fname << "\"" << std::endl;
	std::ofstream os(fname);

	int xres = RES_X,
		yres = RES_Y,
		zres = RES_Z;

	float scale = 1.0f / std::max(std::max(xres, yres), zres);

	os.write("VOL", 3);
	char version = 3;
	os.write((char*)&version, sizeof(char));
	int value = 1;
	os.write((char*)&value, sizeof(int));
	os.write((char*)&xres, sizeof(int));
	os.write((char*)&yres, sizeof(int));
	os.write((char*)&zres, sizeof(int));
	value = 1;
	os.write((char*)&value, sizeof(int));

	float minX = -xres / 2.0f * scale;
	float minY = -yres / 2.0f * scale;
	float minZ = -zres / 2.0f * scale;
	float maxX = xres / 2.0f * scale;
	float maxY = yres / 2.0f * scale;
	float maxZ = zres / 2.0f * scale;

	os.write((char*)&minX, sizeof(float));
	os.write((char*)&minY, sizeof(float));
	os.write((char*)&minZ, sizeof(float));
	os.write((char*)&maxX, sizeof(float));
	os.write((char*)&maxY, sizeof(float));
	os.write((char*)&maxZ, sizeof(float));

	for (int i = 0; i < RES_Z; i++)
		for (int j = 0; j < RES_Y; j++)
			for (int k = 0; k < RES_X; k++)
			{
				int idx = k * (RES_Y * RES_Z) + j * RES_Z + i;
				float value = _density[idx];
				os.write((char*)&value, sizeof(float));
			}
		/*int x = i / RES_Y / RES_Z;
        int mod = i % (RES_Y * RES_Z);
        int y = mod / RES_Z;
        int z = mod % RES_Z;*/
	os.close();

	//system(cmd);
}