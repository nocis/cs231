#pragma once
#include "core.h"
#include "arcball.h"
#include "core.h"
#include "renderer.h"
#include "Torus.h"

class Renderer;

class ISF : public EventListener 
{
public:
	ISF();
	~ISF();
	void RegisterParentWindow(GLFWwindow* windowHandle);
	void Reset();

	virtual void MouseButton(GLFWwindow* window, int button, int action, int mods);
	virtual void MouseMotion(GLFWwindow* window, double nx, double ny);
	virtual void MouseScroll(GLFWwindow* window, double nx, double ny);
	virtual void Keyboard(GLFWwindow* window, int key, int scancode, int action, int mods) {}
	virtual void Resize(GLFWwindow* window, int x, int y);
	void ToFile(int numF);

	Torus _Torus;
	fftw_complex* SchrodingerMask, * psi1, * psi2, * _DIV, * PoissonSolutions;
	float jetVel[3];
	float kvec[3];
	float Omega;

	float* _density;
	float* _velX;
	float* _velY;
	float* _velZ;
	float* Phase;
	float* PoissonParameters;

	int itrNum;
	double particleX[PARTICLE_NUM];
	double particleY[PARTICLE_NUM];
	double particleZ[PARTICLE_NUM];
	int particleGridX[PARTICLE_NUM];
	int particleGridY[PARTICLE_NUM];
	int particleGridZ[PARTICLE_NUM];
	bool particleDie[PARTICLE_NUM];

	void SimulateStep();
	void Show();

	void fftShift(fftw_complex* complexArray);
	void Setup();
	void InitSmoke();
	void BuildISF();
	void constrainVel(float totaltime);
	void pressureProject();
	void velocityOneForm(float localhbar);
	void divergence();
	void poissonSolve();
	void initPoisson();
	void GaugeTransform();
	void psiNormalize();


	void SchrodingerFlowIntegrator();
	void StaggeredAdvect();
	void StaggeredVelocity(double ox, double oy, double oz, double& kx, double& ky, double& kz);


	Renderer* _renderer;
	bool _isRendering;

	GLfloat _rotX;
	GLfloat _rotY;
	GLfloat _depth;
	bool _isLeftKeyPressed;

	Arcball _arcball;

	GLFWwindow* _windowHandle;
	int _winX, _winY;
};