#pragma once
#include "core.h"
#include "ISF.h"
#include "arcball.h"
#include "camera.h"
#include "timer.h"

class Window
{
public:
	Window(int argc, char** argv, const char* windowName);
	~Window();
	void Quit();
	void Init();
	void InitCamera();
	void FPS();
	void BeginLoop();
	void Reset();
	void Render();
	
	void Resize(GLFWwindow* window, int x, int y);
	void MouseButton(GLFWwindow* window, int btn, int action, int mods);
	void MouseMotion(GLFWwindow* window, double x, double y);
	void MouseScroll(GLFWwindow* window, double x, double y);

private:
	std::string _windowName;
	std::stringstream _title;
	GLFWwindow* _windowHandle;

	Camera* _camera;
	Timer _timer;
	double _ts;
	int _frames;
	double _fps;
	int _winX, _winY;

	bool _isLeftKeyPressed;
	double _prevCursorX, _prevCursorY;

	ISF * _isf;
};

