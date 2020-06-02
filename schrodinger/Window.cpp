#include "Window.h"
#include "core.h"

extern Window* _Window;
static void error_callback(int error, const char* description)
{
	fputs(description, stderr);
}

static void resize(GLFWwindow* window, int x, int y)
{
	_Window->Resize(window, x, y);
}

static void mousebutton(GLFWwindow* window, int button, int action, int mods)
{
	_Window->MouseButton(window, button, action, mods);
}

static void mousemotion(GLFWwindow* window, double x, double y)
{
	_Window->MouseMotion(window, x, y);
}

static void mousescroll(GLFWwindow* window, double x, double y)
{
	_Window->MouseScroll(window, x, y);
}

Window::Window(int argc, char** argv, const char* windowName)
{
	_windowName = windowName;

	_ts = 0;
	_frames = 0;
	_fps = 0;

	_winX = 1366;
	_winY = 768;

	_isf = new ISF;

	_prevCursorX = 0;
	_prevCursorY = 0;

	_isLeftKeyPressed = 0;

	if ( !glfwInit() ) 
	{
		std::cerr << "glfwInit() failed!" << std::endl;
		exit(-1);
	}

	_windowHandle = glfwCreateWindow(_winX, _winY, windowName, NULL, NULL);
	if (!_windowHandle) 
	{
		std::cerr << "Create Window failed" << std::endl;
		exit(-1);
	}
	glfwMakeContextCurrent(_windowHandle);
	glfwSetWindowPos(_windowHandle, 0, 0);
	glClearColor(0., 0., 0., 1.);

	glfwSetErrorCallback(error_callback);
	glfwSetMouseButtonCallback(_windowHandle, mousebutton);
	glfwSetScrollCallback(_windowHandle, mousescroll);
	glfwSetCursorPosCallback(_windowHandle, mousemotion);
	glfwSetWindowSizeCallback(_windowHandle, resize);

	if (glewInit() != GLEW_OK)
	{
		std::cerr << "glewInit() failed!" << std::endl;
		exit(-1);
	}
}

void Window::Init()
{
	_isf->RegisterParentWindow(_windowHandle);
	_isf->Reset();
	_isf->Setup();
	_isf->Show();
}

void Window::InitCamera()
{
	_camera = new Camera(_windowHandle);
}

void Window::BeginLoop()
{
	while (!glfwWindowShouldClose(_windowHandle))
	{
		_Window->Render();
		glfwSwapBuffers(_windowHandle);
		glfwPollEvents();
	}
}

Window::~Window() 
{
	glFinish();
	glfwDestroyWindow(_windowHandle);
	glfwTerminate();
}

void Window::Reset() 
{
	_camera->Reset();
	_isf->Reset();
}

void Window::FPS()
{
	_frames++;
	_ts += _timer.StopTimer();
	_timer.StartTimer();
	if (_ts > 1) 
	{
		_fps = _frames / _ts;
		_ts = 0;
		_frames = 0;
	}
	_title.str("");
	_title << _windowName;
	_title << "     FPS: ";
	_title << std::setprecision(4) << _fps;
	_title.width(2);
}

void Window::Render() 
{
	FPS();
	glfwSetWindowTitle(_windowHandle, _title.str().c_str());
	_camera->Reset();

	_isf->SimulateStep();
	_isf->Show();
}

void Window::Quit() {
	glFinish();
	glfwDestroyWindow(_windowHandle);
	exit(0);
}

void Window::Resize(GLFWwindow* window, int x, int y) 
{
	_winX = x;
	_winY = y;
	_isf->Resize(window, x, y);
}

void Window::MouseButton(GLFWwindow* window, int button, int action, int mods)
{
	_isf->MouseButton(window, button, action, mods);
}

void Window::MouseMotion(GLFWwindow* window, double nx, double ny)
{
	_isf->MouseMotion(window, nx, ny);
}

void Window::MouseScroll(GLFWwindow* window, double nx, double ny)
{
	_isf->MouseScroll(window, nx, ny);
}


