#include "core.h"
#include "Window.h"

Window* _Window = NULL;

int main(int argc, char **argv) 
{
	_Window = new Window(argc, argv, "ISF");
	ISF* _isf = new ISF;
	_Window->Init();
	_Window->InitCamera();
	_Window->BeginLoop();
	return 0;
}
