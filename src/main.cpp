#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <tuple>
#include "Fluid.h"
 
 
#include <iostream>
typedef float Real;
int N = 100;
GLFWwindow *gWindow = nullptr;
int gWindowWidth = N;
int gWindowHeight = N;
Fluid* fluid;

// global options
bool gPause = true;
bool gSaveFile = false;
bool gShowGrid = false;
bool gShowVel = false;
int gSavedCnt = 0;
bool leftMouseButtonPressed = false;

// Global variables to store last known mouse position
double lastXPos = 0.0;
double lastYPos = 0.0;

int SCALE = 10;
GLuint VBO, VAO;



// Executed each time the window is resized. Adjust the aspect ratio and the rendering viewport to the current window.
void windowSizeCallback(GLFWwindow *window, int width, int height)
{
  gWindowWidth = width;
  gWindowHeight = height;
  glViewport(0, 0, static_cast<GLint>(gWindowWidth), static_cast<GLint>(gWindowHeight));
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, N, 0, N, 0, 1);
}

void initGLFW()
{
  // Initialize GLFW, the library responsible for window management
  if(!glfwInit()) {
    std::cerr << "ERROR: Failed to init GLFW" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Before creating the window, set some option flags
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  // glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // only if requesting 3.0 or above
  // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_ANY_PROFILE); // for OpenGL below 3.2
  glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

  // Create the window
  gWindowWidth = N*SCALE;
  gWindowHeight = N*SCALE;
  gWindow = glfwCreateWindow(
    N*SCALE, N*SCALE,
    "Basic SPH Simulator", nullptr, nullptr);
  if(!gWindow) {
    std::cerr << "ERROR: Failed to open window" << std::endl;
    glfwTerminate();
    std::exit(EXIT_FAILURE);
  }

  // Load the OpenGL context in the GLFW window
  glfwMakeContextCurrent(gWindow);

  // not mandatory for all, but MacOS X
  glfwGetFramebufferSize(gWindow, &gWindowWidth, &gWindowHeight);

  // Connect the callbacks for interactive control
  glfwSetWindowSizeCallback(gWindow, windowSizeCallback);
  //glfwSetKeyCallback(gWindow, keyCallback);

  std::cout << "Window created: " <<
    gWindowWidth << ", " << gWindowHeight << std::endl;
}

void clear()
{
  glfwDestroyWindow(gWindow);
  glfwTerminate();
}

void exitOnCriticalError(const std::string &message)
{
  std::cerr << "> [Critical error]" << message << std::endl;
  std::cerr << "> [Clearing resources]" << std::endl;
  clear();
  std::cerr << "> [Exit]" << std::endl;
  std::exit(EXIT_FAILURE);
}

void initOpenGL()
{
  // Load extensions for modern OpenGL
  if(!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    exitOnCriticalError("[Failed to initialize OpenGL context]");

  glDisable(GL_CULL_FACE);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  glViewport(0, 0, static_cast<GLint>(gWindowWidth), static_cast<GLint>(gWindowHeight));
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, N, 0, N, 0, 1);
}

void init()
{
  //gSolver.initScene(48, 32, 16, 16);

  initGLFW();                   // Windowing system
  initOpenGL();
}




void setup(){
    fluid = new Fluid(N ,0.1, 0, 0);
}


void FluidCubeStep(Fluid *cube)
{
    int N = cube->size;
    float visc = cube->visc;
    float diff = cube->diff;
    float dt = cube->dt;
    float *Vx = cube->Vx;
    float *Vy = cube->Vy;
    float *Vx0 = cube->Vx0;
    float *Vy0 = cube->Vy0;
    float *s = cube->s;
    float *density = cube->density;

    Fluid::Diffuse(1, Vx0, Vx, visc, dt, 4, N);
    Fluid::Diffuse(2, Vy0, Vy, visc, dt, 4, N);

    Fluid::project(Vx0, Vy0, Vx, Vy, 4, N);

    Fluid::advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
    Fluid::advect(2, Vy, Vy0, Vx0, Vy0, dt, N);

    Fluid::project(Vx, Vy, Vx0, Vy0, 4, N);

    Fluid::Diffuse(0, s, density, diff, dt, 4, N);
    Fluid::advect(0, density, s, Vx, Vy, dt, N);
}


void render()
{
    glClearColor(.4f, .4f, .4f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBegin(GL_QUADS);
    float squareSize = 1.0; // Set the size of each square

    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            // Define the vertices of the square
            float d = fluid->density[fluid->IX(i,j, N)];
            glColor3f(1.0f - d, 1.0f - d, 1.0f - d);

            glVertex2f(static_cast<Real>(i), static_cast<Real>(j));
            glVertex2f(static_cast<Real>(i + squareSize), static_cast<Real>(j));
            glVertex2f(static_cast<Real>(i + squareSize), static_cast<Real>(j + squareSize));
            glVertex2f(static_cast<Real>(i), static_cast<Real>(j + squareSize));
        }
    }
    glEnd();

    if(gShowGrid) {
        glBegin(GL_LINES);
        for(int i=1; i<N; ++i) {
            glColor3f(0.3, 0.3, 0.3);
            glVertex2f(static_cast<Real>(i), 0.0);
            glColor3f(0.3, 0.3, 0.3);
            glVertex2f(static_cast<Real>(i), static_cast<Real>(N));
        }
        for(int j=1; j<N; ++j) {
            glColor3f(0.3, 0.3, 0.3);
            glVertex2f(0.0, static_cast<Real>(j));
            glColor3f(0.3, 0.3, 0.3);
            glVertex2f(static_cast<Real>(N), static_cast<Real>(j));
        }
        glEnd();
    }   
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            leftMouseButtonPressed = true;
        } else if (action == GLFW_RELEASE) {
            leftMouseButtonPressed = false;
        }
    }
}

int main() {
    setup();
    init();

    glfwSetMouseButtonCallback(gWindow, mouseButtonCallback);
    while(!glfwWindowShouldClose(gWindow)) {
        if (leftMouseButtonPressed) {
            double xpos, ypos;
            glfwGetCursorPos(gWindow, &xpos, &ypos);
            //fluid->printVelocity();
            // Convert screen coordinates to OpenGL coordinates
            int x = static_cast<float>(xpos / gWindowWidth * N);
            int y = static_cast<float>((gWindowHeight - ypos) / gWindowHeight * N);
            fluid->AddDensity(x, y, 1);
            fluid->AddVelocity(x, y, 0,  0.5);
            //fluid->AddVelocity(x, y, 0,  0.2);
            //std::cout << "Mouse position: (" << xpos - lastXPos << ", " << ypos - lastYPos << ")" << std::endl;
            lastXPos = xpos;
            lastYPos = ypos;

        
        }
        FluidCubeStep(fluid);
        render();
        glfwSwapBuffers(gWindow);
        glfwPollEvents();
    }
    clear();
    std::cout << " > Quit" << std::endl;
    return EXIT_SUCCESS;
}