#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include <tuple>
#include "Fluid.h"
 
#include <iostream>
typedef float Real;

// Size of square grid
int N = 80;
int SCALE = 10;

// GLFW variables
GLFWwindow *gWindow = nullptr;
int gWindowWidth = 2*N;
int gWindowHeight = 2*N;

// Fluid
Fluid* fluid;

// global options
bool gPause = true;
bool gSaveFile = false;
bool gShowGrid = false;
bool gShowVel = false;
int gSavedCnt = 0;
bool leftMouseButtonPressed = false;

bool gReset = false;
bool gReflect = false;

// Global variables to store last known mouse position
double lastXPos = 0.0;
double lastYPos = 0.0;

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
  // Initialize GLFW
  if(!glfwInit()) {
    std::cerr << "ERROR: Failed to init GLFW" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  // glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // only if requesting 3.0 or above
  // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_ANY_PROFILE); 
  glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

  // Create the window
  gWindowWidth = N;
  gWindowHeight = N;

  gWindow = glfwCreateWindow(
    N*SCALE, N*SCALE,
    "Advection-Reflection-Simulator", nullptr, nullptr);

  if(!gWindow) {
    std::cerr << "ERROR: Failed to open window" << std::endl;
    glfwTerminate();
    std::exit(EXIT_FAILURE);
  }

  // Load the OpenGL context in the GLFW window
  glfwMakeContextCurrent(gWindow);
  glfwGetFramebufferSize(gWindow, &gWindowWidth, &gWindowHeight);
  glfwSetWindowSizeCallback(gWindow, windowSizeCallback);

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

void setup(){
    fluid = new Fluid(N ,0.1, 0, 0);
    // int centerX = gWindowWidth / 2;  // Center X coordinate of the window
    // int centerY = gWindowHeight / 2; // Center Y coordinate of the window

    // int halfWidth = gWindowWidth / 2; // Half width of the window
    // int halfHeight = gWindowHeight / 2; // Half height of the window

    // int startX = static_cast<float>(centerX / gWindowWidth * N);
    // int startY = static_cast<float>((gWindowHeight - centerY) / gWindowHeight * N);

    // // Iterate through the first half of the screen and add density
    // for (int y = 0; y < gWindowHeight; ++y) {
    //     for (int x = 0; x < gWindowWidth; ++x) {
    //         // Ensure that the coordinates are within bounds
    //         if (x >= 0 && x < N/2 && y >= 0 && y < N/2) {
    //             fluid->AddDensity(x, y, 1);
    //         }
    //         if (x >= N/2 && x < N && y >= N/2 && y < N) {
    //             fluid->AddDensity(x, y, 1);
    //         }
    //     }
    // }      
}

void init()
{
  initGLFW();                  
  initOpenGL();
  setup(); // Create Fluid
}

void FluidCubeStepReflect(Fluid *cube)
{
    int N = cube->size;
    float visc = cube->visc;
    float diff = cube->diff;
    float dt = cube->dt;
    float *Vx = cube->Vx;
    float *Vy = cube->Vy;
    float *Vx0 = cube->Vx0;
    float *Vy0 = cube->Vy0;
    float *s = cube->density0;
    float *density = cube->density;

    // Create temporary arrays to save Vx and Vy before modifying them
    float *u12x = new float[N * N];
    float *u12y = new float[N * N];

    Fluid::buoyancy(Vy0, density, N);

    for (int i = 0; i < N*N; i++)
    {
        Vy[i] -= dt * Vy0[i];
    }

    Fluid::Diffuse(1, Vx0, Vx, visc, dt/2, 4, N);
    Fluid::Diffuse(2, Vy0, Vy, visc, dt/2, 4, N);

    Fluid::project(Vx0, Vy0, Vx, Vy, 4, N);
    
  // u~1/2 = A(u0, u0, 1/2dt)
  
    Fluid::advect(1, Vx, Vx0, Vx0, Vy0, dt/2, N);
    Fluid::advect(2, Vy, Vy0, Vx0, Vy0, dt/2, N);

    // to do SAVE VX AND VY
    std::copy(Vx, Vx + N * N, u12x);
    std::copy(Vy, Vy + N * N, u12y);
        
    // u1/2 = Pu~1/2
    Fluid::project(Vx, Vy, Vx0, Vy0, 4, N);

    // u^1/2 = 2*u1/2 - u~1/2
    for(int i = 0; i < N*N; i++){
      u12x[i] = 2 * Vx[i] - u12x[i];
      u12y[i] = 2 * Vy[i] - u12y[i];
    }
  
    // u~1 = A(u^1/2, u1/2, 1/2dt)
    Fluid::advect(1, u12x, Vx, Vx0, Vy0, dt/2, N);
    Fluid::advect(2, u12y, Vy, Vx0, Vy0, dt/2, N);

    // u1 = Pu~1
    Fluid::project(u12x, u12y, Vx0, Vy0, 4, N);

    Fluid::Diffuse(0, s, density, diff, dt/2, 4, N);
    Fluid::advect(0, density, s, u12x, u12y, dt/2, N);
    
    delete[] u12x;
    delete[] u12y;
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
    float *s = cube->density0;
    float *density = cube->density;

    Fluid::buoyancy(Vy0, density, N);

    for (int i = 0; i < N*N; i++)
    {
        Vy[i] -= dt/2 * Vy0[i];
    }

    Fluid::Diffuse(1, Vx0, Vx, visc, dt, 4, N);
    Fluid::Diffuse(2, Vy0, Vy, visc, dt, 4, N);

    Fluid::project(Vx0, Vy0, Vx, Vy, 4, N);

  // u~1/2 = A(u0, u0, 1/2dt)
    Fluid::advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
    Fluid::advect(2, Vy, Vy0, Vx0, Vy0, dt, N);
        
    Fluid::project(Vx, Vy, Vx0, Vy0, 4, N);

    Fluid::Diffuse(0, s, density, diff, dt, 4, N);
    Fluid::advect(0, density, s, Vx, Vy, dt, N);
}


void render()
{
    //Objstacles
    int lower_bound = (N - 1) / 2 - 2;
    int upper_bound = (N - 1) / 2 + 2;

    glClearColor(.4f, .4f, .4f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBegin(GL_QUADS);
    float squareSize = 1.0; // Set the size of each square

    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            // Define the vertices of the square
            float d = fluid->density[fluid->IX(i,j, N)];
            float vx = fluid->Vx[fluid->IX(i,j, N)];
            float vy = fluid->Vy[fluid->IX(i,j, N)];
            glColor3f(d, d, d);
            if(i == 0 || i == N-1 || j == 0 || j == N-1)
              glColor3f(.0, .0, 1.0);

            if(i >= lower_bound && i < upper_bound && j >= 10 && j <= 11)
              glColor3f(.0, .0, 1.0);

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

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
  if (action == GLFW_PRESS) {
      if (key == GLFW_KEY_M) {
          std::cout << "Key 'm' pressed!" << std::endl;
          
      } else if (key == GLFW_KEY_N) {
          // Perform action for 'n'

            for(int j = 0; j < N; ++j){
              if(fluid->Vy[(N-1)*(N-1)+j] > 0.1)
                std::cout << fluid->Vy[(N-1)*(N-1)+j] << " ";
            }


          std::cout << "Key 'n' pressed!" << std::endl;
      } else if (key == GLFW_KEY_C) {
          std::cout << "Key 'c' pressed!" << std::endl;
          if(gReflect)
            gReflect = false;
          else
            gReflect = true;

          std::cout << "The reflection is " << gReflect << std::endl;

      } else if (key == GLFW_KEY_R) {
          std::cout << "Key 'r' pressed!" << std::endl;
          gReset = true;
      }
  }
}



int main() {
    init();


    glfwSetMouseButtonCallback(gWindow, mouseButtonCallback);
    glfwSetKeyCallback(gWindow, keyCallback);
    while(!glfwWindowShouldClose(gWindow)) {
        // for(int i = N/2-25; i <= N/2 - 20; i++){
        //   for(int j = 3; j < 4; j++){
        //     fluid->AddDensity(i, j, 0.1);
        //     //fluid->AddVelocity(i, j, 0,  0.02)V;
        //   }
        //  }

        // for(int i = N/2+20; i <= N/2+25; i++){
        //   for(int j = 3; j < 4; j++){
        //     fluid->AddDensity(i, j, 0.1);
        //     //fluid->AddVelocity(i, j, 0,  0.02);
        //   }
        // }
        if (leftMouseButtonPressed) {
            double xpos, ypos;
            glfwGetCursorPos(gWindow, &xpos, &ypos);
        
            int x = static_cast<float>(xpos / gWindowWidth * N);
            int y = static_cast<float>((gWindowHeight - ypos) / gWindowHeight * N);
            fluid->AddDensity(x, y, 0.5);
            //fluid->AddVelocity(x, y, 0,  0.01);
        }

        // fluid->ApplyBuoyancyForce(fluid->Vx, fluid->Vy, fluid->density, 0.01
        // , 0.01, fluid->dt, N);
        //fluid->ApplyVortexForce(fluid->Vx, fluid->Vy, fluid->density, N-5, 5, 0.0001, fluid->dt, N);

        if(gReflect)
          FluidCubeStepReflect(fluid);
        else
          FluidCubeStep(fluid);

        for (int i = 0; i < N*N; i++){ fluid->Vx0[i] = 0; fluid->Vy0[i] = 0; }


        render();
        glfwSwapBuffers(gWindow);
        glfwPollEvents();

        if(gReset){
          delete fluid;
          setup();
          gReset = false;
        }

        static auto lastPrintTime = std::chrono::steady_clock::now();
        auto currentTime = std::chrono::steady_clock::now();
        auto elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(currentTime - lastPrintTime).count();

        if (elapsedSeconds >= 1) {
            float acum = 0;
            for (int i = 0; i < N*N; i++){
              acum += (fluid->Vx[i]*fluid->Vx[i] + fluid->Vy[i]*fluid->Vy[i])/2;
            }
            std::cout << acum << std::endl;
            lastPrintTime = currentTime;
        }
    }
    clear();
    delete fluid;
    std::cout << " > Quit" << std::endl;
    return EXIT_SUCCESS;
}