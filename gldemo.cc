//
// gldemo.cc
//
// OpenGL demo.
//
// CPSC 484, CSU Fullerton, Spring 2016, Prof. Kevin Wortman
// Project 3
//
// Name:
// 	* Joe Greene
// 	* Kyle Terrien
// 	* Adam Beck
//
// In case it ever matters, this file is hereby placed under the MIT
// License:
//
// Copyright (c) 2016, Kevin Wortman
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#ifdef __APPLE__ // for compiling on Mac OS X
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <unistd.h>

// custom header
#include "gmath.hh"

// Frame rate (Stupid way of limiting frame rate, but it works for us.)
#define FRAME_RATE 60.0

// Multiple cubes
#define NUM_CUBES 50
#define SIDE_LEN_MAX 4
#define SIDE_LEN_MIN 1
#define COORD_MAX 50
#define AMP_MAX 20
#define PERIOD_MIN 5
#define PERIOD_MAX 60

// Viewer
#define VIEWER_Y 10.0

// Whether or not to print messages to stderr when keys are pressed.
const bool DO_LOGGING = true;

// The viewer "orbits" the origin at discrete VIEWER_POSITIONS angles,
// e.g. 1 o'clock, 2 o'clock, etc. The viewer can also move closer or
// farther from the origin in discrete steps.
const int VIEWER_POSITIONS = 96,
  VIEWER_STEPS_FROM_ORIGIN = 24;

const double WINDOW_WIDTH = 640,
  WINDOW_HEIGHT = 480,
  FIELD_OF_VIEW_DEGREES = 90,
  ASPECT_RATIO = (WINDOW_WIDTH / WINDOW_HEIGHT),

  // how many world coordinates to move when the viewer moves
  // closer/farther from the origin
  VIEWER_STEP = 1.0;

// ASCII for escape
const char ESCAPE_CHARACTER = 27;

// Title of the GUI window.
const std::string WINDOW_NAME = "gldemo";

// Convenience typedefs and helper functions.

typedef gmath::Vector<double, 3> Color;
typedef gmath::Vector<double, 4> Vector4;

Color make_color(double r, double g, double b) {
  Color c;
  c[0] = r;
  c[1] = g;
  c[2] = b;
  return c;
}

Vector4 make_vector4(double x, double y, double z, double w) {
  Vector4 v;
  v[0] = x;
  v[1] = y;
  v[2] = z;
  v[3] = w;
  return v;
}

Vector4 make_point(double x, double y, double z) {
  auto result = make_vector4(x, y, z, 1.0);
  assert(result.is_homogeneous_point());
  return result;
}

Vector4 make_translation(double dx, double dy, double dz) {
  auto result = make_vector4(dx, dy, dz, 0.0);
  assert(result.is_homogeneous_translation());
  return result;
}

bool is_color_intensity(double x) {
  return ((x >= 0.0) && (x <= 1.0));
}

bool is_color(const Color& c) {
  return (is_color_intensity(c[0]) &&
          is_color_intensity(c[1]) &&
          is_color_intensity(c[2]));
}

Color web_color(uint_fast32_t hex) {
  assert(hex <= 0xFFFFFF);
  Color color;
  color[0] = (hex >> 16) / 255.0;
  color[1] = ((hex >> 8) & 0xFF) / 255.0;
  color[2] = (hex & 0xFF) / 255.0;
  assert(is_color(color));
  return color;
}

// Colors of the cube faces.

const Color CUBE_FRONT_COLOR = web_color(0x808000), // olive
  CUBE_TOP_COLOR = web_color(0xADFF2F), // green yellow
  CUBE_RIGHT_COLOR = web_color(0xFFD700), // gold
  CUBE_LEFT_COLOR = web_color(0xFFFFFF), // white
  CUBE_BOTTOM_COLOR = web_color(0x8A2BE2), // blue violet
  CUBE_BACK_COLOR = web_color(0xFF8C00); // dark orange

// Write str to stderr when DO_LOGGING is true.
void log(const char* str) {
  if (DO_LOGGING) {
    std::cerr << str << std::endl;
  }
}

// Class for a triangular mesh face. Each face is defined by three 3D
// points, and a color for the entire triangle.
class Face {
private:
  Color _color;
  Vector4 _v0, _v1, _v2;

public:
  Face(const Color& color,
       const Vector4 v0,
       const Vector4 v1,
       const Vector4 v2)
    : _color(color),
      _v0(v0),
      _v1(v1),
      _v2(v2) { }

  const Color& color() const { return _color; }

  const Vector4& v0() const { return _v0; }
  const Vector4& v1() const { return _v1; }
  const Vector4& v2() const { return _v2; }
};

// Data structure for an entire mesh. This is a very simple data
// structure that stores a vector of separate Face objects. It is not
// very space-efficient since coordinates are duplicated for multiple
// faces.
class Mesh {
private:
  std::vector<Face> _faces;

public:
  Mesh() { }

  void add_face(const Color& color,
                const Vector4 v0,
                const Vector4 v1,
                const Vector4 v2) {
    _faces.push_back(Face(color, v0, v1, v2));
  }

  int face_count() const {
    return _faces.size();
  }
  
  bool is_index(int i) const {
    return ((i >= 0) && (i < face_count()));
  }

  const Face& face(int index) const {
    assert(is_index(index));
    return _faces[index];
  }
};

// A floating cube is a cube, with different colors on each face, that
// oscillates vertically on a sine pattern.
class FloatingCube {
private:
  Mesh _mesh;
  Vector4 _origin;
  double _amplitude, _period;
  
public:
  FloatingCube(double side_length,
               const Vector4& origin,
               double amplitude,
               double period)
    : _origin(origin),
      _amplitude(amplitude),
      _period(period) {

    assert(side_length > 0.0);
    assert(origin.is_homogeneous_point());
    assert(amplitude >= 0.0);
    assert(period > 0.0);

    // Create the mesh geometry.

    double pos = side_length / 2.0, // offset of corners relative to the origin
           neg = -pos;

    // Enumerate the 8 combinations in the same order as counting in
    // binary.
    Vector4 left_bottom_far = make_point(neg, neg, neg),
      left_bottom_near =      make_point(neg, neg, pos),
      left_top_far =          make_point(neg, pos, neg),
      left_top_near =         make_point(neg, pos, pos),
      right_bottom_far =      make_point(pos, neg, neg),
      right_bottom_near =     make_point(pos, neg, pos),
      right_top_far =         make_point(pos, pos, neg),
      right_top_near =        make_point(pos, pos, pos);

    add_square(CUBE_FRONT_COLOR,
               left_bottom_near,
               left_top_near,
               right_top_near,
               right_bottom_near);
    add_square(CUBE_TOP_COLOR,
               left_top_near,
               left_top_far,
               right_top_far,
               right_top_near);
    add_square(CUBE_RIGHT_COLOR,
               right_bottom_near,
               right_top_near,
               right_top_far,
               right_bottom_far);
    add_square(CUBE_LEFT_COLOR,
               left_bottom_far,
               left_top_far,
               left_top_near,
               left_bottom_near);
    add_square(CUBE_BOTTOM_COLOR,
               left_bottom_far,
               left_bottom_near,
               right_bottom_near,
               right_bottom_far);
    add_square(CUBE_BACK_COLOR,
               right_bottom_far,
               right_top_far,
               left_top_far,
               left_bottom_far);
  }

  const Mesh& mesh() const { return _mesh; }

  Vector4 location(int time) const {
    
    double t = static_cast<double>(time) / _period,
           y = _amplitude * sin(t);
    if (y < 0) y *= -1;
    
    auto result = _origin;
    
    result[1] = y;

    return result;
  }

private:
  void add_square(const Color& color,
                  const Vector4& v0,
                  const Vector4& v1,
                  const Vector4& v2,
                  const Vector4& v3) {
    _mesh.add_face(color, v0, v1, v2);
    _mesh.add_face(color, v2, v3, v0);
  }  
};

// The Viewer is at a scalar distance from the origin, and a discrete
// "o'clock" position.
class Viewer {
private:
  double _distance;
  int _position;

public:
  Viewer(double distance,
         int position)
    : _distance(distance),
      _position(position) {

    assert(distance > 0.0);
    assert((position >= 0) && (position < VIEWER_POSITIONS));
  }

  double distance() const { return _distance; }
  int position() const { return _position; }

  void move_clockwise() {
    _position--;
    if (_position == -1) {
      _position = VIEWER_POSITIONS - 1;
    }
  }

  void move_counterclockwise() {
    _position++;
    if (_position == VIEWER_POSITIONS) {
      _position = 0;
    }
  }

  bool can_move_closer() {
    return (_distance > VIEWER_STEP);
  }
  
  void move_closer() {
    assert(can_move_closer());
    _distance -= VIEWER_STEP;
  }

  void move_farther() {
    _distance += VIEWER_STEP;
  }
  
  Vector4 location() const {
    double fraction_of_circle = static_cast<double>(_position) / static_cast<double>(VIEWER_POSITIONS),
      radians = fraction_of_circle * 2.0 * M_PI;
    return make_point(_distance * cos(radians),
                      VIEWER_Y,
                      _distance * sin(radians));
  }
};  

// A Scene is defined by one floating cube, one viewer, and a time
// value used to control the floating oscillation.
class Scene {
private:
  std::vector<FloatingCube> _cubes;
  Viewer _viewer;
  int _time;

public:
  Scene()
    : _viewer(VIEWER_STEP * VIEWER_STEPS_FROM_ORIGIN, 2),
      _time(0) {
    Vector4 origin;
    double side_len, x_coord, z_coord, amp, period;
    for(int i = 0; i < NUM_CUBES; i++) {
      side_len = rand() % (SIDE_LEN_MAX - SIDE_LEN_MIN) + SIDE_LEN_MIN;
      x_coord = rand() % COORD_MAX - COORD_MAX / 2.0;
      z_coord = rand() % COORD_MAX - COORD_MAX / 2.0;
      amp = rand() % AMP_MAX;
      period = rand() % (PERIOD_MAX - PERIOD_MIN) + PERIOD_MIN;
      origin = make_point(x_coord,
                          0.0,
                          z_coord);
      _cubes.push_back(FloatingCube(side_len, origin, amp, period));
    }
  }

  void left() {
    _viewer.move_counterclockwise();
  }
  
  void right() {
    _viewer.move_clockwise();
  }
  
  void up() {
    if (_viewer.can_move_closer()) {
      _viewer.move_closer();
    }
  }
  
  void down() {
    _viewer.move_farther();
  }

  void idle() {
    display();
    _time++;
  }
  
  void display() {
    // Render the cube.

    // Call glCLear() to clear the color buffers and depth
    // buffer.
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Call the following functions to set up a model-view
    // transformation matrix corresponding to the camera positioned at
    // _viewer.location() and pointing toward the (0, 0, 0) world
    // coordinate origin, with the standard up-vector (0, 1, 0).
    //
    // glMatrixMode
    // glLoadIdentity
    // gluLookAt

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(_viewer.location()[0], // eye x,y,z 
              _viewer.location()[1], 
              _viewer.location()[2],
              0, 0, 0,               // center x, y, z
              0, 1, 0);              // up x,y,z

    for (int j = 0; j < NUM_CUBES; ++j) {
      auto location = _cubes[j].location(_time);
      // Call glTranslated() so that the cube is moved to its
      // proper location. More specifically, the cube's coordinates,
      // which are currently in a cube-local object space, need to be
      // translated by the _cube.location(_time) vector, so that they
      // are in world space coordinates. Without this step, the cube
      // will stay stuck at the origin instead of bouncing.
      glTranslated(location[0], location[1], location[2]);

      // Call glBegin() to start drawing triangles.
      auto mesh = _cubes[j].mesh();
      glBegin(GL_TRIANGLES);
      for (int i = 0; i < mesh.face_count(); ++i) {
        auto face = mesh.face(i);
        
        // Call glColor3d to set the GL color to match
        // face.color().
        glColor3d(face.color()[0], face.color()[1], face.color()[2]);

        // Call glVertex3d() three times to register the three
        // vertices of this face. The vertices can be obtained wtih
        // face.v0(), face.v1(), and face.v2().
        glVertex3d(face.v0()[0], face.v0()[1], face.v0()[2]);
        glVertex3d(face.v1()[0], face.v1()[1], face.v1()[2]);
        glVertex3d(face.v2()[0], face.v2()[1], face.v2()[2]);
      }
      // Call glEnd() to stop drawing the triangles.
      glEnd();
      glTranslated(-location[0], -location[1], -location[2]);
    }

    // Swap the back buffer and display buffer.
    glutSwapBuffers();
  }

};

// Global scene object. Global variables are distasteful, but OpenGL
// forces us to use one since the keypress callback functions need to
// somehow influence the scene.
std::unique_ptr<Scene> global_scene;

// Callback, called by exit() when the program ends.
void exit_handler() {
  // Make certain to delete the global scene object.
  global_scene.reset();
}

// Callback, called when an ordinary (ASCII) key is pressed.
void ordinary_key_handler(unsigned char key, int x, int y) {
  switch (toupper(key)) {
  case 'Q':
  case ESCAPE_CHARACTER:
    log("quit");
    exit(0);
    break;
  }
}

// Callback, called when a special (e.g. arrow) key is pressed.
void special_key_handler(int key, int x, int y) {
  switch (key) {
  case GLUT_KEY_LEFT:
    log("left");
    global_scene->left();
    break;

  case GLUT_KEY_RIGHT:
    log("right");
    global_scene->right();
    break;

  case GLUT_KEY_UP:
    log("up");
    global_scene->up();
    break;

  case GLUT_KEY_DOWN:
    log("down");
    global_scene->down();
    break;
  }
}

// Callback, called when the GLUT loop is idle. Essentially this is
// called once per frame, and this is where we update our scene state
// and render one frame.
void idle_handler() {
  global_scene->idle();
  // Stupid way of limiting frame rate, but it works for us.
  // NOTE: usleep is a unix-specific function
  usleep(1000000.0 / FRAME_RATE);
}

// Callback, called whenever the window needs to be
// redrawn. Ordinarily this happens once when the window is first
// opened, then whenever the window re-appears after being hidden or
// minimized (if that ever happens).
void display_handler() {
  log("display");
  global_scene->display();
}

int main(int argc, char** argv) {
  // Set the seed
  srand(time(NULL));

  // Allocate the Scene object.
  global_scene.reset(new Scene());

  // Initialize OpenGL and GLUT.
  glutInit(&argc, argv);
  
  // Call the following functions to set up an OpenGL window
  // with RGB colors, double buffering, and depth buffering.
  //
  // glutInitDisplayMode
  // glutInitWindowSize
  // glutCreateWindow

  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH); // int mode as input
  glutInitWindowSize(WINDOW_WIDTH,WINDOW_HEIGHT);           // width by height
  glutCreateWindow(WINDOW_NAME.c_str());                    // program name (expects a char*)

  // Set up GLUT callbacks.
  atexit(exit_handler);
  glutKeyboardFunc(ordinary_key_handler);
  glutSpecialFunc(special_key_handler);
  glutIdleFunc(idle_handler);
  glutDisplayFunc(display_handler);

  // Call the following functions to set the clear color to
  // black; use flat shading; enable depth buffering; and set up
  // perspective projection based on FIELD_OF_VIEW_DEGREES and
  // ASPECT_RATIO.
  //
  // glClearColor
  // glShadeModel
  // glEnable
  // glDepthFunc
  // glMatrixMode
  // glLoadIdentity
  // glPerspective (is this supposed to be gluPerspective?)

  glClearColor(0,0,0,1);       // black color
  glShadeModel(GL_FLAT);       // flat shading
  glEnable(GL_DEPTH_TEST);     // enable depth buffering
  glDepthFunc(GL_LEQUAL);      // buffer when less or equal
  glMatrixMode(GL_PROJECTION); // perspective projection (I believe)
  glLoadIdentity();            // load the identity matrix

  // fovy, aspect, znear, zfar
  gluPerspective(FIELD_OF_VIEW_DEGREES, ASPECT_RATIO, 
                  0.01,   // misc small value 
                  1024);  // misc large value

  // Start the GLUT loop, which will run indefinitely until the
  // program exits.
  glutMainLoop();
  
  return 0;
}

// vi: et ts=2 sw=2
