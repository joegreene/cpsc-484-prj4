Enhancements to Project 3
===========================

We tweaked the glDemo in project 3 to support rendering many cubes (instead of
one) onto a scene.  We also adjusted the camera angle so it points down
slightly at the center of the field.

Getting this right involved a lot of monkeying with the transformation matrix
to figure out how to do use tranformation matrices correctly.

Steps to compile/run:

	$ make
	$ ./gldemo

A window should appear with many bouncing cubes (50) of various sizes, each
bouncing at different rates.  You can use the left/right arrow keys to spin
around the field of cubes and the up/down arrow keys to pan in/out.
