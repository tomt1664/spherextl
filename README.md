# spherextl
Spherical crystallite polycrystal atomistic system constructor

C++ program to create atomistc system configurations for polycrystalline models of sintered granular materials


<img src="http://i66.tinypic.com/24opzic.png">


Spherical crystallites are randomly packed into a periodic simulation cell and then filled with a crystalline lattice. The lattice orientation can be completely random, or there can be a prefered orientation with a Gaussian spread determined by an anisotropy factor. Any number and spheres with any combination of radius and lattice vectors can be created, and there are options for padding and excuded separations. At the moment the lattice is set up to be the graphite structure, although any other crystal structure can be employed by simply modifying the unit cell definition (stored in the lattice array) which is defined in the Spheretype class constructor. 

The program requires an input file (`beads.dat`) containing the required parameters. 

For example:

```
4653
200.0  200.0  200.0
2
8  50.0 0.0 0.0 6.70 2.44 0.0
20 30.0 0.0 0.0 6.70 2.44 0.0
40 20.0 0.0 0.0 6.70 2.44 0.0
```

The first line specifies a seed for the random number generator. 
The second line specifies the x, y and z dimensions of the orthorhobic periodic cell. 
The third line specifies the number of types of spherical crystallites, and the number of following lines to define each one. 
Then follows a line for each type of sphere: the first number is the number of that type of sphere, the 2nd is the sphere radius, the 3rd the sphere padding, the 4th the excluded radius, the 5th the c lattice parameter for that sphere type, the 6th the a lattice parameter and the 7th the anisotropy factor. If this is non-zero, another 3 numbers are expected specifying the prefered orietation (as rotation angles around the x,y and z axes). 

The progam the output two files: `spheres.out` which lists the radius and coordinates of each placed sphere, and `structure.xyz` which is an xyz format file with all the atomic positions. 

The program is unit agnostic (angles are input in degrees).