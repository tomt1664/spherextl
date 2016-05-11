/* ------------------------------------------------------------
SphereXtl: Spherical polycrystal atomistic system constructor

Ang, Atom and Box struct definitions
GBead and Spheretype class definitions
--------------------------------------------------------------*/

// 3D rotation angle
struct Ang
{
    double x; // rotation about the x-axis
    double y; // y-axis
    double z; // z-axis
};

//atomic coordinates
struct Atom
{
    double x;
    double y;
    double z;
};

//simulation box dimensions
struct Box
{
    double x;
    double y;
    double z;
};

//graphite sphere: it has a position, radius and type
class Gbead
{
    double m_r;

public:
    Gbead(double xin, double yin, double zin, int stin, double rin);
    double radius() { return m_r; };
    double dist(Gbead bead, Box box);

    double x;
    double y; //coordinates
    double z;
    int st; // sphere type
    double ax;
    double ay;
    double az;
};

// class to define bead type properties
class Spheretype
{
    int m_num; //target number of spheres
    double m_rad; //the sphere radius
    double m_pad; //the padding: the excluded area at the sphere border
    double m_rep; //the excluded separation between spheres of the same type
    Ang m_ang; //the preferential lattice oritentation
    double m_spread; //the anisotropy factor
    double m_aparam; //a-axis lattice parameter
    double m_cparam; //c-axis lattice parameter

    std::vector<Atom> lattice; //the graphite lattice unit cell
    Atom orlat; //orthorhombic lattice vectors

public:
    Spheretype(int,double,double,double,double,double,double);  //constructor
    void setAng(double xa, double ya, double za); //set prefered angle
    int create(long numits,Box box,std::vector<Gbead>& beads,int stp,double rad); //create bead position and orientation in the simulation box
    long fill(std::vector<Gbead>& beads,std::vector<Atom>& atoms); //fill the bead with graphite lattice positions
    void rotate(std::vector<Gbead>& beads,std::vector<Atom>& atoms,long numat); //rotate the bead and set position
    void periodic(Box box,std::vector<Atom>& atoms,long numat); //apply the periodic boundary conditions

    double radius() { return m_rad; } //return radius
    int number() { return m_num; } //return target number
    double drand(void); // produce a uniform random number between 0 and 1
};
