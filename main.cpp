/*--------------------------------------------------------------
SphereXtl: Spherical polycrystal atomistic system constructor

C++ program to create atomistc system configurations for
polycrystalline models of sintered granular materials

T.Trevethan 2015
---------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <graph_beads.h>

using namespace std;

void readparams(Box &box,int &nsph,vector<Spheretype> &spheres,unsigned& seed);

int main()
{
    // the simulation cell size
    Box box;

    // number of bead types
    int nsph = 0;
    // total number of atoms
    long totat = 0;

    // the sphere types
    vector<Spheretype> spheres;

    // placed beads
    vector<Gbead> beads;

    //lattice atoms - temporary array for each bead
    vector<Atom> atoms;

    // the number of attempts to place a sphere
    long piter = 1000000;

    //seed for PRNG
    unsigned seed;

    // read input file
    readparams(box,nsph,spheres,seed);

    //initialise PRNG with seed
    srand(seed);

    // write input summary
    cout << "gBeads" << endl;
    cout << "t.trevethan" << endl;
    cout << "Cell dimensions: x = " << box.x << " y = " << box.y << " z = " << box.z << endl;
    cout << "Sphere types: " << nsph << endl;

    for(int i=0; i < nsph; i++)
    {
        cout << "    " << spheres[i].number() << " x " << spheres[i].radius() << " A " << endl;
    }

    // opening output files
    // xyz format file for the atomistic structure
    ofstream outconfg("structure.xyz");
    // list of placed spheres
    ofstream outsphr("spheres.out");

    if (!outconfg)
    {
        cerr << "Error: structure.xyz could not be opened for writing" << endl;
        exit(1);
    }
    if (!outsphr)
    {
        cerr << "Error: spheres.out could not be opened for writing" << endl;
        exit(1);
    }

    //loop over the sphere types
    for(int ist=0; ist < nsph; ist++)
    {
        //create a sphere in the periodic box for the number of each sphere
        for(int isn = 0; isn < spheres[ist].number(); isn++)
        {
            //create sphere
            int csuc;
            csuc = spheres[ist].create(piter,box,beads,ist,spheres[ist].radius());
            if(!csuc)
            {
                cout << "Sphere " << isn << " of radius " << spheres[ist].radius() << " cannot be placed "<< endl;
                cout << "Exiting ..." << endl;
                exit(1);
            }
            //fill sphere with lattice structure
            long numat;
            atoms.clear(); // clear temporary vector
            numat = spheres[ist].fill(beads,atoms);
            //keep runnning total of the number of atoms
            totat += numat;
            spheres[ist].rotate(beads,atoms,numat);
            spheres[ist].periodic(box,atoms,numat);
            //write atomic positions to output configuration file
//            cout << " filled sphere " << isn << " with " << numat << endl;
//            cout << "vec size " << atoms.size() << endl;
            for(int irc = 0; irc < numat; irc++)
            {
                outconfg << "C " << atoms[irc].x << " " << atoms[irc].y << " " << atoms[irc].z << endl;
            }

            outsphr << beads.back().radius() << " " << beads.back().x << " " << beads.back().y << " "<< beads.back().z << endl; //writing placed sphere positions to file

        }
        cout << "Created " << spheres[ist].number() << " spheres of radius " << spheres[ist].radius() << " Total atoms: " << totat << endl;
    }

    //print finalise
    cout << "Completed" << endl;

    return 0;
}

void readparams(Box &box,int &nsph,vector<Spheretype> &spheres,unsigned& seed)
{
    // function to read input parameters from file
    ifstream inf("beads.dat");
    if (!inf)
    {
        cerr << "Error: Input file beads.dat could not be opened" << endl;
        exit(1);
    }

    inf >> seed;
    inf >> box.x >> box.y >> box.z;
    inf >> nsph;

    for(int i=0; i < nsph; i++)
    {
        int ntsp;
        double rsph, rpad, rrep, clat, alat, aspr;
        inf >> ntsp >> rsph >> rpad >> rrep >> clat >> alat >> aspr;
        spheres.push_back(Spheretype(ntsp,rsph,rpad,rrep,clat,alat,aspr));
        // read data for anisotropy factors if present
        if(aspr > 0.0001)
        {
            int ax,ay,az;
            inf >> ax;
            inf >> ay;
            inf >> az;
            spheres[i].setAng(ax,ay,az);
        }
    }
}
