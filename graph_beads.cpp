#include <iostream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <graph_beads.h>

const double pi = 3.1415926535897;

// spheretype input constructor
Spheretype::Spheretype(int numb, double rad, double pad, double rep, double cax,double aax,double aspr)
    {
        if(rad < cax || rad < aax)
        {
            std::cerr << "Error: Radius smaller than lattice spacing" << std::endl;
            exit(1);
        }
        if(pad < 0.0 || rep < 0.0)
        {
            std::cerr << "Error: Padding/Repulsion cannot be less than zero" << std::endl;
            exit(1);
        }
        if(cax < 0.0 || aax < 0.0)
        {
            std::cerr << "Error: Lattice constants cannot be less than zero" << std::endl;
            exit(1);
        }
        if(numb < 0)
        {
            std::cerr << "Zero beads of radius " << rad << ": delete line from input" << std::endl;
            exit(1);
        }
        m_num = numb;
        m_rad = rad;
        m_pad = pad;
        m_rep = rep;
        m_aparam = aax;
        m_cparam = cax;
        m_spread = aspr;

        //set up the graphite unit cell for the spheretype
        Atom point;

        point.x = 0.0;   //atom 1
        point.y = 0.0;
        point.z = 0.0;
        lattice.push_back(point);

        point.x = m_aparam*0.5;   //atom 2
        point.y = m_aparam*0.28868;
        point.z = 0.0;
        lattice.push_back(point);

        point.x = m_aparam*0.5;   //atom 3
        point.y = m_aparam*0.86603;
        point.z = 0.0;
        lattice.push_back(point);

        point.x = 0.0;   //atom 4
        point.y = m_aparam*1.15470;
        point.z = 0.0;
        lattice.push_back(point);

        point.x = 0.0;   //atom 5
        point.y = 0.0;
        point.z = m_cparam*0.5;
        lattice.push_back(point);

        point.x = 0.0;   //atom 6
        point.y = m_aparam*0.57735;
        point.z = m_cparam*0.5;
        lattice.push_back(point);

        point.x = m_aparam*0.5;   //atom 7
        point.y = m_aparam*0.86603;
        point.z = m_cparam*0.5;
        lattice.push_back(point);

        point.x = m_aparam*0.5;   //atom 8
        point.y = m_aparam*1.44338;
        point.z = m_cparam*0.5;
        lattice.push_back(point);

        //set orthorhombic lattice vectors
        orlat.x = m_aparam;
        orlat.y = m_aparam*1.7321;
        orlat.z = m_cparam;
    }

double Spheretype::drand(void)
{
    //function to return a uniform deviate (between 0 and 1) based on rand()
    double udev = (rand()*1.0)/(RAND_MAX*1.0);
    return udev;
}

//setting anisotropy factors
void Spheretype::setAng(double xa, double ya, double za)
    {   // set anisotropy prefered orientation in radians
        m_ang.x = xa;
        m_ang.y = ya;
        m_ang.z = za;
    }

//create spheres in the simulation box
int Spheretype::create(long numits,Box box,std::vector<Gbead>& beads,int stp,double rad)
    {
        int tooclose = 0;

        double bx,by,bz;
        //generate random position in box
        bx = box.x*drand();
        by = box.y*drand();
        bz = box.z*drand();

        Gbead curb(bx,by,bz,stp,rad); //current placed bead

        //generate random angles in radians for each axis
        curb.ax = drand()*pi;
        curb.ay = drand()*pi;
        curb.az = drand()*pi;

        //max number of tries: numits
        for(long n = 0; n < numits; n++)
        {
            for(unsigned i = 0; i < beads.size(); i++)
            {
                if(beads[i].dist(curb,box) < (beads[i].radius() + rad))
                {
                    tooclose = 1;
//                    std::cout << "too close break " << beads[i].dist(curb,box) << " br " << beads[i].radius() << " r " << rad << std::endl;
                    break;
                }
            }

            if(tooclose == 0)
            {
                beads.push_back(curb);
//                std::cout << "n its " << n << std::endl;
                return 1;
            }

            //generate new random position
            bx = box.x*drand();
            by = box.y*drand();
            bz = box.z*drand();

            curb.x = bx;
            curb.y = by;
            curb.z = bz;

            tooclose = 0;
        }
        return 0;
    }

//member function to fill the placed beads with graphite lattice atoms
long Spheretype::fill(std::vector<Gbead>& beads,std::vector<Atom>& atoms)
    {
        Atom atm;
        long natm = 0;

        double r = beads.back().radius();
        double tdist;

        int xli = static_cast<int>(r/orlat.x) + 1; // lattice repetitions to fill sphere
        int yli = static_cast<int>(r/orlat.y) + 1;
        int zli = static_cast<int>(r/orlat.z) + 1;

        for(int ix = -xli; ix <= xli; ix++)
        {
            for(int iy = -yli; iy <= yli; iy++)
            {
                for(int iz = -zli; iz <= zli; iz++)
                {
                    for(int il = 0; il < 8; il++)
                    {
                        atm.x = ix*orlat.x + lattice[il].x;
                        atm.y = iy*orlat.y + lattice[il].y;
                        atm.z = iz*orlat.z + lattice[il].z;

                        //exclude atoms falling outside the radius minus padding
                        tdist = sqrt(atm.x*atm.x + atm.y*atm.y + atm.z*atm.z);

                        if(tdist < (m_rad - m_pad))
                        {
                            atoms.push_back(atm);
                            natm++;
                        }
                    }
                }
            }
        }
        return natm;
    }

//function to rotate each bead to the set orientation and then centre the bead on the set position
void Spheretype::rotate(std::vector<Gbead>& beads,std::vector<Atom>& atoms,long numat)
    {
        Atom cntr; //centre point of the sphere
        cntr.x = beads.back().x;
        cntr.y = beads.back().y;
        cntr.z = beads.back().z;

        double xr,yr,zr;

        double sinx = sin(beads.back().ax);
        double siny = sin(beads.back().ay);
        double sinz = sin(beads.back().az);

        double cosx = cos(beads.back().ax);
        double cosy = cos(beads.back().ay);
        double cosz = cos(beads.back().az);

       for(int i = 0; i < numat; i++)
       {
           //rotate about x
           yr = cosx*atoms[i].y + sinx*atoms[i].z;
           zr = cosx*atoms[i].z - sinx*atoms[i].y;
           atoms[i].y = yr;
           atoms[i].z = zr;

           //rotate about y
           xr = cosy*atoms[i].x - siny*atoms[i].z;
           zr = cosy*atoms[i].z + siny*atoms[i].x;
           atoms[i].x = xr;
           atoms[i].z = zr;

           //rotate about z
           xr = cosz*atoms[i].x + sinz*atoms[i].y;
           yr = cosz*atoms[i].y - sinz*atoms[i].x;
           atoms[i].x = xr;
           atoms[i].y = yr;

           //place bead at set point
           atoms[i].x = atoms[i].x + cntr.x;
           atoms[i].y = atoms[i].y + cntr.y;
           atoms[i].z = atoms[i].z + cntr.z;
       }
    }

//function to apply periodic cell conditions to the filled bead
void Spheretype::periodic(Box box,std::vector<Atom>& atoms,long numat)
{
    for(int i = 0; i < numat; i++)
    {
        if(atoms[i].x < 0.0) atoms[i].x = atoms[i].x + box.x;
        if(atoms[i].x > box.x) atoms[i].x = atoms[i].x - box.x;

        if(atoms[i].y < 0.0) atoms[i].y = atoms[i].y + box.y;
        if(atoms[i].y > box.y) atoms[i].y = atoms[i].y - box.y;

        if(atoms[i].z < 0.0) atoms[i].z = atoms[i].z + box.z;
        if(atoms[i].z > box.z) atoms[i].z = atoms[i].z - box.z;
    }
}

//constructor for gbead
Gbead::Gbead(double xin, double yin, double zin, int stin, double rin)
{
    x = xin;
    y = yin;
    z = zin;
    st = stin;
    m_r = rin;
}

//calcualte the distance between a bead and another
// including periodic images
// return the smallest distance found
double Gbead::dist(Gbead bead,Box box)
{
    double xd,yd,zd;
    double xp,yp,zp;
    double tmpdist,mindist;

    //initialise minium separation
    mindist = box.x + box.y + box.z;

    //loop over all periodic images
    for(int ip = -1; ip < 2 ; ip++)
    {
        for(int jp = -1; jp < 2 ; jp++)
        {
            for(int kp = -1; kp < 2 ; kp++)
            {
                xp = ip*box.x;
                yp = jp*box.y;
                zp = kp*box.z;

                xd = x - bead.x + xp;
                yd = y - bead.y + yp;
                zd = z - bead.z + zp;

                tmpdist = sqrt(xd*xd + yd*yd + zd*zd);
                if(tmpdist < mindist)
                {
                    mindist = tmpdist;
                }
            }
        }
    }
    return mindist;
}
