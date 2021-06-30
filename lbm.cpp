#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include "cmdparser.h"

using namespace std;

struct Parameters
{
    int n_x;
    int n_y;
    int timesteps;
    double u_in;
    double RE;
    int sphere_x;
    int sphere_y;
    int diameter;
    string vtk_file_name;
    int vtk_step;

};

//Cell Lookup map in order to more transparently access the positions on the cells.
map<string, int> cell_lookup = {{"C",0},{"N",1},{"S",2},{"W",3},{"E",4},{"NW",5},{"NE",6},{"SW",7},{"SE",8}};
map<string,vector<int>> c_q={{"C",vector<int>{0,0}},{"N",vector<int>{0,1}},{"S",vector<int>{0,-1}},{"W",vector<int>{1,0}},{"E",vector<int>{1,0}},{"NW",vector<int>{-1,1}},{"NE",vector<int>{1,1}},{"SW",vector<int>{-1,-1}},{"SE",vector<int>{1,-1}}};
map<string,double> w_q ={{"C",4.0/9.0},{"N",1.0/9.0},{"S",1.0/9.0},{"W",1.0/9.0},{"E",1.0/9.0},{"NW",1.0/36.0},{"NE",1.0/36.0},{"SW",1.0/36.0},{"SE",1.0/36.0}};
map<string,int> cellType={{"fluid",0},{"boundary",1},{"velocity boundary",2},{"density boundary",3}};

//Method definition for distance of two points 
double distance(double , double , double , double);

template<typename T> class Lattice
{
    int  n_cells_x, n_cells_y;

    vector<vector<vector<T>>> domain;

    vector<vector<int>> flag_field;

public:
    Lattice(int N_x, int N_y):  n_cells_x(N_x), n_cells_y(N_y)
    {
        domain = vector<vector<vector<T>>>(n_cells_x, vector<vector<T>>(n_cells_y, vector<T>(9)));
        flag_field = vector<vector<int>>(n_cells_x, vector<int>(n_cells_y, 0));
    }

    void initialize_Latice(Parameters *);

    void put(int, int, string, T);
    T get(int, int, string);
    int getX();
    int getY();
    int get_flag(int, int);
    T &operator()(const int, const int, string);
};

//@Patrik
template<typename T> void Lattice<T>::initialize_Latice(Parameters *param)
{
    //Calculate the "bounding box of the circle"
    double radius = param->diameter / 2.0;
    int bounding_upper_x = ceil(param->sphere_x + radius), bounding_upper_y = ceil(param->sphere_y + radius);
    int bounding_lower_x = floor(param->sphere_x - radius), bounding_lower_y = floor(param->sphere_y - radius);
    
    for (int x = 0; x < n_cells_x; x++) 
    {
        for (int y = 0; y < n_cells_y; y++) 
        {
            //Cells on the upper and lower bound
            if(y == 0 || y == n_cells_y-1) {
                //domain remains untouched
                //Flag is set to no-slip
                flag_field[x][y] = cellType["boundary"];
            }
            //Cells which are on the left bound
            else if (x == 0) 
            {
                //domain remains untouched
                //Flag is set to velocity
                flag_field[x][y] = cellType["velocity boundary"];
            }
            //Cells which are on the right bound
            else if (x == n_cells_x-1)
            {
                 //domain remains untouched
                //Flag is set to desity
                flag_field[x][y] = cellType["density boundary"];
            }
            else
            {
                //All remaining cells are either regular fluid cells or part of the obstacle
                //TODO: include the cells where the radius equals the distance or not???? ---In examples the visulalisation without looks better
                //Cells that are inside the bounding box of the circle and  need checking
                if(!(x-1 > bounding_upper_x || y-1 > bounding_upper_y || x-1 < bounding_lower_x || y-1 < bounding_lower_y))
                {
                    if(distance((double) x-1, (double) y-1,  (double) param->sphere_x, (double) param->sphere_y) < radius){
                        //domain remains untouched
                        //Flag is set to no-slip
                        flag_field[x][y] = cellType["boundary"];
                        continue;
                    }
                }

                //The remaining  cells are all fluid cells
                for(auto q = cell_lookup.begin(); q != cell_lookup.end(); q++)
                {
                    domain[x][y][q->second] = w_q[q->first];
                }
                flag_field[x][y] = cellType["fluid"];
            }
        }
    } 
    return;
}

template<typename T> void Lattice<T>::put(int x, int y, string q, T value){
    domain[x][y][cell_lookup[q]] = value;
}

template<typename T> T Lattice<T>::get(int x, int y, string q)
{
    return domain[x][y][cell_lookup[q]];
}

template<typename T> int Lattice<T>::getX()
{
    return n_cells_x;
}

template<typename T> int Lattice<T>::getY()
{
    return n_cells_y;
}

template<typename T> int Lattice<T>::get_flag(int x, int y)
{
    return flag_field[x][y];
}

template<typename T> T &Lattice<T>::operator()(const int x, const int y, string q)
{
    return domain[x][y][cell_lookup[q]];
}

//TODO: Maybe write as Grid functions -> may be better option
//Method definitions for density and velocity
double density(Lattice<double> &, int, int );
void velocity(Lattice<double> &,int,int ,vector<double> &);


//----------------------------------------------------------------
//Write the given data to a VTK file with the specified filename and the specified number.
//----------------------------------------------------------------
void write_VTK_file(string filename, int number, Lattice<double> &outGrid,  Parameters *param)
{
    //"Calculate the proper filename"
    filename.append(to_string(number));
    filename.append(".vtk");

    ofstream outstream;
    outstream.open(filename);

    //Write the header
    outstream << "# vtk DataFile Version 4.0" << endl
              << "SiWiRVisFile" << endl
              << "ASCII" << endl
              << "DATASET STRUCTURED_POINTS" << endl
              << "DIMENSIONS " << param->n_x << " " << param->n_y << " 1" << endl
              << "ORIGIN 0 0 0" << endl
              << "SPACING 1 1 1" << endl
              << "POINT_DATA " << param->n_x * param->n_y << endl;
    outstream << endl;

    outstream << "SCALARS flags unsigned_int 1" << endl;
    outstream << "LOOKUP_TABLE default" << endl;
    //Write all flags
    for (int y = 1; y < outGrid.getY()-1; y++)
    {
       for(int x = 1; x < outGrid.getX()-1; x++)
       {
           outstream << outGrid.get_flag(x,y) << endl;
       }
    }
    outstream << endl;

    outstream << "SCALARS density double 1" << endl;
    outstream << "LOOKUP_TABLE default" << endl;
    //Write all densities
    for (int y = 1; y < outGrid.getY()-1; y++)
    {
       for(int x = 1; x < outGrid.getX()-1; x++)
       {
           outstream << density(outGrid, x, y) << endl;
       }
    }
    outstream << endl;
    
    outstream << "VECTORS velocity double" << endl;
    for (int y = 1; y < outGrid.getY()-1; y++)
    {
       for(int x = 1; x < outGrid.getX()-1; x++)
       {
           vector<double> u(2);
           velocity(outGrid, x, y, u);
           outstream << u[0] << " " << u[1] << " 0" << endl;
       }
    }
    return;
}


double distance(double x1, double y1, double x2, double y2)
{
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

//TODO @Lukas
double density(Lattice<double> &grid, int x, int y){
    return 1;
}

void velocity(Lattice<double> &grid,int x,int y,vector<double> &u){
    u[0] = 1;
    u[1] = 0;
    return;
}

void collideStep(Lattice<double> &grid){
    return;
}

//@Patrik
void borderUpdate(Lattice<double> &grid, double u_in){
    return;
}

void streamPullStep(Lattice<double> &inGrid,Lattice<double> &outGrid){
    return;
}

//@Lukas
void readInParameters(Parameters &p){
    return;
}

//@Lukas
//TODO: Read in the parameters and call a to be defined latice method which takes care of the simulation.
int main(int argc, char* argv[])
{
    cout << "Very quiet here..." << endl;
    Parameters p;
    p.n_x = 400;
    p.n_y = 80;
    p.diameter = 20;
    p.sphere_x = 100;
    p.sphere_y = 40;

    Lattice<double> l(p.n_x+2,p.n_y+2);
    l.initialize_Latice(&p);
    write_VTK_file("test", 1, l, &p);
    return 0;
}
