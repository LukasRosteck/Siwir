#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <iomanip>
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
    double ny;
    double tau;
};


//Framework to identify the type and directions of a single cell
enum cellType {fluid, boundary, velocity_boundary, density_boundary};

enum direction {C, N, E, S, W, NW, NE, SW, SE};

map<direction, vector<int>> c_q = {{C, vector<int>{0, 0}}, {N, vector<int>{0, 1}}, {S, vector<int>{0, -1}}, {W, vector<int>{1, 0}}, {E, vector<int>{-1, 0}}, {NW, vector<int>{1, 1}}, {NE, vector<int>{-1, 1}}, {SW, vector<int>{1, -1}}, {SE, vector<int>{-1, -1}}};
map<direction, double> w_q = {{C, 4.0 / 9.0}, {N, 1.0 / 9.0}, {S, 1.0 / 9.0}, {W, 1.0 / 9.0}, {E, 1.0 / 9.0}, {NW, 1.0 / 36.0}, {NE, 1.0 / 36.0}, {SW, 1.0 / 36.0}, {SE, 1.0 / 36.0}};

//Method definition for distance of two points
double distance(double, double, double, double);

//Class representing teh Lattice Grid and all its neccessary information
class Lattice
{
    int n_cells_x, n_cells_y;

    vector<vector<vector<double>>> domain;

    vector<vector<int>> flag_field;

public:
    Lattice(int N_x, int N_y) : n_cells_x(N_x), n_cells_y(N_y)
    {
        domain = vector<vector<vector<double>>>(n_cells_x, vector<vector<double>>(n_cells_y, vector<double>(9)));
        flag_field = vector<vector<int>>(n_cells_x, vector<int>(n_cells_y, 0));
    }

    void initialize_Latice(Parameters *);

    void put(int, int, direction, double);
    double get(int, int, direction);
    int getX();
    int getY();
    int get_flag(int, int);
    double density(int, int);
    void velocity(int, int, vector<double> &);
    double &operator()(const int, const int, direction);
};


void Lattice::initialize_Latice(Parameters *param)
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
            if (y == 0 || y == n_cells_y - 1)
            {
                //domain remains untouched
                //Flag is set to no-slip
                flag_field[x][y] = boundary;
                for (int d = 0; d < 9; d++)
                {
                    domain[x][y][d] =.0100;
                }
            }
            //Cells which are on the left bound
            else if (x == 0)
            {
                //domain remains untouched
                //Flag is set to velocity
                flag_field[x][y] = velocity_boundary;
            }
            //Cells which are on the right bound
            else if (x == n_cells_x - 1)
            {
                //domain remains untouched
                //Flag is set to desity
                flag_field[x][y] = density_boundary;
                for(int d = 0; d < 9; d++) 
                {
                    domain[x][y][d] = w_q[(direction) d];
                }
            }
            else
            {
                //All remaining cells are either regular fluid cells or part of the obstacle
                //Cells that are inside the bounding box of the circle and  need checking
                if (!(x - 1 > bounding_upper_x || y - 1 > bounding_upper_y || x - 1 < bounding_lower_x || y - 1 < bounding_lower_y))
                {
                    if (distance((double) x - 0.5, (double) y - 0.5, (double) param->sphere_x, (double) param->sphere_y) <= radius)
                    {
                        //Flag is set to no-slip
                        flag_field[x][y] = boundary;
                        for (int i = 0; i < 9; i++)
                        {
                            domain[x][y][i] = .0100;
                        }
                        continue;
                    }
                }

                //The remaining  cells are all fluid cells
                 for(int d = 0; d < 9; d++) 
                {
                    domain[x][y][d] = w_q[(direction) d];
                }
                flag_field[x][y] = fluid;
            }
        }
    }
    return;
}


void Lattice::put(int x, int y, direction q, double value)
{
    domain[x][y][q] = value;
}

double Lattice::get(int x, int y, direction q)
{
    return domain[x][y][q];
}

int Lattice::getX()
{
    return n_cells_x;
}

int Lattice::getY()
{
    return n_cells_y;
}

int Lattice::get_flag(int x, int y)
{
    return flag_field[x][y];
}

double &Lattice::operator()(const int x, const int y, direction q)
{
    return domain[x][y][q];
}

double Lattice::density(int x, int y)
{
    double val = 0;
    for (int i = 0; i < 9; i++)
    {
        val += domain[x][y][i];
    }
    return val;
}

void Lattice::velocity(int x, int y, vector<double> &u)
{
    //Initialize the given velocity vector to 0
    u[0] = 0;
    u[1] = 0;

    //Calculate the velocity
    for (int d = 0; d < 9; d++)
    {
        u[0] += domain[x][y][d] * c_q[(direction) d][0];
        u[1] += domain[x][y][d] * c_q[(direction) d][1];
    }
}

//----------------------------------------------------------------
//Write the given data to a VTK file with the specified filename and the specified number.
//----------------------------------------------------------------
void write_VTK_file(string filename, int number, Lattice &outGrid, Parameters *param)
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
              << "DIMENSIONS " << param->n_x + 2 << " " << param->n_y + 2 << " 1" << endl
              << "ORIGIN 0 0 0" << endl
              << "SPACING 1 1 1" << endl
              << "POINT_DATA " << (param->n_x + 2) * (param->n_y + 2) << endl;
    outstream << endl;

    outstream << "SCALARS flags unsigned_int 1" << endl;
    outstream << "LOOKUP_TABLE default" << endl;
    //Write all flags
    for (int y = 0; y < outGrid.getY(); y++)
    {
        for (int x = 0; x < outGrid.getX(); x++)
        {
            outstream << outGrid.get_flag(x, y) << endl;
        }
    }
    outstream << endl;

    outstream << "SCALARS density double 1" << endl;
    outstream << "LOOKUP_TABLE default" << endl;
    //Write all densities
    for (int y = 0; y < outGrid.getY(); y++)
    {
        for (int x = 0; x < outGrid.getX(); x++)
        {
            outstream << outGrid.density(x, y) << endl;
        }
    }
    outstream << endl;

    outstream << "VECTORS velocity double" << endl;
    for (int y = 0; y < outGrid.getY(); y++)
    {
        for (int x = 0; x < outGrid.getX(); x++)
        {
            vector<double> u(2);
            outGrid.velocity(x, y, u);
            outstream << u[0] << " " << u[1] << " 0" << endl;
        }
    }
    return;
}

double distance(double x1, double y1, double x2, double y2)
{
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

void collideStep(Lattice &grid, Parameters *p)
{
    for (int x = 1; x < grid.getX(); x++)
    {
        for (int y = 1; y < grid.getY(); y++)
        {
            if (grid.get_flag(x, y) == fluid)
            {
                double density = grid.density(x, y);
                vector<double> u(2);
                grid.velocity(x, y, u);
                for (int d = 0; d < 9; d++)
                {
                    double vectorProduct = u[0] * c_q[(direction) d][0] + u[1] * c_q[(direction) d][1];

                    double uSquared = u[0] * u[0] + u[1] * u[1];
                    double f_eq = w_q[(direction) d] * (density + 3.0 * vectorProduct + (9.0 / 2.0) * vectorProduct * vectorProduct - (3.0 / 2.0) * uSquared);
                    double oldVal = grid.get(x, y, (direction) d);
                    double newVal = oldVal - (1.0 / p->tau) * (oldVal - f_eq);
                    grid.put(x, y, (direction) d, newVal);
                }
            }
        }
    }
    return;
}

//Method to update the helping border cells and the obstacle cells of the domain
void borderUpdate(Lattice &grid, double u_in, Parameters *param)
{
    //Step 1: north and south boundary
    for (int x = 0; x < grid.getX(); x++)
    {
        //Seperate Corner points
        if (x == 0)
        {
            //Reflect the values
            grid.put(x, 0, NE, grid(x + 1, 1, SW));

            grid.put(x, grid.getY() - 1, SE, grid(x + 1, grid.getY() - 2, NW));
        }
        else if (x == grid.getX() - 1)
        {
            //Reflect the values
            grid.put(x, 0, NW, grid(x - 1, 1, SE));

            grid.put(x, grid.getY() - 1, SW, grid(x - 1, grid.getY() - 2, NE));
        }
        else
        {
            //Reflect the values of the fluid cells.
            grid.put(x, 0, N, grid(x, 1, S));
            grid.put(x, 0, NW, grid(x - 1, 1, SE));
            grid.put(x, 0, NE, grid(x + 1, 1, SW));

            grid.put(x, grid.getY() - 1, S, grid(x, grid.getY() - 2, N));
            grid.put(x, grid.getY() - 1, SW, grid(x - 1, grid.getY() - 2, NE));
            grid.put(x, grid.getY() - 1, SE, grid(x + 1, grid.getY() - 2, NW));
        }
    }

    //Step 2: Iterate over domain, find the obstacle cells, and update them
    double radius = param->diameter / 2.0;
    int bounding_upper_x = ceil(param->sphere_x + radius), bounding_upper_y = ceil(param->sphere_y + radius);
    int bounding_lower_x = floor(param->sphere_x - radius), bounding_lower_y = floor(param->sphere_y - radius);
     for (int x = bounding_lower_x-1; x < bounding_upper_x+1; x++)
    {
        for (int y = bounding_lower_y-1; y < bounding_upper_y+1; y++)
        {
            //obstacle cell found
            if (grid.get_flag(x, y) == boundary)
            {
                grid.put(x, y, N, grid(x, y + 1, S));
                grid.put(x, y, E, grid(x + 1, y, W));
                grid.put(x, y, S, grid(x, y - 1, N));
                grid.put(x, y, W, grid(x - 1, y, E));

                grid.put(x, y, NW, grid(x - 1, y + 1, SE));
                grid.put(x, y, NE, grid(x + 1, y + 1, SW));
                grid.put(x, y, SE, grid(x + 1, y - 1, NW));
                grid.put(x, y, SW, grid(x - 1, y - 1, NE));
            }
        }
    }

    //Step 3: west (velocity) and east boundary
    for (int y = 1; y < grid.getY() - 1; y++)
    {
        grid.put(0, y, E, grid(1, y, W) + 6 * w_q[E] * u_in);
        grid.put(0, y, NE, grid(1, y + 1, SW) + 6 * w_q[NE] * u_in);
        grid.put(0, y, SE, grid(1, y - 1, NW) + 6 * w_q[SE] * u_in);

        vector<double> u(2);
        grid.velocity(grid.getX() - 1, y, u);

        double val = -grid(grid.getX() - 1, y, W) + 2 * w_q[W] * (1 + (9.0 / 2.0) * (c_q[W][0] * u[0] + c_q[W][1] * u[1]) * (c_q[W][0] * u[0] + c_q[W][1] * u[1]) - (3.0 / 2.0) * (u[0] * u[0] + u[1] * u[1]));
        grid.put(grid.getX() - 1, y, W, val);

        val = -grid(grid.getX() - 1, y, NW) + 2 * w_q[NW] * (1 + (9.0 / 2.0) * (c_q[NW][0] * u[0] + c_q[NW][1] * u[1]) * (c_q[NW][0] * u[0] + c_q[NW][1] * u[1]) - (3.0 / 2.0) * (u[0] * u[0] + u[1] * u[1]));
        grid.put(grid.getX() - 1, y, NW, val);

        val = -grid(grid.getX() - 1, y, SW) + 2 * w_q[SW] * (1 + (9.0 / 2.0) * (c_q[SW][0] * u[0] + c_q[SW][1] * u[1]) * (c_q[SW][0] * u[0] + c_q[SW][1] * u[1]) - (3.0 / 2.0) * (u[0] * u[0] + u[1] * u[1]));
        grid.put(grid.getX() - 1, y, SW, val);
        
    }

    return;
}

void streamPullStep(Lattice &inGrid, Lattice &outGrid)
{
    for (int x = 1; x < inGrid.getX() - 1; x++)
    {
        for (int y = 1; y < inGrid.getY() - 1; y++)
        {
            //Only stream the fluid cells
            if (inGrid.get_flag(x, y) == fluid)
            {
                //Pull the values frome the correct positions of the Grid
                outGrid.put(x, y, C, inGrid(x, y, C));
                outGrid.put(x, y, N, inGrid(x, y - 1, N));
                outGrid.put(x, y, E, inGrid(x - 1, y, E));
                outGrid.put(x, y, S, inGrid(x, y + 1, S));
                outGrid.put(x, y, W, inGrid(x + 1, y, W));
                outGrid.put(x, y, NW, inGrid(x + 1, y - 1, NW));
                outGrid.put(x, y, NE, inGrid(x - 1, y - 1, NE));
                outGrid.put(x, y, SE, inGrid(x - 1, y + 1, SE));
                outGrid.put(x, y, SW, inGrid(x + 1, y + 1, SW));
            }
        }
    }
    return;
}


void readInParameters(Parameters &p, string fileName)
{
    ifstream file;
    file.open(fileName, ios::in);
    string s;
    file >> s >> p.n_x;
    file >> s >> p.n_y;
    file >> s >> p.timesteps;
    file >> s >> p.u_in;
    file >> s >> p.RE;
    file >> s >> p.sphere_x;
    file >> s >> p.sphere_y;
    file >> s >> p.diameter;

    file >> s >> p.vtk_file_name;
    file >> s >> p.vtk_step;
    p.ny = p.u_in * p.n_y / p.RE;
    p.tau = 3 * (p.ny) + 0.5;
    file.close();
}
void swap_Lattice(Lattice &A, Lattice &B)
{
    int lenX = A.getX();
    int lenY = A.getY();
    for (int x = 0; x < lenX; x++)
    {
        for (int y = 0; y < lenY; y++)
        {
            for (int d = 0; d < 9; d++)
            {
                A.put(x,y,(direction) d,B.get(x,y,(direction) d));
            }
        }
    }
}
void printGrid(Lattice &l){
    cout << fixed;
    cout << setprecision(5);
    for(int y=l.getY()-1;y>=0;y--){
        for(int i=0;i<3;i++){
            for(int x =0;x<l.getX();x++){
                if(i==0){
                    cout <<"|"<< l(x,y,NW) << " "<< l(x,y,N) << " "<< l(x,y,NE) << "|";
                }
                if(i==1){
                    cout <<"|"<< l(x,y,W) << " "<<l(x,y,C) << " "<<l(x,y,E) << "|";
                }
                
                if(i==2){
                    cout <<"|"<< l(x,y,SW) << " "<<l(x,y,S) << " "<<l(x,y,SE) << "|";
                }
            }
            cout <<endl;
        }
        cout <<endl;
    }
}


int main(int argc, char *argv[])
{
    string fileName = argv[1];
    Parameters p;
    readInParameters(p, fileName);

    Lattice l(p.n_x + 2, p.n_y + 2);
    l.initialize_Latice(&p);
    Lattice l_cpy(p.n_x + 2, p.n_y + 2);
    l_cpy.initialize_Latice(&p);
    //Just for testig
    int index = 0;
    write_VTK_file(p.vtk_file_name, index, l, &p);
    index++;
    for (int t = 0, vtk = 1; t < p.timesteps; t++, vtk++)
    {
        collideStep(l,&p);
        
        borderUpdate(l, p.u_in, &p);
        streamPullStep(l, l_cpy);

        swap_Lattice(l, l_cpy);
        
        if (t % p.vtk_step == 0 || index == 0)
        {
            vtk = 1;
            write_VTK_file(p.vtk_file_name, index, l, &p);
            index++;
            cout << index << endl;
        }
    }
    return 0;
}
