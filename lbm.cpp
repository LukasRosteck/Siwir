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



//TODO: Implement the Latice class and replace the Grid clas with it
typedef struct Grid
{
    int num_points;
    int num_cells;

    //Stores the point coordinates x,y,z and the value v. (z coordinate is unused in this programm but is stored for completion.)
    vector<vector<double>> points;
    vector<vector<int>> cells;
} Grid;

//Cell Lookup map in order to more transparently access the positions on the cells.
map<string, int> cell_lookup = {{"C",0},{"N",1},{"S",2},{"W",3},{"E",4},{"NW",5},{"NE",6},{"SW",7},{"SE",8}};


template<typename T> class Lattice
{
    int n_cells_y, n_cells_x;

    vector<vector<vector<T>>> domain;

    vector<vector<int>> flag_field;

public:
    Lattice(int N_y, int N_x): n_cells_y(N_y), n_cells_x(N_x)
    {
        domain = vector<vector<vector<T>>>(n_cells_x, vector<vector<T>>(n_cells_y, vector<T>(9)));
        flag_field = vector<vector<int>>(n_cells_x, vector<int>(n_cells_y, 0));
    }

    //TODO: Introduce parameters for the circle and write method to initialize the Lattice Grid
    void initialize_Latice();

    void put(int, int, string, T);
    T get(int, int, string);
    int getX();
    int getY();
    T &operator()(const int, const int, string);
};


template<typename T> void Lattice<T>::initialize_Latice()
{
    //TODO: Introduce parameters for the circle and write method to initialize IMPLEMENTATION
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

template<typename T> T &Lattice<T>::operator()(const int x, const int y, string q)
{
    return domain[x][y][cell_lookup[q]];
}


//TODO: Use the Latice class to adapt this vtk file writing mehthod. This method should give a basic overview of the writing structure.
//----------------------------------------------------------------
//Write the given data to a legacy VTK file with the specified filename.
//In this method the given outfile is completely overwritten.
//----------------------------------------------------------------
void write_VTK_file(string &filename, Grid *outgrid)
{
    ofstream outstream;
    outstream.open(filename);

    //Write the header
    outstream << "# vtk DataFile Version 2.0" << endl
              << "Output triangulation" << endl
              << "ASCII" << endl
              << "DATASET UNSTRUCTURED_GRID" << endl;
    outstream << endl;
    outstream << "POINTS " << outgrid->num_points << " float" << endl;
    //Write all points
    for (int i = 0; i < outgrid->num_points; i++)
    {
        outstream << outgrid->points[i][0] << " " << outgrid->points[i][1] << " " << outgrid->points[i][2] << endl;
    }
    outstream << endl;

    //Write all cells
    outstream << "CELLS " << outgrid->num_cells << " " << outgrid->num_cells * 4 << endl;
    for (int i = 0; i < outgrid->num_cells; i++)
    {
        outstream << 3 << " " << outgrid->cells[i][0] << " " << outgrid->cells[i][1] << " " << outgrid->cells[i][2] << endl;
    }
    outstream << endl;

    //Write cell type which is fixed as the triangle (5) in this programm
    outstream << "CELL_TYPES " << outgrid->num_cells << endl;
    for (int i = 0; i < outgrid->num_cells; i++)
    {
        outstream << 5 << endl;
    }
    outstream << endl;

    outstream << "POINT_DATA " << outgrid->num_points << endl
              << "SCALARS value float" << endl
              << "LOOKUP_TABLE default" << endl;

    //Write the values at each point
    for (int i = 0; i < outgrid->num_points; i++)
    {
        outstream << outgrid->points[i][3] << endl;
    }

    return;
}



//TODO: Read in the parameters and call a to be defined latice method which takes care of the simulation.
int main(int argc, char* argv[])
{
    cout << "Very quiet here..." << endl;

    return 0;
}
