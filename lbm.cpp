#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
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

template<typename T>
class Lattice
{
    int n_cells_y, n_cells_x;

    vector<vector<vector<T>>> domain;

    vector<vector<int>> cell_types;

public:
    Lattice(int N_y, int N_x): n_cells_y(N_y), n_cells_x(N_x)
    {
        domain.reserve(n_cells_x);
        cell_types.reserve(n_cells_x );

        for (int i = 0; i < n_cells_x; i++)
        {
            vector<vector<T>> tmp1;
            tmp1.reserve(n_cells_y);
            domain[i].push_back(tmp1);

            vector<int> tmp2;
            tmp2.reserve(n_cells_y);
            cell_types.push_back(tmp2);

            for(int j = 0; j < n_cells_y; j++)
            {
                 vector<T> tmp3(9);
                 domain[i][j].push_back(tmp3);
            }
        }
    }

    //TODO: Introduce parameters for the circle and write method to initialize the Lattice Grid
    void initialize_Latice();

    void put(int, int, T);
    T get(int, int);
    int getX();
    int getY();
    T &operator()(const int, const int);
};



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
