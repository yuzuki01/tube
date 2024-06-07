#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <cmath>
#include <cstring>
#include <sstream>
#include <iomanip>
#include "omp.h"


namespace Tube {
    const int N_MESH = 50;
    const int N_DVS = 100;

    extern int ncell, nface;
    const double L0 = 1.0;
    const double R = 0.5;
    const int K = 0;
    extern double Ma, Kn, Pr;
    extern double gamma, Cv;
    extern double miu0, dt, half_dt;
    extern double solution_time;

    class Cell {
    public:
        int id;
        double x;
        double volume;
        std::vector<int> faces;
        std::vector<int> neighbors;

        Cell(int cell_id, double x, double volume) : id(cell_id), x(x), volume(volume),
                                                     neighbors({cell_id, cell_id + 1}) {};
    };

    class Face {
    public:
        int id;
        double x;
        std::array<int, 2> neighbors;

        Face(int face_id, double x) : id(face_id), x(x), neighbors({face_id, face_id}) {};
    };

    extern std::vector<Cell> cells;
    extern std::vector<Face> faces;

    typedef std::vector<double> Field;
    typedef std::vector<Field> DistributionFunction;

    extern DistributionFunction g_cell, h_cell;
    extern DistributionFunction g_face, h_face;
    extern Field rho_cell, rho_cell_n, rho_face;
    extern Field vel_cell, vel_cell_n, vel_face;
    extern Field T_cell, T_cell_n, T_face;
    extern Field tau_cell, tau_cell_n;
    extern Field q_cell, q_face;
    extern Field dvs_wgt, dvs_vel;


    void initial();
    inline double g_maxwell(double rho, double t, double cc);
    inline double h_maxwell(double t, double gm);
    inline double g_shakhov(double rho, double t, double cc, double cq, double gm);
    inline double h_shakhov(double rho, double t, double cc, double cq, double gm);

    void init_mesh(int N);
    void init_dvs(double scale);
    void output();
}
