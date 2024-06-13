#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "omp.h"


/// type def
typedef int ObjectId;
typedef std::vector<ObjectId> ObjectIdList;
typedef std::array<ObjectId, 2> ObjectIdSet;
typedef double Scalar;
typedef std::array<Scalar, 2> VectorSet;
typedef std::vector<Scalar> Field;
typedef std::vector<Field> DistributionFunction;

/// tube head
#include "mesh.h"
#include "config.h"

class Tube {
private:
    bool run_state = false;
public:
    Config config;
    Mesh mesh, dvs;

    Scalar Kn, Pr;
    Scalar R, Rho0, T0, L0;
    Scalar vhs_omega, vhs_index;
    Scalar gamma{}, Cv{};
    Scalar miu0{}, dt, half_dt, solution_time{};
    Scalar stop_time;
    int K;
    int step{};

    Field rho_cell, T_cell;
    Field vel_cell, q_cell;

    DistributionFunction g_cell, h_cell;
    DistributionFunction g_face, h_face;
    DistributionFunction flux_g, flux_h;

    explicit Tube(const String& config_file);

    void initial();

    inline Scalar tau_f(Scalar rho, Scalar t) const;

    inline Scalar g_maxwell(Scalar rho, Scalar t, Scalar cc) const;

    inline Scalar g_shakhov(Scalar rho, Scalar t, Scalar cc, Scalar cq, Scalar gm) const;

    inline Scalar h_maxwell(Scalar t, Scalar gm) const;

    inline Scalar h_shakhov(Scalar rho, Scalar t, Scalar cc, Scalar cq, Scalar gm) const;

    void reconstruct();

    void fvm_update();

    void do_step();

    void output() const;

    bool run_status() const { return run_state; };
};

Field gradient(Field &field, Mesh &mesh);

void output(const String &file_name, Field &field);
