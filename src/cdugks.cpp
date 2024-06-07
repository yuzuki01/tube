#include "tube.h"


using namespace Tube;

double Tube::gamma, Tube::Cv;

DistributionFunction Tube::g_cell, Tube::h_cell;
DistributionFunction Tube::g_face, Tube::h_face;
Field Tube::rho_cell, Tube::rho_cell_n, Tube::rho_face;
Field Tube::vel_cell, Tube::vel_cell_n, Tube::vel_face;
Field Tube::T_cell, Tube::T_cell_n, Tube::T_face;
Field Tube::tau_cell, Tube::tau_cell_n;
Field Tube::q_cell, Tube::q_face;

void Tube::initial() {
    init_mesh(N_MESH);
    init_dvs(5.0);

    g_cell.resize(N_DVS, Field(ncell));
    h_cell.resize(N_DVS, Field(ncell));
    g_face.resize(N_DVS, Field(nface));
    h_face.resize(N_DVS, Field(nface));

    rho_cell = Field(ncell);
    rho_cell_n = Field(ncell);
    rho_face = Field(nface);

    vel_cell = Field(ncell);
    vel_cell_n = Field(ncell);
    vel_face = Field(nface);

    T_cell = Field(ncell);
    T_cell_n = Field(ncell);
    T_face = Field(ncell);

    tau_cell = Field(ncell);
    tau_cell_n = Field(ncell);

    q_cell = Field(ncell);

    rho_face = Field(nface);
    vel_face = Field(nface);
    T_face = Field(nface);

    q_face = Field(nface);

    gamma = (K + 5.0) / (K + 3.0);
    Cv = R / (gamma - 1.0);

    for (auto &cell : cells) {
        double u = 0.0;
        double rho = (2.0 * cell.id < N_MESH) ? 1.0 : 0.125;
        double T = (2.0 * cell.id < N_MESH) ? 1.0 : 0.8;
        double m0 = 0.0, m1 = 0.0, m2 = 0.0, m3 = 0.0;
        for (int p = 0; p < N_DVS; ++p) {
            auto w = dvs_wgt[p];
            auto k = dvs_vel[p];
            auto c = k - u;
            auto cc = c * c;
            auto g = g_maxwell(rho, T, cc);
            auto h = h_maxwell(T, g);

            g_cell[p][cell.id] = g;
            h_cell[p][cell.id] = h;

            m0 += w * g;
            m1 += w * k * g;
            m2 += w * (k * k * g + h);
            m3 += w * c * (cc * g + h);
        }
        u = m1 / m0;
        rho_cell[cell.id] = m0;
        vel_cell[cell.id] = u;
        T_cell[cell.id] = (m2 / m0 - u * u) / (2.0 * Cv);
        q_cell[cell.id] = 0.5 * m3;
    }
    output();
}


inline double Tube::g_maxwell(double rho, double t, double cc) {
    double RT = R * t;
    return rho / sqrt(2.0 * M_PI * RT) * exp(-cc / (2.0 * RT));
}

inline double Tube::h_maxwell(double t, double gm) {
    return (K + 2.0) * R * t * gm;
}

inline double Tube::g_shakhov(double rho, double t, double cc, double cq, double gm) {
    /// g_shakhov = g_maxwell + g_pr
    double RT = R * t;
    double p = rho * RT;
    return (1.0 - Pr) * (cq / (5.0 * p * RT)) * (cc / RT - 3.0) * gm + gm;
}

inline double Tube::h_shakhov(double rho, double t, double cc, double cq, double gm) {
    /// h_shakhov = h_maxwell + h_pr
    double RT = R * t;
    double p = rho * RT;
    double h_m = h_maxwell(t, gm);
    return (1.0 - Pr) * (cq / (5.0 * p)) * ((cc / RT - 1.0) * (K + 2.0) - 2.0 * K) * gm + h_m;
}
