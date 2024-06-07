#include "tube.h"


using namespace Tube;

int Tube::ncell, Tube::nface;

void Tube::init_mesh(int N) {
    ncell = N;
    nface = N + 1;
    for (int i = 0; i < N + 1; ++i) {
        auto pos = L0 * (double(i) / double(N));
        faces.emplace_back(i, pos);
    }
    for (int i = 0; i < N; ++i) {
        auto pos = (faces[i].x + faces[i + 1].x) / 2.0;
        auto vol = faces[i + 1].x - faces[i].x;
        cells.emplace_back(i, pos, vol);
    }
    for (int i = 1; i < N + 1; ++i) {
        faces[i].neighbors[0] = i - 1;
    }
}


Field Tube::dvs_wgt(N_DVS), Tube::dvs_vel(N_DVS);

void Tube::init_dvs(double scale) {
    for (int i = 0; i < N_DVS; ++i) {
        dvs_vel[i] = scale * ((2.0 * i + 1.0 - N_DVS) / (N_DVS - 1.0));
        dvs_wgt[i] = 2.0 * scale / (N_DVS - 1.0);
    }
}

void Tube::output() {
    std::ofstream fp;
}
