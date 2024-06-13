#include "tube.h"


Cell::Cell(ObjectId id) : id(id), face_id({id, id + 1}) {}

Face::Face(ObjectId id) : id(id) {}

Field nc_weight = {1. / 8., 3. / 8., 3. / 8., 1. / 8.};

void Mesh::generate_dvs(int mount, double scale) {
    const int n = 3;
    NCELL = n * mount + 1;
    double dh = 2.0 * scale / mount;
    Field coordinate(NCELL), weight(NCELL);
    for (int i = 0; i < NCELL; i++) {
        coordinate[i] = scale * double(2 * i + 1 - NCELL) / double(NCELL - 1);
        weight[i] = 0.0;
    }
    for (int i = 0; i < mount; ++i) {
        for (int j = 0; j < n + 1; ++j) {
            weight[i * n + j] += nc_weight[j] * dh;
        }
    }
    max_cell_mag = 0.0;
    for (int i = 0; i < NCELL; ++i) {
        Cell cell(i);
        cell.position = coordinate[i];
        cell.volume = weight[i];
        cells.push_back(cell);
        double cell_mag = cell.position;
        if (max_cell_mag < cell_mag) max_cell_mag = cell_mag;
    }
}

void Mesh::generate_mesh(int ncell) {
    NCELL = ncell;
    NFACE = ncell + 1;
    cells.reserve(NCELL);
    faces.reserve(NFACE);
    /// generate face
    for (int i = 0; i < NFACE; ++i) {
        Face face(i);
        face.position = double(i) / double(NFACE - 1);
        if (i == 0) {
            // left
            face.cell_id = {i, i};
            face.normal_vector = {-1, 1};
        } else if (i == NFACE - 1) {
            // right
            face.cell_id = {i - 1, i - 1};
            face.normal_vector = {1, -1};
        } else {
            // interior
            face.cell_id = {i - 1, i};
            face.normal_vector = {1, -1};
        }
        faces.push_back(face);
    }
    min_cell_size = faces[1].position - faces[0].position;
    /// generate cell
    for (int i = 0; i < NCELL; ++i) {
        Cell cell(i);
        auto &f0 = faces[cell.face_id[0]];
        auto &f1 = faces[cell.face_id[1]];
        cell.position = (f0.position + f1.position) * 0.5;
        cell.volume = f1.position - f0.position;
        auto cell_size = 0.5 * cell.volume;
        if (min_cell_size > cell_size) min_cell_size = cell_size;
        if (i == 0) {
            // left
            cell.neighbors = {i + 1};
        } else if (i == NCELL - 1) {
            // right
            cell.neighbors = {i - 1};
        } else {
            // interior
            cell.neighbors = {i - 1, i + 1};
        }
        cells.push_back(cell);
    }
}

Field Mesh::zero_field(int flag) const {
    size_t size;
    if (flag == cell_field_flag) {
        size = NCELL;
    } else {
        size = NFACE;
    }
    return Field(size, 0.0);
}

inline int sgn(Scalar _x) {
    if (_x == 0.0) return 0;
    if (_x < 0.0) return -1;
    return 1;
}

Field gradient(Field &field, Mesh &mesh) {
    Field result(mesh.NCELL);
    for (int i = 0; i < mesh.NCELL; ++i) {
        if (i == 0 or i == mesh.NCELL - 1) {
            result[i] = 0.0;
        } else {
            auto &c0 = mesh.cells[i - 1];
            auto &c1 = mesh.cells[i];
            auto &c2 = mesh.cells[i + 1];
            auto &f0 = field[i - 1];
            auto &f1 = field[i];
            auto &f2 = field[i + 1];
            auto s1 = (f1 - f0) / (c1.position - c0.position);
            auto s2 = (f2 - f1) / (c2.position - c1.position);
            auto s1_abs = std::fabs(s1);
            auto s2_abs = std::fabs(s2);
            if ((s1_abs + s2_abs) == 0.0) {
                result[i] = 0.0;
            } else {
                result[i] = (sgn(s1) + sgn(s2)) * (s1_abs * s2_abs) / (s1_abs + s2_abs);
            }
        }
    }
    return result;
}
