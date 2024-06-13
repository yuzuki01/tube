#ifndef TUBE_MESH_H
#define TUBE_MESH_H


class Cell {
public:
    ObjectId id;
    ObjectId group_id{};
    Scalar position{};
    Scalar volume{};
    ObjectIdSet face_id;
    ObjectIdList neighbors;

    explicit Cell(ObjectId id);
};


class Face {
public:
    ObjectId id;
    Scalar position{};
    ObjectIdSet cell_id{};
    VectorSet normal_vector{};

    explicit Face(ObjectId id);
};


typedef std::vector<Cell> CellList;
typedef std::vector<Face> FaceList;


enum {cell_field_flag, face_field_flag};

class Mesh {
public:
    int NCELL{}, NFACE{};
    CellList cells;
    FaceList faces;

    double min_cell_size{};
    double max_cell_mag{};

    Mesh() = default;
    void generate_mesh(int ncell);
    void generate_dvs(int mount, double scale);
    [[nodiscard]] Field zero_field(int flag=cell_field_flag) const;
};

#endif //TUBE_MESH_H
