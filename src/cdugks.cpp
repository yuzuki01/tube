#include "tube.h"


Tube::Tube(const String &config_file) : config(config_file) {
    int omp_threads = config.get("omp-threads", 4);
    std::cout << "Set omp threads num = " << omp_threads << std::endl;
    omp_set_num_threads(omp_threads);

    Kn = config.get("Kn", 1e-7);
    Pr = config.get("Pr", 1.0);
    R = config.get("gas-constant", 1.0);
    K = config.get("gas-k", 0);
    Rho0 = config.get("ref-density", 1.0);
    L0 = config.get("ref-length", 1.0);
    T0 = config.get("ref-temperature", 1.0);
    vhs_omega = config.get("vhs-omega", 0.5);
    vhs_index = config.get("vhs-index", 0.81);
    stop_time = config.get("stop-time", 0.14);
    dt = config.get("time-step", 0.0002);
    half_dt = 0.5 * dt;

    mesh.generate_mesh(config.get("ncell", 500));
    dvs.generate_dvs(config.get("dvs-mount", 33), config.get("dvs-scale", 8.0));

    std::cout << "Solver loaded." << std::endl;
}

void Tube::initial() {
    run_state = true;
    step = 0;
    solution_time = 0.0;
    gamma = (K + 5.0) / (K + 3.0);
    Cv = R / (gamma - 1.0);
    double RT = R * T0;
    miu0 = (Kn * L0) * (15 * Rho0 * sqrt(2.0 * M_PI * RT)) / ((7.0 - 2.0 * vhs_omega) * (5.0 - 2.0 * vhs_omega));

    rho_cell = mesh.zero_field(cell_field_flag);
    T_cell = mesh.zero_field(cell_field_flag);
    vel_cell = mesh.zero_field(cell_field_flag);
    q_cell = mesh.zero_field(cell_field_flag);

    g_cell.resize(dvs.NCELL, mesh.zero_field(cell_field_flag));
    h_cell.resize(dvs.NCELL, mesh.zero_field(cell_field_flag));
    g_face.resize(dvs.NCELL, mesh.zero_field(face_field_flag));
    h_face.resize(dvs.NCELL, mesh.zero_field(face_field_flag));
    /// flux
    flux_g.resize(dvs.NCELL, mesh.zero_field(cell_field_flag));
    flux_h.resize(dvs.NCELL, mesh.zero_field(cell_field_flag));

    for (auto &cell: mesh.cells) {
        double m0, m1, m2, m3;
        m0 = m1 = m2 = m3 = 0.0;
        {
            Scalar rho, u, T;
            if (cell.position <= 0.5) {
                rho = 1.0;
                u = 0.0;
                T = 1.0;
            } else {
                rho = 0.125;
                u = 0.0;
                T = 0.8;
            }
#pragma omp parallel for shared(rho, u, T, cell) reduction(+:m0) reduction(+:m1) reduction(+:m2) reduction(+:m3) default(none)
            for (int p = 0; p < dvs.NCELL; ++p) {
                auto &particle = dvs.cells[p];
                auto c = particle.position - u;
                auto cc = c * c;
                auto kk = particle.position * particle.position;
                auto g = g_maxwell(rho, T, cc);
                auto h = h_maxwell(T, g);
                g_cell[p][cell.id] = g;
                h_cell[p][cell.id] = h;
                m0 += particle.volume * g;
                m1 += particle.volume * g * particle.position;
                m2 += particle.volume * (kk * g + h);
                m3 += particle.volume * c * (cc * g + h);
            }
        }
        {
            auto rho = m0;
            auto rhoU = m1;
            auto rhoE = m2 * 0.5;
            auto q = m3 * 0.5;
            auto u = rhoU / rho;
            rho_cell[cell.id] = rho;
            vel_cell[cell.id] = u;
            T_cell[cell.id] = (rhoE / rho - 0.5 * u * u) / Cv;
            q_cell[cell.id] = q;
        }
    }
    std::cout << "Initialization finished." << std::endl;
}

inline Scalar Tube::tau_f(double rho, double t) const {
    return miu0 * pow(t / T0, vhs_index) / (rho * R * t);
}

inline Scalar Tube::g_maxwell(double rho, double t, double cc) const {
    double RT = R * t;
    return rho / sqrt(2.0 * M_PI * RT) * exp(-cc / (2.0 * RT));
}

inline Scalar Tube::h_maxwell(double t, double gm) const {
    return (K + 2.0) * R * t * gm;
}

inline Scalar Tube::g_shakhov(double rho, double t, double cc, double cq, double gm) const {
    /// g_shakhov = g_maxwell + g_pr
    double RT = R * t;
    double p = rho * RT;
    return (1.0 - Pr) * (cq / (5.0 * p * RT)) * (cc / RT - 3.0) * gm + gm;
}

inline Scalar Tube::h_shakhov(double rho, double t, double cc, double cq, double gm) const {
    /// h_shakhov = h_maxwell + h_pr
    double RT = R * t;
    double p = rho * RT;
    double h_m = h_maxwell(t, gm);
    return (1.0 - Pr) * (cq / (5.0 * p)) * ((cc / RT - 1.0) * (K + 2.0) - 2.0 * K) * gm + h_m;
}

void Tube::reconstruct() {
    /// get f_bar on face
#pragma omp parallel for default(none)
    for (int p = 0; p < dvs.NCELL; ++p) {
        auto &particle = dvs.cells[p];
        /// cell gradient - van leer
        Field grad_g = gradient(g_cell[p], mesh);
        Field grad_h = gradient(h_cell[p], mesh);
        /// interp to face
        for (auto &face: mesh.faces) {
            Scalar &nv = face.normal_vector[0];
            auto &cell = (nv * particle.position >= 0.0) ? mesh.cells[face.cell_id[0]] : mesh.cells[face.cell_id[1]];
            Scalar dr_ij = face.position - cell.position;
            g_face[p][face.id] = g_cell[p][cell.id]
                                 + (dr_ij - particle.position * half_dt) * (grad_g[cell.id]);
            h_face[p][face.id] = h_cell[p][cell.id]
                                 + (dr_ij - particle.position * half_dt) * (grad_h[cell.id]);
        }
    }
    {
        auto m0_face = mesh.zero_field(face_field_flag);
        auto m1_face = mesh.zero_field(face_field_flag);
        auto m2_face = mesh.zero_field(face_field_flag);
        /// face macro vars
        for (auto &face: mesh.faces) {
            double m0, m1, m2;
            m0 = m1 = m2 = 0.0;
#pragma omp parallel for shared(face) reduction(+:m0) reduction(+:m1) reduction(+:m2) default(none)
            for (int p = 0; p < dvs.NCELL; ++p) {
                auto &particle = dvs.cells[p];
                auto kk = particle.position * particle.position;
                auto g = g_face[p][face.id];
                auto h = h_face[p][face.id];
                m0 += particle.volume * g;
                m1 += particle.volume * g * particle.position;
                m2 += particle.volume * (kk * g + h);
            }
            m0_face[face.id] = m0;
            m1_face[face.id] = m1;
            m2_face[face.id] = m2;
        }
        auto m3_face = mesh.zero_field(face_field_flag);
        for (auto &face: mesh.faces) {
            double m3;
            m3 = 0.0;
            auto rho = m0_face[face.id];
            auto rhoU = m1_face[face.id];
            auto u = rhoU / rho;
#pragma omp parallel for shared(face, u) reduction(+:m3) default(none)
            for (int p = 0; p < dvs.NCELL; ++p) {
                auto &particle = dvs.cells[p];
                auto c = particle.position - u;
                auto cc = c * c;
                auto g = g_face[p][face.id];
                auto h = h_face[p][face.id];
                m3 += particle.volume * c * (cc * g + h);
            }
            m3_face[face.id] = m3;
        }
        /// get original f on face
        for (auto &face: mesh.faces) {
            auto rho = m0_face[face.id];
            auto rhoU = m1_face[face.id];
            auto rhoE = 0.5 * m2_face[face.id];
            auto u = rhoU / rho;
            auto T = (rhoE / rho - 0.5 * u * u) / Cv;
            auto tau = tau_f(rho, T);
            auto q = (tau / (2.0 * tau + half_dt * Pr)) * m3_face[face.id];
            auto C_s = half_dt / (2.0 * tau + half_dt);
#pragma omp parallel for shared(face, rho, u, T, q, C_s) default(none)
            for (int p = 0; p < dvs.NCELL; ++p) {
                auto &particle = dvs.cells[p];
                auto c = particle.position - u;
                auto cc = c * c;
                auto cq = c * q;
                auto g_m = g_maxwell(rho, T, cc);
                auto g_s = g_shakhov(rho, T, cc, cq, g_m);
                auto h_s = h_shakhov(rho, T, cc, cq, g_m);
                g_face[p][face.id] = (1.0 - C_s) * g_face[p][face.id] + C_s * g_s;
                h_face[p][face.id] = (1.0 - C_s) * h_face[p][face.id] + C_s * h_s;
            }
        }
    }
    /// boundary
    {
        auto &face = mesh.faces[0];
        auto &neighbor = mesh.cells[face.cell_id[0]];
#pragma omp parallel for shared(face, neighbor) default(none)
        for (int p = 0; p < dvs.NCELL; ++p) {
            auto &particle = dvs.cells[p];
            auto kn = particle.position * face.normal_vector[1];
            if (kn >= 0.0) {
                auto c = particle.position - vel_cell[neighbor.id];
                auto cc = c * c;
                auto rho = rho_cell[neighbor.id];
                auto T = T_cell[neighbor.id];
                auto g_m = g_maxwell(rho, T, cc);
                g_face[p][face.id] = g_m;
                h_face[p][face.id] = h_maxwell(T, g_m);
            }
        }
    }
    {
        auto &face = mesh.faces[mesh.NFACE - 1];
        auto &neighbor = mesh.cells[face.cell_id[0]];
#pragma omp parallel for shared(face, neighbor) default(none)
        for (int p = 0; p < dvs.NCELL; ++p) {
            auto &particle = dvs.cells[p];
            auto kn = particle.position * face.normal_vector[1];
            if (kn >= 0.0) {
                auto c = particle.position - vel_cell[neighbor.id];
                auto cc = c * c;
                auto rho = rho_cell[neighbor.id];
                auto T = T_cell[neighbor.id];
                auto g_m = g_maxwell(rho, T, cc);
                g_face[p][face.id] = g_m;
                h_face[p][face.id] = h_maxwell(T, g_m);
            }
        }
    }
}

void Tube::fvm_update() {
    auto flux_m0 = mesh.zero_field();
    auto flux_m1 = mesh.zero_field();
    auto flux_m2 = mesh.zero_field();
    for (auto &cell: mesh.cells) {
        double flux_m0_c, flux_m1_c, flux_m2_c;
        flux_m0_c = flux_m1_c = flux_m2_c = 0.0;
#pragma omp parallel for shared(cell) reduction(+:flux_m0_c) reduction(+:flux_m1_c) reduction(+:flux_m2_c) default(none)
        for (int p = 0; p < dvs.NCELL; ++p) {
            auto &particle = dvs.cells[p];
            double flux_g_tmp = 0.0, flux_h_tmp = 0.0;
            for (auto face_id: cell.face_id) {
                auto &face = mesh.faces[face_id];
                auto &nv = (face.cell_id[0] == cell.id) ? face.normal_vector[0] : face.normal_vector[1];
                auto knA = (particle.position * nv);
                flux_g_tmp += knA * g_face[p][face.id];
                flux_h_tmp += knA * h_face[p][face.id];
            }
            flux_g[p][cell.id] = flux_g_tmp;
            flux_h[p][cell.id] = flux_h_tmp;
            flux_m0_c += particle.volume * flux_g_tmp;
            flux_m1_c += particle.volume * flux_g_tmp * particle.position;
            flux_m2_c += particle.volume * ((particle.position * particle.position) * flux_g_tmp + flux_h_tmp);
        }
        flux_m0[cell.id] = flux_m0_c;
        flux_m1[cell.id] = flux_m1_c;
        flux_m2[cell.id] = flux_m2_c;
    }

    auto m3_cell = mesh.zero_field();
    for (auto &cell: mesh.cells) {
        auto dt_v = dt / cell.volume;
        /// tn = n
        auto rho_n = rho_cell[cell.id];
        auto u_n = vel_cell[cell.id];
        auto T_n = T_cell[cell.id];
        auto q_n = q_cell[cell.id];
        auto rhoU_n = rho_n * u_n;
        auto rhoE_n = rho_n * (0.5 * (u_n * u_n) + Cv * T_n);
        auto tau_n = tau_f(rho_n, T_n);
        /// tn = n + 1
        auto rho = rho_n - dt_v * flux_m0[cell.id];
        auto rhoU = rhoU_n - dt_v * flux_m1[cell.id];
        auto rhoE = rhoE_n - dt_v * 0.5 * flux_m2[cell.id];
        auto u = rhoU / rho;
        auto T = (rhoE / rho - (u * u) * 0.5) / Cv;
        auto tau = tau_f(rho, T);
        /// Cell Macro Vars
        rho_cell[cell.id] = rho;
        vel_cell[cell.id] = u;
        T_cell[cell.id] = T;
        /// Evolution
        auto cm = half_dt / (half_dt - 2.0 * tau_n);
        auto C = tau / (tau + half_dt);
        auto C_s = half_dt / (2.0 * tau);
        double m3;
        m3 = 0.0;
#pragma omp parallel for shared(cell, dt_v, cm, C, C_s, rho_n, u_n, T_n, q_n, tau_n, rho, u, T, tau) reduction(+:m3) default(none)
        for (int p = 0; p < dvs.NCELL; ++p) {
            auto &particle = dvs.cells[p];
            /// tn = n
            auto c_n = particle.position - u_n;
            auto cc_n = c_n * c_n;
            auto cq_n = c_n * q_n;
            auto g_m_n = g_maxwell(rho_n, T_n, cc_n);
            auto g_s_n = g_shakhov(rho_n, T_n, cc_n, cq_n, g_m_n);
            auto h_s_n = h_shakhov(rho_n, T_n, cc_n, cq_n, g_m_n);
            auto g_n = cm * g_s_n + (1.0 - cm) * g_cell[p][cell.id];
            auto h_n = cm * h_s_n + (1.0 - cm) * h_cell[p][cell.id];
            /// tn = n + 1
            auto c = particle.position - u;
            auto cc = c * c;
            auto cq = c * q_n;
            auto g_m = g_maxwell(rho, T, cc);
            auto g_s = g_shakhov(rho, T, cc, cq, g_m);
            auto h_s = h_shakhov(rho, T, cc, cq, g_m);
            /// Evolution
            auto g = C * (g_n + half_dt * (g_s / tau + (g_s_n - g_n) / tau_n) - dt_v * flux_g[p][cell.id]);
            auto h = C * (h_n + half_dt * (h_s / tau + (h_s_n - h_n) / tau_n) - dt_v * flux_h[p][cell.id]);
            /// heat-flux
            m3 += particle.volume * c * (cc * g + h);
            /// f -> f_bar_plus
            g_cell[p][cell.id] = (1.0 - C_s) * g + C_s * g_s;
            h_cell[p][cell.id] = (1.0 - C_s) * h + C_s * h_s;
        }
        m3_cell[cell.id] = 0.5 * m3;
    }
    /// update heat-flux
    q_cell = m3_cell;
}

void Tube::output() const {
    system("if not exist result mkdir result");
    std::stringstream ss;
    ss << "./result/step-" << step << ".dat.plt";
    std::ofstream file;
    file.open(ss.str());
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << ss.str() << std::endl;
        throw std::invalid_argument("Failed to open file!");
    }
    /// Write
    file << R"(VARIABLES= "X", "Rho", "U", "T", "q")" << std::endl;
    file << "ZONE T=\"Step-" << step << "\", I=" << mesh.NCELL << ", DATAPACKING=POINT, SOLUTIONTIME=" << solution_time
         << std::endl << std::endl;
    for (auto &cell: mesh.cells) {
        file << " " << cell.position << " " << rho_cell[cell.id]
             << " " << vel_cell[cell.id] << " " << T_cell[cell.id]
             << " " << q_cell[cell.id] << std::endl;
    }
}

void Tube::do_step() {
    reconstruct();
    fvm_update();
    step++;
    solution_time += dt;
    if (solution_time + dt > stop_time and solution_time < stop_time) {
        dt = stop_time - solution_time;
        half_dt = 0.5 * dt;
        reconstruct();
        fvm_update();
        step++;
        solution_time = stop_time;
        run_state = false;
    }
}
