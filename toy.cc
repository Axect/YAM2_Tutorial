#include <iostream>
#include "cnpy.h"
#include "YAM2/yam2.h"
#include <eigen3/Eigen/Dense>
#include <h5pp/h5pp.h>
#include "progressbar.hpp"

using namespace std;
using namespace Eigen;

static const int NROWS = 1080000;
static const int NCOLS = 4;
static const int NDBG = 10;

int main() {
    // =========================================================================
    // Read npz file
    // =========================================================================
    auto data = cnpy::npz_load("data/toy_array_new.npz");
    
    auto b1 = data["b1"];
    double* pa1 = b1.data<double>();
 
    auto b2 = data["b2"];
    double* pa2 = b2.data<double>();

    auto l2 = data["l2"];
    double* pb1 = l2.data<double>();

    auto l1 = data["l1"];
    double* pb2 = l1.data<double>();

    // =========================================================================
    // netCDF to Eigen
    // =========================================================================
    auto pa1_e = Map<MatrixXd>(pa1, b1.shape[1], b1.shape[0]).adjoint();
    auto pa2_e = Map<MatrixXd>(pa2, b2.shape[1], b2.shape[0]).adjoint();
    auto pb1_e = Map<MatrixXd>(pb1, l2.shape[1], l2.shape[0]).adjoint();
    auto pb2_e = Map<MatrixXd>(pb2, l1.shape[1], l1.shape[0]).adjoint();

    cout << "pa1_e.rows(): " << pa1_e.rows() << " pa1_e.cols(): " << pa1_e.cols() << endl;

    auto pt_x = -(pa1_e.col(1) + pa2_e.col(1) + pb1_e.col(1) + pb2_e.col(1));
    auto pt_y = -(pa1_e.col(2) + pa2_e.col(2) + pb1_e.col(2) + pb2_e.col(2));

    // =========================================================================
    // Eigen to FourMomentum
    // =========================================================================
    vector<yam2::FourMomentum> pa1_v(NROWS);
    vector<yam2::FourMomentum> pa2_v(NROWS);
    vector<yam2::FourMomentum> pb1_v(NROWS);
    vector<yam2::FourMomentum> pb2_v(NROWS);
    const yam2::TransverseMomentum ptmiss{0.0, 0.0};
    vector<yam2::TransverseMomentum> pt_miss(NROWS, ptmiss);

    for (int i = 0; i < NROWS; i++) {
        pa1_v[i] = yam2::FourMomentum(pa1_e(i, 0), pa1_e(i, 1), pa1_e(i, 2), pa1_e(i, 3));
        pa2_v[i] = yam2::FourMomentum(pa2_e(i, 0), pa2_e(i, 1), pa2_e(i, 2), pa2_e(i, 3));
        pb1_v[i] = yam2::FourMomentum(pb1_e(i, 0), pb1_e(i, 1), pb1_e(i, 2), pb1_e(i, 3));
        pb2_v[i] = yam2::FourMomentum(pb2_e(i, 0), pb2_e(i, 1), pb2_e(i, 2), pb2_e(i, 3));
        pt_miss[i] = yam2::TransverseMomentum(pt_x(i), pt_y(i));
    }

    cout << pt_miss[0] << endl;
    cout << pa1_v[0].m() << endl;
    cout << pa2_v[0].m() << endl;
    cout << pb1_v[0].m() << endl;
    cout << pb2_v[0].m() << endl;

    // =========================================================================
    // M2CC
    // =========================================================================
    const yam2::Mass m_c{700.0};

    vector<double> m2cc_dbg(NDBG, 0.0);
    vector<yam2::FourMomentum> k1_dbg(NDBG);
    vector<yam2::FourMomentum> k2_dbg(NDBG);

    // For Debug
    for (int i = 0; i < NDBG; i++) {
        const auto input = yam2::mkInput(
            {pa1_v[i], pa2_v[i]}, {pb1_v[i], pb2_v[i]}, pt_miss[i], m_c
        );
        const auto m2sol = yam2::m2CC(input);

        if (!m2sol) {
            cerr << "m2CC failed for " << i << endl;
            cout << pa1_v[i].m() << endl;
            cout << pa2_v[i].m() << endl;
            cout << pb1_v[i].m() << endl;
            cout << pb2_v[i].m() << endl;
            // return 1;
        } else {
            m2cc_dbg[i] = m2sol.value().m2();
            k1_dbg[i] = m2sol.value().k1();
            k2_dbg[i] = m2sol.value().k2();
        }
    }

    vector<double> m2cc;
    vector<yam2::FourMomentum> k1;
    vector<yam2::FourMomentum> k2;
    vector<vector<double>> pa1_new;
    vector<vector<double>> pa2_new;
    vector<vector<double>> pb1_new;
    vector<vector<double>> pb2_new;

    // For Running
    progressbar bar(NROWS);
    for (int i = 0; i < NROWS; i++) {
        bar.update();
        const auto input = yam2::mkInput(
            {pa1_v[i], pa2_v[i]}, {pb1_v[i], pb2_v[i]}, pt_miss[i], m_c
        );
        const auto m2sol = yam2::m2CC(input);

        if (m2sol) {
            m2cc.push_back(m2sol.value().m2());
            k1.push_back(m2sol.value().k1());
            k2.push_back(m2sol.value().k2());
            pa1_new.push_back({pa1_e(i, 0), pa1_e(i, 1), pa1_e(i, 2), pa1_e(i, 3)});
            pa2_new.push_back({pa2_e(i, 0), pa2_e(i, 1), pa2_e(i, 2), pa2_e(i, 3)});
            pb1_new.push_back({pb1_e(i, 0), pb1_e(i, 1), pb1_e(i, 2), pb1_e(i, 3)});
            pb2_new.push_back({pb2_e(i, 0), pb2_e(i, 1), pb2_e(i, 2), pb2_e(i, 3)});
        }
    }

    cout << "m2cc[0] = " << m2cc[0] << endl;
    cout << "k1[0] = " << k1[0] << endl;
    cout << "k2[0] = " << k2[0] << endl;

    cout << "Run Finish" << endl;

    // =========================================================================
    // Write hdf5 file
    // =========================================================================
    auto row = m2cc.size();
    VectorXd m2cc_e(row);
    MatrixX4d k1_e(row, 4);
    MatrixX4d k2_e(row, 4);
    MatrixX4d pa1_new_e(row, 4);
    MatrixX4d pa2_new_e(row, 4);
    MatrixX4d pb1_new_e(row, 4);
    MatrixX4d pb2_new_e(row, 4);

    for (int i = 0; i < row; i++) {
        m2cc_e(i) = m2cc[i];

        k1_e(i, 0) = k1[i].e();
        k1_e(i, 1) = k1[i].px();
        k1_e(i, 2) = k1[i].py();
        k1_e(i, 3) = k1[i].pz();

        k2_e(i, 0) = k2[i].e();
        k2_e(i, 1) = k2[i].px();
        k2_e(i, 2) = k2[i].py();
        k2_e(i, 3) = k2[i].pz();

        pa1_new_e(i, 0) = pa1_new[i][0];
        pa1_new_e(i, 1) = pa1_new[i][1];
        pa1_new_e(i, 2) = pa1_new[i][2];
        pa1_new_e(i, 3) = pa1_new[i][3];

        pa2_new_e(i, 0) = pa2_new[i][0];
        pa2_new_e(i, 1) = pa2_new[i][1];
        pa2_new_e(i, 2) = pa2_new[i][2];
        pa2_new_e(i, 3) = pa2_new[i][3];

        pb1_new_e(i, 0) = pb1_new[i][0];
        pb1_new_e(i, 1) = pb1_new[i][1];
        pb1_new_e(i, 2) = pb1_new[i][2];
        pb1_new_e(i, 3) = pb1_new[i][3];

        pb2_new_e(i, 0) = pb2_new[i][0];
        pb2_new_e(i, 1) = pb2_new[i][1];
        pb2_new_e(i, 2) = pb2_new[i][2];
        pb2_new_e(i, 3) = pb2_new[i][3];
    }

    cout << "Total data size: " << row << endl;

    // Write to hdf5 file
    h5pp::File file("data/m2cc.h5");
    file.writeDataset(m2cc_e, "m2cc");
    file.writeDataset(k1_e, "k1");
    file.writeDataset(k2_e, "k2");
    file.writeDataset(pa1_new_e, "pa1");
    file.writeDataset(pa2_new_e, "pa2");
    file.writeDataset(pb1_new_e, "pb1");
    file.writeDataset(pb2_new_e, "pb2");

    cout << "Write Finish" << endl;

    return 0;
}
