#include <iostream>
#include "cnpy.h"
#include "YAM2/yam2.h"
#include <eigen3/Eigen/Dense>
#include <fstream>
#include "progressbar.hpp"

using namespace std;
using namespace Eigen;

static const int NROWS = 1080000;
static const int NCOLS = 4;
static const int NDBG = 10;

int main() {
    std::ofstream output{"data/m2ccb.dat"};
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
    int n = 0;

    // For Running
    progressbar bar(NROWS);
    for (int i = 0; i < NROWS; i++) {
        bar.update();
        const auto input = yam2::mkInput(
            {pa1_v[i], pa2_v[i]}, {pb1_v[i], pb2_v[i]}, pt_miss[i], m_c
        );
        const auto m2sol = yam2::m2CCB(input);

        if (m2sol) {
            const auto m2 = m2sol.value().m2();
            const auto k1 = m2sol.value().k1();
            const auto k2 = m2sol.value().k2();

            output << m2 << '\t' 
                << k1.e() << '\t'
                << k1.px() << '\t' 
                << k1.py() << '\t'
                << k1.pz() << '\t'
                << k2.e() << '\t'
                << k2.px() << '\t'
                << k2.py() << '\t'
                << k2.pz() << endl;
            n++;
        }
    }

    cout << "Run Finish" << endl;

    // =========================================================================
    // Write hdf5 file
    // =========================================================================

    cout << "Total data size: " << n << endl;

    output.close();

    cout << "Write Finish" << endl;

    return 0;
}
