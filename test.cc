#include <iostream>
#include "YAM2/yam2.h"

int main() {
    const yam2::FourMomentum v1{4.53329, -0.23492, 0.28285, 0.23056};
    const yam2::FourMomentum v2{3.47236, -1.19597, -0.065699, -1.41944};

    const yam2::TransverseMomentum ptmiss{1.43089, -0.21715};

    const yam2::Mass m_invis{0.0};

    const yam2::Mass m_parent{5.279};

    const double sqrt_s = 10.583;

    const double ptot_z = 0.0;

    const auto input = yam2::mkInput(v1, v2, ptmiss, m_invis, {m_parent}, sqrt_s, {ptot_z});

    if (!input) {
        std::cerr << "Input is invalid: " << std::endl;
        return 1;
    }
    std::cout << input.value() << std::endl;

    // const auto m2sol = yam2::m2CConsSQP(input, 1.0e-10, 1000);
    const auto m2sol = yam2::m2CCons(input, 1.0e-10, 1000);

    if (!m2sol) {
        std::cerr << "Failed to find minimum.\n";
        return 1;
    } else {
        // std::cout << m2sol.value() << '\n';

        /* k1 and k2 are the M2 solutions to the momenta of the invisible
         * particles X1 and X2.
         */
        std::cout << "M2CCons = " << m2sol.value().m2() << '\n'
                  << "where \n"
                  << "  k1: " << m2sol.value().k1() << '\n'
                  << "  k2: " << m2sol.value().k2() << '\n'
                  << "found after " << m2sol.value().neval_objf()
                  << " evaluations.\n";
    }
}

