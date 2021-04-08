//
// Created by john on 1/5/21.
//

#ifndef CHITECH_CUSTOM_QUADRATURE_H
#define CHITECH_CUSTOM_QUADRATURE_H

#include "angular_quadrature_base.h"

namespace chi_math
{
    class CustomQuadrature : public AngularQuadrature
    {
    private:
        const std::string quad_file;
        const std::string M2D_file;
        const std::string D2M_file;

    public:
        explicit
        CustomQuadrature(std::string& in_filename1, std::string& in_filename2, std::string& in_filename3) :
        AngularQuadrature(AngularQuadratureType::Arbitrary),
        quad_file(in_filename1),
        M2D_file(in_filename2),
        D2M_file(in_filename3)
        {}

        void BuildDiscreteToMomentOperator(int scatt_order, bool oneD) override;
        void BuildMomentToDiscreteOperator(int scatt_order, bool oneD) override;
    };
}

#endif //CHITECH_CUSTOM_QUADRATURE_H
