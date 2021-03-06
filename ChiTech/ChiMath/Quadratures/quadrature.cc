#include "quadrature.h"

#include <cassert>

/**Scales 1D quadratures from an old range to a new range.*/
void chi_math::Quadrature::Scale(std::pair<double,double> old_range,
                                 std::pair<double,double> new_range)
{
  double h_new = new_range.second - new_range.first;
  double h_old = old_range.second - old_range.first;

  if (h_new<0.0 or h_old<0.0)
    throw std::invalid_argument("Quadrature::"+std::string(__FUNCTION__)+
                                ": called with negative ranges.");

  if (qpoints.empty())
    throw std::invalid_argument("Quadrature::"+std::string(__FUNCTION__)+
                                ": called with no abscissae initialized.");

  double scale_factor = h_new/h_old;

  for (unsigned int i=0; i < qpoints.size(); ++i)
  {
    qpoints[i](0) = new_range.first +
                    (qpoints[i][0] - old_range.first) * scale_factor;

    weights[i] *= scale_factor;
  }
}



