#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "ChiMesh/chi_mesh.h"

#include <vector>

namespace chi_math
{
  enum class QuadratureOrder : int {
    CONSTANT = 0, FIRST = 1, SECOND = 2, THIRD = 3,
    FOURTH = 4, FIFTH = 5, SIXTH = 6, SEVENTH = 7,
    EIGHTH = 8, NINTH = 9, TENTH = 10, ELEVENTH = 11,
    TWELFTH = 12, THIRTEENTH = 13, FOURTEENTH = 14, FIFTEENTH = 15,
    SIXTEENTH = 16, SEVENTEENTH = 17, EIGHTTEENTH = 18, NINETEENTH = 19,
    TWENTIETH = 20, TWENTYFIRST = 21, TWENTYSECOND = 22, TWENTYTHIRD = 23,
    TWENTYFOURTH = 24, TWENTYFIFTH = 25, TWENTYSIXTH = 26, TWENTYSEVENTH = 27,
    TWENTYEIGHTH = 28, TWENTYNINTH = 29, THIRTIETH = 30, THIRTYFIRST = 31,
    THIRTYSECOND = 32, THIRTYTHIRD = 33, THIRTYFOURTH = 34, THIRTYFIFTH = 35,
    THIRTYSIXTH = 36, THIRTYSEVENTH = 37, THIRTYEIGHTH = 38, THIRTYNINTH = 39,
    FORTIETH = 40, FORTYFIRST = 41, FORTYSECOND = 42, FORTYTHIRD = 43,
    INVALID_ORDER
  };
  typedef chi_mesh::Vector3 QuadraturePointXYZ;
  class Quadrature;
}

//######################################################### Class def
/**Parent class for quadratures.*/
class chi_math::Quadrature
{
public:
  const QuadratureOrder order;
  std::vector<chi_math::QuadraturePointXYZ> qpoints;
  std::vector<double> weights;

protected:
  explicit
  Quadrature(QuadratureOrder in_order) : order(in_order) {}

public:
  void Scale(std::pair<double,double> old_range,
             std::pair<double,double> new_range);

};


#endif