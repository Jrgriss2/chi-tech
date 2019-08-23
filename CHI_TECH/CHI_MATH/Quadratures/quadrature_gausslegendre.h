#ifndef _quadrature_gausslegendre_h
#define _quadrature_gausslegendre_h

#include "quadrature.h"
#include <stdio.h>

//######################################################### Class Def
/**Gauss-Legendre quadrature.*/
class CHI_QUADRATURE_GAUSSLEGENDRE : public CHI_QUADRATURE
{
public:
  //01
  void Initialize(int N, int maxiters=1000,
                  double tol=1.0e-10, bool verbose=false);

};

#endif