#include "custom_quadrature.h"
#include "ChiLua/chi_lua.h"

#include <fstream>
#include <iomanip>

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiMath/chi_math.h"
extern ChiMath& chi_math_handler;

void chi_math::CustomQuadrature::
MakeHarmonicIndices(int scatt_order, int dimension)
{
  if (m_to_ell_em_map.empty())
  {
    if (dimension == 1)
      for (int ell=0; ell<=scatt_order; ell++)
        m_to_ell_em_map.emplace_back(ell,0);
    else if (dimension == 2)
      for (int ell=0; ell<=scatt_order; ell++)
        if (ell == scatt_order)
        {
          for (int m=-ell; m < 0; m+=2)
          {
            m_to_ell_em_map.emplace_back(ell,m);
          }
        }
        else
        {
          for (int m=-ell; m<=ell; m+=2)
          {
            m_to_ell_em_map.emplace_back(ell,m);
          }
        }
//    else if (dimension == 2)
//      for (int ell=0; ell<=scatt_order; ell++)
//        for (int m=-ell; m<=ell; m+=2)
//        {
//          if (ell == 0 or m != 0)
//            m_to_ell_em_map.emplace_back(ell,m);
//        }
    else if (dimension == 3)
      for (int ell=0; ell<=scatt_order; ell++)
        for (int m=-ell; m<=ell; m++)
          m_to_ell_em_map.emplace_back(ell,m);
  }
}

void chi_math::CustomQuadrature::
  BuildDiscreteToMomentOperator(int scatt_order, int dimension)
{
  chi_log.Log() << "Hello";

//  MakeHarmonicIndices(scatt_order,dimension);
//
//  int num_angles = abscissae.size();
//  int num_moms = 0;
  float a;
//
//  d2m_op.clear();
//
  std::ifstream file;
  file.open(D2M_file);
  if (not file.is_open())
  {
    chi_log.Log(LOG_ALLERROR) << "Failed to open Quadrature File in call to " << __FUNCTION__ << ".";
    exit(EXIT_FAILURE);
  }
//
  file >> a;
  double moments = a;
//
//  int mc=-1; //moment count
//  for (int ell=0; ell < moments; ell++)
//  {
//      std::vector<double> cur_mom; mc++;
//      num_moms++;
//
//      for (int n = 0; n < num_angles; n++)
//      {
//        file >> a;
//        double value = a;
//        cur_mom.push_back(value);
////      chi_log.Log() << value;
//      }
//
//      d2m_op.push_back(cur_mom);
//  }//for ell

  if (d2m_op_built) return;

  d2m_op.clear();
  MakeHarmonicIndices(scatt_order,dimension);

  int num_angles = abscissae.size();
  int num_moms = m_to_ell_em_map.size();

  for (const auto& ell_em : m_to_ell_em_map)
  {
    std::vector<double> cur_mom;
    cur_mom.reserve(num_angles);

    for (int n=0; n<num_angles; n++)
    {
      const auto& cur_angle = abscissae[n];

      file >> a;
      double value = a;
      double w = weights[n];
      cur_mom.push_back(value);
    }

    d2m_op.push_back(cur_mom);
  }
  d2m_op_built = true;

  std::stringstream outs;
  outs
          << "\nQuadrature d2m operator:\n";
  for (int n=0; n<num_angles; n++)
  {
    outs << std::setw(5) << n;
    for (int m=0; m<num_moms; m++)
    {
      outs
              << std::setw(15) << std::left << std::fixed
              << std::setprecision(10) << d2m_op[m][n] << " ";
    }
    outs << "\n";
  }
  chi_log.Log() << outs.str();

  file.close();

  chi_log.Log() << "Goodbye";
}
void chi_math::CustomQuadrature::BuildMomentToDiscreteOperator(int scatt_order, int dimension)
{
    chi_log.Log() << "Hello";

//    MakeHarmonicIndices(scatt_order,dimension);
//
//    int num_angles = abscissae.size();
//    int num_moms = 0;
    float a;
//
//    m2d_op.clear();
//
    std::ifstream file;
    file.open(M2D_file);
    if (not file.is_open())
    {
        chi_log.Log(LOG_ALLERROR) << "Failed to open Quadrature File in call to " << __FUNCTION__ << ".";
        exit(EXIT_FAILURE);
    }

    file >> a;
    double moments = a;
//
//    int mc=-1; //moment count
//    for (int ell=0; ell<moments; ell++)
//    {
//            std::vector<double> cur_mom; mc++;
//            num_moms++;
//
//            for (int n = 0; n < num_angles; n++)
//            {
//                file >> a;
//                double value = a;
//                cur_mom.push_back(value);
////                chi_log.Log() << value;
//            }
//
//            m2d_op.push_back(cur_mom);
//    }//for ell


  if (m2d_op_built) return;

  m2d_op.clear();
  MakeHarmonicIndices(scatt_order,dimension);

  int num_angles = abscissae.size();
  int num_moms = m_to_ell_em_map.size();

  double normalization = 1.0;
  if (dimension == 1) normalization = 2.0;
  if (dimension == 2) normalization = 4.0*M_PI;
  if (dimension == 3) normalization = 4.0*M_PI;

  for (const auto& ell_em : m_to_ell_em_map)
  {
    std::vector<double> cur_mom;
    cur_mom.reserve(num_angles);

    for (int n=0; n<num_angles; n++)
    {
      const auto& cur_angle = abscissae[n];
      file >> a;
      double value = a;
      cur_mom.push_back(value);
    }

    m2d_op.push_back(cur_mom);
  }//for m
  m2d_op_built = true;

    std::stringstream outs;
    outs
            << "\nQuadrature m2d operator:\n";
    for (int n=0; n<num_angles; n++)
    {
        outs << std::setw(5) << n;
        for (int m=0; m<num_moms; m++)
        {
            outs
                    << std::setw(15) << std::left << std::fixed
                    << std::setprecision(10) << m2d_op[m][n] << " ";
        }
        outs << "\n";
    }
    chi_log.Log() << outs.str();

    file.close();

    chi_log.Log() << "Goodbye";
}

// #########################################################################################
// #########################################################################################
// #########################################################################################

int chiCreateCustomQuadrature(lua_State* L)
{

//    printf("Hello World");

//    int num_args = lua_gettop(L);
//    if (num_args != 1)
//        LuaPostArgAmountError(__FUNCTION__ , 1, num_args);

    const char* file_name_raw1 = lua_tostring(L, 1);
    std::string quad_file(file_name_raw1);
    const char* file_name_raw2 = lua_tostring(L, 2);
    std::string M2D_file(file_name_raw2);
    const char* file_name_raw3 = lua_tostring(L, 3);
    std::string D2M_file(file_name_raw3);

    std::ifstream file;
    file.open(quad_file);
    if (not file.is_open())
    {
        chi_log.Log(LOG_ALLERROR) << "Failed to open Quadrature File in call to " << __FUNCTION__ << ".";
        exit(EXIT_FAILURE);
    }

    auto new_quad = std::make_shared<chi_math::CustomQuadrature>(quad_file, M2D_file, D2M_file);
    chi_math::QuadraturePointPhiTheta q_points;

    float a;
    chi_log.Log() << "Quad Points (w, phi, theta):";
    while(file >> a)
    {
        double weight = a;
        file >> a;
        q_points.phi = a;
        file >> a;
        q_points.theta = a;
        chi_log.Log() << weight << "   "  << q_points.phi << "   " << q_points.theta;
//        chi_log.Log() << q_points.theta;
        double x = sin(q_points.theta)*cos(q_points.phi);
        double y = sin(q_points.theta)*sin(q_points.phi);
        double z = cos(q_points.theta);

        new_quad->abscissae.push_back(q_points);
        new_quad->weights.push_back(weight);
        new_quad->omegas.emplace_back(x, y, z);
    }
//    q_points.phi = 0.5;
//    q_points.theta = 0.5;
//    double x = sin(q_points.theta)*cos(q_points.phi);
//    double y = sin(q_points.theta)*sin(q_points.phi);
//    double z = cos(q_points.theta);
//
//    new_quad->abscissae.push_back(q_points);
//    new_quad->weights.push_back(1.0);
//    new_quad->omegas.emplace_back(x, y, z);

//    std::string line;
//    bool not_eof = bool(std::getline(file, line));
//
//    while (not_eof)
//    {
//        chi_log.Log() << line << "\n";
//        not_eof = bool(std::getline(file, line));
//    }

    file.close();

    chi_math_handler.angular_quadratures.push_back(new_quad);
    lua_pushnumber(L, chi_math_handler.angular_quadratures.size()-1);

    return 1;

}