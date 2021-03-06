#ifndef SPATIAL_DISCRETIZATION_PWLC_H
#define SPATIAL_DISCRETIZATION_PWLC_H

#include "ChiMath/SpatialDiscretization/CellMappings/FE_PWL/pwl_cellbase.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/spatial_discretization_FE.h"
#include "ChiMath/Quadratures/quadrature_gausslegendre.h"
#include "ChiMath/Quadratures/quadrature_triangle.h"
#include "ChiMath/Quadratures/quadrature_quadrilateral.h"
#include "ChiMath/Quadratures/quadrature_tetrahedron.h"
#include "ChiMath/Quadratures/quadrature_hexahedron.h"

//######################################################### Class def
/**Generalization of the Galerkin Finite Element Method
 * with piecewise linear basis functions
 * for use by either a Continues Finite Element Method (CFEM)
 * or a Discontinuous Finite Element Method (DFEM). */
class SpatialDiscretization_PWLC : public SpatialDiscretization_FE
{
public:
  std::vector<std::shared_ptr<CellMappingFE_PWL>> cell_mappings;

private:
  bool                     mapping_initialized=false;
public:
  chi_math::QuadratureGaussLegendre line_quad_order_second;
  chi_math::QuadratureTriangle      tri_quad_order_second;
  chi_math::QuadratureQuadrilateral quad_quad_order_second;
  chi_math::QuadratureTetrahedron   tet_quad_order_second;
  chi_math::QuadratureHexahedron    hex_quad_order_second;

  chi_math::QuadratureGaussLegendre line_quad_order_arbitrary;
  chi_math::QuadratureTriangle      tri_quad_order_arbitrary;
  chi_math::QuadratureQuadrilateral quad_quad_order_arbitrary;
  chi_math::QuadratureTetrahedron   tet_quad_order_arbitrary;
  chi_math::QuadratureHexahedron    hex_quad_order_arbitrary;

  std::map<int,int> node_mapping;

  int local_block_address = 0;
//  std::vector<int> cell_local_block_address;
//  std::vector<std::pair<int,int>> neighbor_cell_block_address;

  std::vector<int> locJ_block_address;
  std::vector<int> locJ_block_size;

  unsigned int local_base_block_size=0;
  unsigned int globl_base_block_size=0;

private:
//  std::vector<chi_mesh::Cell*> neighbor_cells;
//  std::vector<CellPWLFEValues*> neighbor_cell_fe_views;

private:
  chi_math::finite_element::UnitIntegralData            scratch_intgl_data;
  chi_math::finite_element::InternalQuadraturePointData scratch_vol_qp_data;
  chi_math::finite_element::FaceQuadraturePointData     scratch_face_qp_data;

private:
  //00
  explicit
  SpatialDiscretization_PWLC(chi_mesh::MeshContinuumPtr& in_grid,
                             chi_math::finite_element::SetupFlags setup_flags,
                             chi_math::QuadratureOrder qorder);

public:
  //prevent anything else other than a shared pointer
  static
  std::shared_ptr<SpatialDiscretization_PWLC>
  New(chi_mesh::MeshContinuumPtr& in_grid,
      chi_math::finite_element::SetupFlags setup_flags=
      chi_math::finite_element::SetupFlags::NO_FLAGS_SET,
      chi_math::QuadratureOrder qorder =
      chi_math::QuadratureOrder::SECOND)
  { if (in_grid == nullptr) throw std::invalid_argument(
      "Null supplied as grid to SpatialDiscretization_PWLC.");
    return std::shared_ptr<SpatialDiscretization_PWLC>(
      new SpatialDiscretization_PWLC(in_grid,setup_flags,qorder));}

  //01
private:
  std::shared_ptr<CellMappingFE_PWL> MakeCellMappingFE(const chi_mesh::Cell& cell) const;

public:

  void PreComputeCellSDValues() override;
//  void PreComputeNeighborCellSDValues(chi_mesh::MeshContinuumPtr grid);
  std::shared_ptr<CellMappingFE_PWL> GetCellMappingFE(int cell_local_index);

private:
  //02
  void OrderNodes();

public:
  //03
  void BuildSparsityPattern(chi_mesh::MeshContinuumPtr grid,
                            std::vector<int>& nodal_nnz_in_diag,
                            std::vector<int>& nodal_nnz_off_diag,
                            chi_math::UnknownManager& unknown_manager) override;
//  chi_mesh::Cell* MapNeighborCell(int cell_glob_index);
//  CellPWLFEValues* MapNeighborCellFeView(int cell_glob_index);

  //04 Mappings
  int MapDOF(const chi_mesh::Cell& cell,
             int node,
             const chi_math::UnknownManager& unknown_manager,
             unsigned int unknown_id,
             unsigned int component= 0) const;
  int MapDOFLocal(const chi_mesh::Cell& cell,
                  int node,
                  const chi_math::UnknownManager& unknown_manager,
                  unsigned int unknown_id,
                  unsigned int component= 0) const;
  int MapDOF(const chi_mesh::Cell& cell, int node) const
  { return MapDOF(cell,node,ChiMath::UNITARY_UNKNOWN_MANAGER,0,0); }
//  int MapDOFLocal(const chi_mesh::Cell& cell, int node)
//  { return MapDOFLocal(cell,node,ChiMath::UNITARY_UNKNOWN_MANAGER,0,0); }
  //05
  size_t GetNumLocalDOFs(chi_mesh::MeshContinuumPtr grid,
                         chi_math::UnknownManager& unknown_manager) override;
  size_t GetNumGlobalDOFs(chi_mesh::MeshContinuumPtr grid,
                          chi_math::UnknownManager& unknown_manager) override;
//  unsigned int GetNumGhostDOFs(chi_mesh::MeshContinuumPtr grid,
//                               chi_math::UnknownManager* unknown_manager);
//
//  std::vector<int> GetGhostDOFIndices(chi_mesh::MeshContinuumPtr grid,
//                                      chi_math::UnknownManager* unknown_manager,
//                                      unsigned int unknown_id=0);

  void LocalizePETScVector(Vec petsc_vector,
                           std::vector<double>& local_vector,
                           chi_math::UnknownManager& unknown_manager)
                           override;

  //FE-utils
  const chi_math::finite_element::UnitIntegralData&
  GetUnitIntegrals(const chi_mesh::Cell& cell) override
  {
    if (integral_data_initialized)
      return fe_unit_integrals.at(cell.local_id);
    else
    {
      auto cell_fe_view = GetCellMappingFE(cell.local_id);
      scratch_intgl_data.Reset();
      cell_fe_view->ComputeUnitIntegrals(scratch_intgl_data);
      return scratch_intgl_data;
    }
  }

  const chi_math::finite_element::InternalQuadraturePointData&
  GetQPData_Volumetric(const chi_mesh::Cell& cell) override
  {
    if (qp_data_initialized)
      return fe_vol_qp_data.at(cell.local_id);
    else
    {
      auto cell_fe_view = GetCellMappingFE(cell.local_id);
      cell_fe_view->InitializeVolumeQuadraturePointData(scratch_vol_qp_data);
      return scratch_vol_qp_data;
    }
  }

  const chi_math::finite_element::FaceQuadraturePointData&
  GetQPData_Surface(const chi_mesh::Cell& cell,
                    const unsigned int face) override
  {
    if (qp_data_initialized)
    {
      const auto& face_data = fe_srf_qp_data.at(cell.local_id);

      return face_data.at(face);
    }
    else
    {
      auto cell_fe_view = GetCellMappingFE(cell.local_id);
      cell_fe_view->InitializeFaceQuadraturePointData(face, scratch_face_qp_data);
      return scratch_face_qp_data;
    }

  }
};

#endif //SPATIAL_DISCRETIZATION_PWLC_H