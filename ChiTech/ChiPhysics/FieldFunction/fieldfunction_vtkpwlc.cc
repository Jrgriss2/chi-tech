#include "fieldfunction.h"

#include "ChiMesh/Cell/cell_slab.h"
#include "ChiMesh/Cell/cell_polygon.h"
#include "ChiMesh/Cell/cell_polyhedron.h"
#include "ChiPhysics/chi_physics.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"

#include "chi_log.h"
#include "chi_mpi.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;
extern ChiPhysics&  chi_physics_handler;

#include <vtkCellType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>

#include <vtkInformation.h>



//###################################################################
/**Handles the PWLD version of a field function export to VTK.
 *
 * */
void chi_physics::FieldFunction::ExportToVTKPWLC(const std::string& base_name,
                                                 const std::string& field_name)
{
  if (spatial_discretization->type !=
      chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_CONTINUOUS)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Field function spatial discretization"
                                " is not of type "
                                " PIECEWISE_LINEAR_CONTINUOUS.");

  auto& pwl_sdm = static_cast<SpatialDiscretization_PWLC&>(*spatial_discretization);

  std::vector<std::vector<double>>    d_nodes;

  auto points = vtkSmartPointer<vtkPoints>::New();

  //============================================= Init grid and material name
  vtkUnstructuredGrid* ugrid;
  vtkIntArray*      matarray;
  vtkIntArray*      pararray;
  vtkDoubleArray*   phiarray;
  vtkDoubleArray*   phiavgarray;

  ugrid    = vtkUnstructuredGrid::New();
  matarray = vtkIntArray::New();
  matarray->SetName("Material");
  pararray = vtkIntArray::New();
  pararray->SetName("Partition");
  phiarray = vtkDoubleArray::New();
  phiarray->SetName(field_name.c_str());
  phiavgarray = vtkDoubleArray::New();
  phiavgarray->SetName((field_name + std::string("-Avg")).c_str());

  //======================================== Precreate nodes to map
  std::vector<std::tuple<uint64_t,uint,uint>> cell_node_component_tuples;

  for (const auto& cell : grid->local_cells)
  {
    auto cell_mapping = pwl_sdm.GetCellMappingFE(cell.local_id);
    for (unsigned int i=0; i<cell_mapping->num_nodes; ++i)
      cell_node_component_tuples.emplace_back(cell.local_id,i,0);
  }

  std::vector<uint64_t> mapping;
  Vec phi_vec;

  CreateCFEMMappingLocal(phi_vec,
                         cell_node_component_tuples,
                         mapping);


  //======================================== Populate cell information
  int nc=0;
  int counter=-1;
  for (const auto& cell : grid->local_cells)
  {
    int mat_id = cell.material_id;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      auto& slab_cell = cell;

      int num_verts = 2;
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        int vgi = slab_cell.vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi]->x;
        d_node[1] = grid->vertices[vgi]->y;
        d_node[2] = grid->vertices[vgi]->z;


        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      ugrid->
        InsertNextCell(VTK_LINE,2,
                       cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        counter++;
        int ir = mapping[counter];
        double dof_value = 0.0;
        VecGetValues(phi_vec,1,&ir,&dof_value);;
        cell_avg_value+= dof_value;
        phiarray->InsertNextValue(dof_value);
      }
      phiavgarray->InsertNextValue(cell_avg_value/num_verts);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      auto& poly_cell = cell;

      int num_verts = poly_cell.vertex_ids.size();
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        int vgi = poly_cell.vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi]->x;
        d_node[1] = grid->vertices[vgi]->y;
        d_node[2] = grid->vertices[vgi]->z;

        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      ugrid->
        InsertNextCell(VTK_POLYGON,num_verts,
                       cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        counter++;
        int ir = mapping[counter];
        double dof_value = 0.0;
        VecGetValues(phi_vec,1,&ir,&dof_value);;
        cell_avg_value+= dof_value;
        phiarray->InsertNextValue(dof_value);
      }
      phiavgarray->InsertNextValue(cell_avg_value/num_verts);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto& polyh_cell = cell;
      auto cell_fe_view = pwl_sdm.GetCellMappingFE(cell.local_id);

      int num_verts = polyh_cell.vertex_ids.size();
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        int vgi = polyh_cell.vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi]->x;
        d_node[1] = grid->vertices[vgi]->y;
        d_node[2] = grid->vertices[vgi]->z;

        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      vtkSmartPointer<vtkCellArray> faces =
        vtkSmartPointer<vtkCellArray>::New();

      int num_faces = polyh_cell.faces.size();
      for (int f=0; f<num_faces; f++)
      {
        int num_fverts = polyh_cell.faces[f].vertex_ids.size();
        std::vector<vtkIdType> face(num_fverts);
        for (int fv=0; fv<num_fverts; fv++)
        {
          int v = cell_fe_view->face_dof_mappings[f][fv];
          face[fv] = cell_info[v];
        }


        faces->InsertNextCell(num_fverts,face.data());
      }//for f

      ugrid->
        InsertNextCell(VTK_POLYHEDRON,num_verts,
                       cell_info.data(),num_faces,faces->GetPointer());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      double cell_avg_value = 0.0;
      for (int v=0; v<num_verts; v++)
      {
        counter++;
        int ir = mapping[counter];
        double dof_value = 0.0;
        VecGetValues(phi_vec,1,&ir,&dof_value);;
        cell_avg_value+= dof_value;
        phiarray->InsertNextValue(dof_value);
      }
      phiavgarray->InsertNextValue(cell_avg_value/num_verts);
    }//polyhedron
  }//for local cells

  ugrid->SetPoints(points);

  //============================================= Construct file name
  std::string base_filename     = std::string(base_name);
  std::string location_filename = base_filename +
                                  std::string("_") +
                                  std::to_string(chi_mpi.location_id) +
                                  std::string(".vtu");

  //============================================= Serial Output each piece
  vtkXMLUnstructuredGridWriter* grid_writer =
    vtkXMLUnstructuredGridWriter::New();

  ugrid->GetCellData()->AddArray(matarray);
  ugrid->GetCellData()->AddArray(pararray);
  ugrid->GetPointData()->AddArray(phiarray);
  ugrid->GetCellData()->AddArray(phiavgarray);

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  grid_writer->Write();

  //============================================= Parallel summary file
  if (chi_mpi.location_id == 0)
  {
      WritePVTU(base_filename, field_name);
  }
}



//###################################################################
/**Handles the PWLD version of a field function export to VTK.
 *
 * */
void chi_physics::FieldFunction::ExportToVTKPWLCG(const std::string& base_name,
                                                  const std::string& field_name)
{
  if (spatial_discretization->type !=
      chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_CONTINUOUS)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Field function spatial discretization"
                                " is not of type "
                                " PIECEWISE_LINEAR_CONTINUOUS.");

  auto& pwl_sdm = static_cast<SpatialDiscretization_PWLC&>(*spatial_discretization);

  std::vector<std::vector<double>>    d_nodes;

  auto points = vtkSmartPointer<vtkPoints>::New();

  auto& ff_uk = this->unknown_manager.unknowns[ref_variable];

  //============================================= Init grid and material name
  vtkUnstructuredGrid* ugrid;
  vtkIntArray*      matarray;
  vtkIntArray*      pararray;
  std::vector<vtkDoubleArray*>   phiarray(ff_uk.num_components);
  std::vector<vtkDoubleArray*>   phiavgarray(ff_uk.num_components);

  ugrid    = vtkUnstructuredGrid::New();
  matarray = vtkIntArray::New();
  matarray->SetName("Material");
  pararray = vtkIntArray::New();
  pararray->SetName("Partition");

  for (int g=0; g < ff_uk.num_components; g++)
  {
    char group_text[100];
    sprintf(group_text,"%03d",g);
    phiarray[g]    = vtkDoubleArray::New();
    phiavgarray[g] = vtkDoubleArray::New();

    phiarray[g]   ->SetName((field_name +
                             std::string("_g") +
                             std::string(group_text)).c_str());
    phiavgarray[g]->SetName((field_name +
                             std::string("_g") +
                             std::string(group_text) +
                             std::string("_avg")).c_str());
  }

  //======================================== Precreate nodes to map

  std::vector<std::tuple<uint64_t,uint,uint>> cell_node_component_tuples;

  for (unsigned int g=0; g < ff_uk.num_components; g++)
    for (const auto& cell : grid->local_cells)
    {
      auto cell_mapping = pwl_sdm.GetCellMappingFE(cell.local_id);
      for (unsigned int i=0; i<cell_mapping->num_nodes; ++i)
        cell_node_component_tuples.emplace_back(cell.local_id,i,g);
    }

  std::vector<uint64_t> mapping;
  Vec phi_vec;

  CreateCFEMMappingLocal(phi_vec,
                         cell_node_component_tuples,
                         mapping);


  //======================================== Populate cell information
  int nc=0;
  int counter=-1;
  for (const auto& cell : grid->local_cells)
  {
    int mat_id = cell.material_id;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SLAB
    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      auto& slab_cell = cell;

      int num_verts = 2;
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        int vgi = slab_cell.vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi]->x;
        d_node[1] = grid->vertices[vgi]->y;
        d_node[2] = grid->vertices[vgi]->z;


        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      ugrid->InsertNextCell(VTK_LINE,2,cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      for (int g=0; g < ff_uk.num_components; g++)
      {
        double cell_avg_value = 0.0;
        for (int v=0; v<num_verts; v++)
        {
          counter++;
          int ir = mapping[counter];
          double dof_value = 0.0;
          VecGetValues(phi_vec,1,&ir,&dof_value);;
          cell_avg_value+= dof_value;
          phiarray[g]->InsertNextValue(dof_value);
        }
        phiavgarray[g]->InsertNextValue(cell_avg_value/num_verts);
      }

    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYGON
    if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      auto& poly_cell = cell;

      int num_verts = poly_cell.vertex_ids.size();
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        int vgi = poly_cell.vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi]->x;
        d_node[1] = grid->vertices[vgi]->y;
        d_node[2] = grid->vertices[vgi]->z;

        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      ugrid->
        InsertNextCell(VTK_POLYGON,num_verts,
                       cell_info.data());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      int old_counter=counter;
      for (int g=0; g < ff_uk.num_components; g++)
      {
        counter = old_counter;
        double cell_avg_value = 0.0;
        for (int v=0; v<num_verts; v++)
        {
          counter++;
          int ir = (int)mapping[counter]+g;
          double dof_value = 0.0;
          VecGetValues(phi_vec,1,&ir,&dof_value);;
          cell_avg_value+= dof_value;
          phiarray[g]->InsertNextValue(dof_value);
        }
        phiavgarray[g]->InsertNextValue(cell_avg_value/num_verts);
      }
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYHEDRON
    if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto& polyh_cell = cell; //(chi_mesh::CellPolyhedron*)(&cell);
      auto cell_fe_view = pwl_sdm.GetCellMappingFE(cell.local_id);

      int num_verts = polyh_cell.vertex_ids.size();
      std::vector<vtkIdType> cell_info(num_verts);
      for (int v=0; v<num_verts; v++)
      {
        int vgi = polyh_cell.vertex_ids[v];
        std::vector<double> d_node(3);
        d_node[0] = grid->vertices[vgi]->x;
        d_node[1] = grid->vertices[vgi]->y;
        d_node[2] = grid->vertices[vgi]->z;

        points->InsertPoint(nc,d_node.data());
        cell_info[v] = nc; nc++;
      }

      vtkSmartPointer<vtkCellArray> faces =
        vtkSmartPointer<vtkCellArray>::New();

      int num_faces = polyh_cell.faces.size();
      for (int f=0; f<num_faces; f++)
      {
        int num_fverts = polyh_cell.faces[f].vertex_ids.size();
        std::vector<vtkIdType> face(num_fverts);
        for (int fv=0; fv<num_fverts; fv++)
        {
          int v = cell_fe_view->face_dof_mappings[f][fv];
          face[fv] = cell_info[v];
        }


        faces->InsertNextCell(num_fverts,face.data());
      }//for f

      ugrid->
        InsertNextCell(VTK_POLYHEDRON,num_verts,
                       cell_info.data(),num_faces,faces->GetPointer());

      matarray->InsertNextValue(mat_id);
      pararray->InsertNextValue(cell.partition_id);

      for (int g=0; g < ff_uk.num_components; g++)
      {
        double cell_avg_value = 0.0;
        for (int v=0; v<num_verts; v++)
        {
          counter++;
          int ir = mapping[counter];
          double dof_value = 0.0;
          VecGetValues(phi_vec,1,&ir,&dof_value);;
          cell_avg_value+= dof_value;
          phiarray[g]->InsertNextValue(dof_value);
        }
        phiavgarray[g]->InsertNextValue(cell_avg_value/num_verts);
      }
    }//polyhedron
  }//for local cells

  ugrid->SetPoints(points);

  //============================================= Construct file name
  std::string base_filename     = std::string(base_name);
  std::string location_filename = base_filename +
                                  std::string("_") +
                                  std::to_string(chi_mpi.location_id) +
                                  std::string(".vtu");

  //============================================= Serial Output each piece
  vtkXMLUnstructuredGridWriter* grid_writer =
    vtkXMLUnstructuredGridWriter::New();

  ugrid->GetCellData()->AddArray(matarray);
  ugrid->GetCellData()->AddArray(pararray);

  for (int g=0; g < ff_uk.num_components; g++)
  {
    ugrid->GetPointData()->AddArray(phiarray[g]);
    ugrid->GetCellData()->AddArray(phiavgarray[g]);
  }

  grid_writer->SetInputData(ugrid);
  grid_writer->SetFileName(location_filename.c_str());

  /*It seems that cluster systems throw an error when the pvtu file
   * also tries to write to the serial file.*/
//  if (chi_mpi.location_id != 0)
    grid_writer->Write();

  //============================================= Parallel summary file
  if (chi_mpi.location_id == 0)
  {
      WritePVTU(base_filename, field_name, (int)ff_uk.num_components);
  }
}
