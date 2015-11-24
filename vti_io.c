#include <petscdm.h>
#include <petscdmda.h>

typedef struct FD3D_Parameters; 

#undef __FUNCT__
#define __FUNCT__ "VTKIO_vti_appended"
PetscErrorCode VTKIO_vti_appended(FD3D_Parameters *p, DM da,Vec FIELD, const char* field_name, const char file_prefix[])
{
  char           vtk_filename[PETSC_MAX_PATH_LEN];
  PetscMPIInt    rank;
  MPI_Comm       comm;
  FILE           *vtk_fp = NULL;
  PetscInt       si,sj,sk,nx,ny,nz,i;
  PetscInt       f,n_fields,N;
  DM             cda;
  Vec            coords;
  Vec            l_FIELD;
  PetscScalar    *_L_FIELD;
  PetscInt       memory_offset;
  PetscScalar    *buffer;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  /* create file name */
  PetscObjectGetComm((PetscObject)da,&comm);
  MPI_Comm_rank(comm,&rank);
  ierr = PetscSNPrintf(vtk_filename,sizeof(vtk_filename),"subdomain-%s-p%1.4d.vti",file_prefix,rank);CHKERRQ(ierr);

  /* open file and write header */
  vtk_fp = fopen(vtk_filename,"w");
  if (!vtk_fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SYS,"Cannot open file = %s \n",vtk_filename);

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<?xml version=\"1.0\"?>\n");

  /* coords */
  ierr = DMDAGetCorners(da,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
  N    = nx * ny * nz;

#if defined(PETSC_WORDS_BIGENDIAN)
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <ImageData WholeExtent=\"%D %D %D %D %D %D\" Origen=\"0 0 0\" Spacing=%f %f %f\">\n",si,si+nx-1,sj,sj+ny-1,sk,sk+nz-1, p->deltaxyz, p->deltaxyz, p->deltaxyz);
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    <Piece Extent=\"%D %D %D %D %D %D\">\n",si,si+nx-1,sj,sj+ny-1,sk,sk+nz-1);

  memory_offset = 0;

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      <CellData></CellData>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      <PointData Scalars=\" ");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"%s ",field_name);
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"\">\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"        <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%d\"/>\n", field_name,memory_offset);
  memory_offset = memory_offset + sizeof(PetscInt) + sizeof(PetscScalar)*N;

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      </PointData>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    </Piece>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  </ImageData>\n");
  ierr = VecGetArray(FIELD,&_L_FIELD);CHKERRQ(ierr);

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <AppendedData encoding=\"raw\">\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"_");
  int length = sizeof(PetscScalar)*N;
  fwrite(&length,sizeof(int),1,vtk_fp);
  fwrite(_L_FIELD,sizeof(PetscScalar),N,vtk_fp);
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"\n  </AppendedData>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"</VTKFile>\n");
  ierr = VecRestoreArray(FIELD,&_L_FIELD);CHKERRQ(ierr);
 
  if (vtk_fp) {
    fclose(vtk_fp);
    vtk_fp = NULL;
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DAView_3DVTK_PStructuredGrid"
PetscErrorCode VTKIO_vti(FD3D_Parameters *p, DM da, const char fieldname[], const char file_prefix[],const char local_file_prefix[])
{
  MPI_Comm       comm;
  PetscMPIInt    size,rank;
  char           vtk_filename[PETSC_MAX_PATH_LEN];
  FILE           *vtk_fp = NULL;
  PetscInt       M,N,P,si,sj,sk,nx,ny,nz;
  PetscInt       i,dofs;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* only master generates this file */
  PetscObjectGetComm((PetscObject)da,&comm);
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  if (rank != 0) PetscFunctionReturn(0);

  /* create file name */
  ierr   = PetscSNPrintf(vtk_filename,sizeof(vtk_filename),"%s.pvti",file_prefix);CHKERRQ(ierr);
  vtk_fp = fopen(vtk_filename,"w");
  if (!vtk_fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SYS,"Cannot open file = %s \n",vtk_filename);

  /* (VTK) generate pvts header */
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<?xml version=\"1.0\"?>\n");

#if defined(PETSC_WORDS_BIGENDIAN)
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif

  /* define size of the nodal mesh based on the cell DM */
  ierr = DMDAGetInfo(da,0,&M,&N,&P,0,0,0,&dofs,0,0,0,0,0);CHKERRQ(ierr);
  ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <PImageData GhostLevel=\"0\" WholeExtent=\"%d %d %d %d %d %d\">\n",0,M-1,0,N-1,0,P-1); /* note overlap = 1 for Q1 */

  /* DUMP THE CELL REFERENCES */
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    <PCellData>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    </PCellData>\n");

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    <PPointData>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      <PDataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\"/>\n",fieldname);
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    </PPointData>\n");

  /* write out the parallel information */
  ierr = VTKIO_write_PieceExtend(vtk_fp,2,da,local_file_prefix);CHKERRQ(ierr);

  /* close the file */
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  </PStructuredGrid>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"</VTKFile>\n");

  if (vtk_fp) {
    fclose(vtk_fp);
    vtk_fp = NULL;
  }
  PetscFunctionReturn(0);
}
