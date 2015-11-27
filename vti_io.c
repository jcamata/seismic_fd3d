#include <petscdm.h>
#include <petscdmda.h>

#include "seismic.h"

#undef __FUNCT__
#define __FUNCT__ "VTKIO_vti_appended"
PetscErrorCode VTKIO_VTI_appended(FD3D_Parameters *p, DM da,Vec FIELD, const char* field_name, const char file_prefix[])
{
  char           vtk_filename[PETSC_MAX_PATH_LEN];
  PetscMPIInt    rank;
  MPI_Comm       comm;
  FILE           *vtk_fp = NULL;
  PetscInt       si,sj,sk,nx,ny,nz,i;
  PetscInt       N;
  PetscInt       Mx, My, Mz;
  PetscScalar    *_L_FIELD;
  PetscInt       memory_offset;
  PetscScalar    *buffer;
  PetscScalar    x,y,z;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  /* create file name */
  PetscObjectGetComm((PetscObject)da,&comm);
  MPI_Comm_rank(comm,&rank);
  ierr = PetscSNPrintf(vtk_filename,sizeof(vtk_filename),"subdomain-%s-p%1.4d.vti",file_prefix,rank);CHKERRQ(ierr);

  ierr = DMDAGetInfo(da,0,&Mx,&My,&Mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  /* open file and write header */
  vtk_fp = fopen(vtk_filename,"w");
  if (!vtk_fp) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_SYS,"Cannot open file = %s \n",vtk_filename);

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<?xml version=\"1.0\"?>\n");

  /* coords */
  //ierr = DMDAGetCorners(da,&si,&sj,&sk,&nx,&ny,&nz);
  ierr = DMDAGetGhostCorners(da,&si,&sj,&sk,&nx,&ny,&nz);CHKERRQ(ierr);
  N    = nx * ny * nz;

  x = 0.0;
  y = 0.0;
  z = 0.0;

#if defined(PETSC_WORDS_BIGENDIAN)
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"BigEndian\">\n");
#else
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
#endif
  
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <ImageData WholeExtent=\"%D %D %D %D %D %D\" Origen=\"%f %f %f\" Spacing=\"%f %f %f\">\n",si,si+nx-1,sj,sj+ny-1,sk,sk+nz-1, x,y,z, p->delta_xyz, p->delta_xyz, p->delta_xyz);
  //PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <ImageData WholeExtent=\"%D %D %D %D %D %D\" Origen=\"%f %f %f\" Spacing=\"%f %f %f\">\n",0,Mx-1, 0,My-1,0,Mz-1, x,y,z, p->delta_xyz, p->delta_xyz, p->delta_xyz);
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
  PetscMalloc1(N,&buffer);


  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <AppendedData encoding=\"raw\">\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"_");
  int length = sizeof(PetscScalar)*N;
  fwrite(&length,sizeof(int),1,vtk_fp);
  /* load */
  for (i=0; i<N; i++) buffer[i] = _L_FIELD[i];

  fwrite(buffer,sizeof(PetscScalar),N,vtk_fp);
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"\n  </AppendedData>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"</VTKFile>\n");
  ierr = VecRestoreArray(FIELD,&_L_FIELD);CHKERRQ(ierr);

  PetscFree(buffer);
  if (vtk_fp) {
    fclose(vtk_fp);
    vtk_fp = NULL;
  }

  PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "VTKIO_VTI_PieceExtend"
PetscErrorCode VTKIO_VTI_PieceExtend(FILE *vtk_fp,PetscInt indent_level,DM da,const char local_file_prefix[])
{
  PetscMPIInt    size,rank;
  MPI_Comm       comm;
  const PetscInt *lx,*ly,*lz;
  PetscInt       M,N,P,pM,pN,pP,sum,*olx,*oly,*olz;
  PetscInt       *osx,*osy,*osz,*oex,*oey,*oez;
  PetscInt       i,j,k,II,stencil;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* create file name */
  PetscObjectGetComm((PetscObject)da,&comm);
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  ierr = DMDAGetInfo(da,0,&M,&N,&P,&pM,&pN,&pP,0,&stencil,0,0,0,0);CHKERRQ(ierr);
  ierr = DMDAGetOwnershipRanges(da,&lx,&ly,&lz);CHKERRQ(ierr);

  /* generate start,end list */
  ierr = PetscMalloc1(pM+1,&olx);CHKERRQ(ierr);
  ierr = PetscMalloc1(pN+1,&oly);CHKERRQ(ierr);
  ierr = PetscMalloc1(pP+1,&olz);CHKERRQ(ierr);
  sum  = 0;
  for (i=0; i<pM; i++) {
    olx[i] = sum;
    sum    = sum + lx[i];
  }
  olx[pM] = sum;
  sum     = 0;
  for (i=0; i<pN; i++) {
    oly[i] = sum;
    sum    = sum + ly[i];
  }
  oly[pN] = sum;
  sum     = 0;
  for (i=0; i<pP; i++) {
    olz[i] = sum;
    sum    = sum + lz[i];
  }
  olz[pP] = sum;

  ierr = PetscMalloc1(pM,&osx);CHKERRQ(ierr);
  ierr = PetscMalloc1(pN,&osy);CHKERRQ(ierr);
  ierr = PetscMalloc1(pP,&osz);CHKERRQ(ierr);
  ierr = PetscMalloc1(pM,&oex);CHKERRQ(ierr);
  ierr = PetscMalloc1(pN,&oey);CHKERRQ(ierr);
  ierr = PetscMalloc1(pP,&oez);CHKERRQ(ierr);

  for (i=0; i<pM; i++) {
    osx[i] = olx[i] - stencil;
    oex[i] = olx[i] + lx[i] + stencil;
    if (osx[i]<0) osx[i]=0;
    if (oex[i]>M) oex[i]=M;
  }

  for (i=0; i<pN; i++) {
    osy[i] = oly[i] - stencil;
    oey[i] = oly[i] + ly[i] + stencil;
    if (osy[i]<0)osy[i]=0;
    if (oey[i]>M)oey[i]=N;
  }
  for (i=0; i<pP; i++) {
    osz[i] = olz[i] - stencil;
    oez[i] = olz[i] + lz[i] + stencil;
    if (osz[i]<0) osz[i]=0;
    if (oez[i]>P) oez[i]=P;
  }

  for (k=0; k<pP; k++) {
    for (j=0; j<pN; j++) {
      for (i=0; i<pM; i++) {
        char     name[PETSC_MAX_PATH_LEN];
        PetscInt procid = i + j*pM + k*pM*pN; /* convert proc(i,j,k) to pid */
        ierr = PetscSNPrintf(name,sizeof(name),"subdomain-%s-p%1.4d.vti",local_file_prefix,procid);CHKERRQ(ierr);
        for (II=0; II<indent_level; II++) PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  ");

        PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"<Piece Extent=\"%d %d %d %d %d %d\"      Source=\"%s\"/>\n",
                     osx[i],oex[i]-1,
                     osy[j],oey[j]-1,
                     osz[k],oez[k]-1,name);
      }
    }
  }
  ierr = PetscFree(olx);CHKERRQ(ierr);
  ierr = PetscFree(oly);CHKERRQ(ierr);
  ierr = PetscFree(olz);CHKERRQ(ierr);
  ierr = PetscFree(osx);CHKERRQ(ierr);
  ierr = PetscFree(osy);CHKERRQ(ierr);
  ierr = PetscFree(osz);CHKERRQ(ierr);
  ierr = PetscFree(oex);CHKERRQ(ierr);
  ierr = PetscFree(oey);CHKERRQ(ierr);
  ierr = PetscFree(oez);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "VTKIO_VTI_Write"
PetscErrorCode VTKIO_VTI_Write(FD3D_Parameters *p, DM da, const char fieldname[], const char file_prefix[],const char local_file_prefix[])
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
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  <PImageData GhostLevel=\"4\" WholeExtent=\"%d %d %d %d %d %d\" Origin=\"0 0 0\" Spacing=\"%f %f %f\">\n",0,M-1,0,N-1,0,P-1,p->delta_xyz, p->delta_xyz,p->delta_xyz); /* note overlap = 1 for Q1 */

  /* DUMP THE CELL REFERENCES */
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    <PCellData>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    </PCellData>\n");

  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    <PPointData>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"      <PDataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"1\"/>\n",fieldname);
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"    </PPointData>\n");

  /* write out the parallel information */
  ierr = VTKIO_VTI_PieceExtend(vtk_fp,2,da,local_file_prefix);CHKERRQ(ierr);

  /* close the file */
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"  </PImageData>\n");
  PetscFPrintf(PETSC_COMM_SELF,vtk_fp,"</VTKFile>\n");

  if (vtk_fp) {
    fclose(vtk_fp);
    vtk_fp = NULL;
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VTKIO_PVTI_Write"
PetscErrorCode VTKIO_PVTI_Write(FD3D_Parameters *p, Vec x,const char fieldname[], const char filename[])
{
  char           vti_filename[PETSC_MAX_PATH_LEN];
  char           pvti_filename[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PetscSNPrintf(vti_filename,sizeof(vti_filename),"%s.%04d",filename,p->file_id);CHKERRQ(ierr);
  ierr = VTKIO_VTI_appended(p, p->da,x, fieldname, vti_filename);CHKERRQ(ierr);

  ierr = PetscSNPrintf(pvti_filename,sizeof(pvti_filename),"%s.%04d",filename, p->file_id);CHKERRQ(ierr);
  ierr = VTKIO_VTI_Write(p, p->da,fieldname,pvti_filename,vti_filename);CHKERRQ(ierr);
  
  p->file_id++;
  
  PetscFunctionReturn(0);
}

