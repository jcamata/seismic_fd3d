
static char help[] = "SEISMIC 3D.\n\n";

#include <petscdm.h>
#include <petscdmda.h>

#include "seismic.h"



void CreateGrid(FD3D_Parameters *p, int nx, int ny, int nz)
{

   DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
                       nx,ny,nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,HALF_LENGTH,NULL,NULL,NULL,&p->da);
  
   DMGetGlobalVector(p->da,&p->next);
   DMGetGlobalVector(p->da,&p->prev);
   DMGetGlobalVector(p->da,&p->vel);
   p->file_id       = 0;
  
}

void DestroyGrid(FD3D_Parameters *p)
{
    DMDestroy(&p->da)   ;
}


void InitializeArrays(FD3D_Parameters *p)
{
    PetscInt       s,i,j,k,nx,ny,nz,xs,ys,zs,xm,ym,zm;
    VecSet(p->next, 0.0);
    VecSet(p->prev, 0.0);
    VecSet(p->vel , 2250000.0*p->delta_t*p->delta_t);
 
    
    PetscScalar  ***p_prev;
    
    DMDAVecGetArray(p->da,p->prev,&p_prev);
    
    DMDAGetInfo(p->da,PETSC_IGNORE,&nx,&ny,&nz,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
    DMDAGetCorners(p->da,&xs,&ys,&zs,&xm,&ym,&zm);
    
    int nz2 = nz/2;
    int ny4 = ny/2;
    int nx4 = nx/2;
    double val = 1.0;
    
    for(s=5; s >=0; s--) 
    {
        for (k=zs; k<zs+zm; k++) 
        {
            if( (k >= (nz2-s)) && (k < (nz2+s)) ) 
            { 
                for (j=ys; j<ys+ym; j++) {
                    if((j >= (ny4-s)) && (j < (ny4+s)) ) { 
                        for (i=xs; i<xs+xm; i++) {
                            if( (i >= (nx4-s)) && (i < (nx4+s)) ) {
                                p_prev[k][j][i] = val;
                            }
                        } 
                    }
                }
            }
        }
        val*=10.0;
    }
    DMDAVecRestoreArray(p->da,p->prev,&p_prev);

}

void iso_3dfd_it(FD3D_Parameters *p, Vec next, Vec prev, Vec vel, double* coeff)
{
    PetscInt       r,i,j,k,nx,ny,nz,xs,ys,zs,xm,ym,zm;
    Vec   l_prev;
    
    PetscScalar  ***p_vel;
    PetscScalar  ***p_next;
    PetscScalar  ***p_prev;
    PetscScalar  uh2 = 1.0/(p->delta_xyz*p->delta_xyz);
    PetscInt     i_start, i_end, j_start, j_end, k_start, k_end;
    
    DMDAGetInfo(p->da,0,&nx,&ny,&nz,0,0,0,0,0,0,0,0,0);
    DMDAGetCorners(p->da,&xs,&ys,&zs,&xm,&ym,&zm);
     
    DMGetLocalVector(p->da,&l_prev);
    
    DMGlobalToLocalBegin(p->da,prev,INSERT_VALUES,l_prev);
    DMGlobalToLocalEnd(p->da,prev,INSERT_VALUES,l_prev);
    
    DMDAVecGetArray(p->da,vel   ,&p_vel);
    DMDAVecGetArray(p->da,l_prev,&p_prev);
    DMDAVecGetArray(p->da,next  ,&p_next);
    i_start = xs;
    j_start = ys;
    k_start = zs;
    i_end   = xs+xm;
    j_end   = ys+ym;
    k_end   = zs+zm;
    
    if(i_start == 0 )    i_start = HALF_LENGTH;
    if(i_end   == nx)    i_end   = nx-HALF_LENGTH;
    if(j_start == 0 )    j_start = HALF_LENGTH;
    if(j_end   == ny)    j_end   = ny-HALF_LENGTH;
    if(k_start == 0)     k_start = HALF_LENGTH;
    if(k_end == nz )     k_end   = nz-HALF_LENGTH;
       

#pragma omp parallel for private(k,j,i,r)
    for (k=k_start; k<k_end; k++) {
        for (j=j_start; j<j_end; j++) {
            for (i=i_start; i<i_end; i++) {
                //if( i>=HALF_LENGTH && i<(nx-HALF_LENGTH) && j>=HALF_LENGTH && j<(ny-HALF_LENGTH) && k>=HALF_LENGTH && k<(nz-HALF_LENGTH) ) 
                {
                   
                    double value = p_prev[k][j][i]*coeff[0];
                    double FDx = value;
                    double FDy = value;
                    double FDz = value;
                    
                    for(r=1; r <= HALF_LENGTH;r++)
                    {
                        FDx += coeff[r] * (p_prev[k][j][i+r] + p_prev[k][j][i-r]);
                        FDy += coeff[r] * (p_prev[k][j+r][i] + p_prev[k][j-r][i]);
                        FDz += coeff[r] * (p_prev[k+r][j][i] + p_prev[k-r][j][i]);
                    }
                    p_next[k][j][i] = 2.0*p_prev[k][j][i] - p_next[k][j][i] + p_vel[k][j][i]*(uh2*FDx+uh2*FDy+uh2*FDz);
                }
            }
        }
    }
    
    DMDAVecRestoreArray(p->da,vel   ,&p_vel);
    DMDAVecRestoreArray(p->da,next  ,&p_next);
    DMDAVecRestoreArray(p->da,l_prev,&p_prev);
    DMRestoreLocalVector(p->da,&l_prev);    

}


void iso_3dfd( FD3D_Parameters *p, double* coeff)
{
   int it;
   
   VTKIO_PVTI_Write(p, p->prev,"pressure", "teste");
   
   for(it=0; it< p->nreps; it+=2){
       
       iso_3dfd_it(p,p->next,p->prev,p->vel, coeff);
       
       iso_3dfd_it(p,p->prev,p->next,p->vel, coeff);
       
      if(it%50 == 0 ) 
            VTKIO_PVTI_Write(p, p->prev,"pressure", "teste");
       
   }
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{

  FD3D_Parameters p;
  
  double coeff[HALF_LENGTH+1] = {
                        -2.847222222,
                        +1.6,
                        -0.2,
                        +2.53968e-2,
                        -1.785714e-3};
  
  PetscInitialize(&argc,&argv,(char*)0,help);
  
  CreateGrid(&p,256,300,300);
  p.delta_t   = 0.002;
  p.delta_xyz = 50.0;
  p.nreps     = 500;
  
  InitializeArrays(&p);
  
  iso_3dfd(&p,coeff);
  
  //PetscViewer viewer;
  //PetscViewerVTKOpen(PETSC_COMM_WORLD,"fd3d.vts",FILE_MODE_WRITE,&viewer);
  //PetscObjectSetName((PetscObject)p.next,"prev");
  //VecView(p.next,    viewer);
  //PetscViewerDestroy(&viewer);

  //VTKIO_PVTI_Write(&p, p.next,"pressure", "teste");
          
  DestroyGrid(&p);
  PetscFinalize();
  return 0;
}


/*
void WriteSismogram(FD3D_Parameters *p, int xplane, int yplane, const char* filename)
{
    PetscMPIInt    rank, size;
    MPI_Comm       comm;
    PetscInt       r,i,j,k,M,N,P,xs,ys,zs,xm,ym,zm, sum, ierr, pM,pN,pP, stencil;
    const PetscInt *lx,*ly,*lz;
    PetscInt       *osx,*osy,*osz,*oex,*oey,*oez, *olx,*oly,*olz;
    MPI_File fhw;
    
    DMDAGetInfo(p->da,0,&M,&N,&P,&pM,&pN,&pP,0,&stencil,0,0,0,0);
    DMDAGetCorners(p->da,&xs,&ys,&zs,&xm,&ym,&zm);
    ierr = DMDAGetOwnershipRanges(p->da,&lx,&ly,&lz);CHKERRQ(ierr);
    
    PetscObjectGetComm((PetscObject)p->da,&comm);
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);
    
    // generate start,end list 
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

    MPI_File_open(MPI_COMM_WORLD, "sismograma_x.bin",
       MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw);

    if( (xs >= xplane) && (xplane < (xs+xm)) ) {
        
        int local = ym*zm;
        int myoffset = 0;
        for(int i = 0; i <= rank; i++)
        {
            if( osx[i] >= xplane &&  < (oex[i]))
                myoffset += oly[i]*olz[i]*sizeof(double);
        }
        
        MPI_File_write_at(fhw, myoffset, buf, (N/size), MPI_INT, &status);
        
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
    
        
    
  
    
}
*/