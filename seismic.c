
static char help[] = "Tests DMGetGlobalVector() and DMRestoreGlobalVector().\n\n";

/*
Use the options
     -da_grid_x <nx> - number of grid points in x direction, if M < 0
     -da_grid_y <ny> - number of grid points in y direction, if N < 0
     -da_processors_x <MX> number of processors in x directio
     -da_processors_y <MY> number of processors in x direction
*/

#include <petscdm.h>
#include <petscdmda.h>

#define HALF_LENGTH 4

typedef struct{
	int nreps;        	// number of time-steps, over which performance is averaged
	Vec prev;	
	Vec next;
	Vec vel;
        DM  da;
        PetscReal delta_t;
        PetscReal delta_xyz;
} FD3D_Parameters; 

void CreateGrid(FD3D_Parameters *p, int nx, int ny, int nz)
{

   DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
                       nx,ny,nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,HALF_LENGTH,NULL,NULL,NULL,&p->da);
  
   DMGetGlobalVector(p->da,&p->next);
   DMGetGlobalVector(p->da,&p->prev);
   DMGetGlobalVector(p->da,&p->vel);
  
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
    int ny4 = ny/4;
    int nx4 = nx/4;
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
    
    DMDAGetInfo(p->da,PETSC_IGNORE,&nx,&ny,&nz,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE,PETSC_IGNORE);
    DMDAGetCorners(p->da,&xs,&ys,&zs,&xm,&ym,&zm);
     
    DMGetLocalVector(p->da,&l_prev);
    
    DMGlobalToLocalBegin(p->da,prev,INSERT_VALUES,l_prev);
    DMGlobalToLocalEnd(p->da,prev,INSERT_VALUES,l_prev);
    
    DMDAVecGetArray(p->da,vel   ,&p_vel);
    DMDAVecGetArray(p->da,l_prev,&p_prev);
    DMDAVecGetArray(p->da,next  ,&p_next);
    
    for (k=zs; k<zs+zm; k++) {
        for (j=ys; j<ys+ym; j++) {
            for (i=xs; i<xs+xm; i++) {
                if( i>=HALF_LENGTH && i<(nx-HALF_LENGTH) && j>=HALF_LENGTH && j<(ny-HALF_LENGTH) && k>=HALF_LENGTH && k<(nz-HALF_LENGTH) ) 
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
   for(it=0; it< p->nreps; it+=2){
       
       iso_3dfd_it(p,p->next,p->prev,p->vel, coeff);
       
       iso_3dfd_it(p,p->prev,p->next,p->vel, coeff);
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
  p.nreps     = 100;
  
  InitializeArrays(&p);
  iso_3dfd(&p,coeff);
  
  PetscViewer viewer;
  PetscViewerVTKOpen(PETSC_COMM_WORLD,"fd3d.vts",FILE_MODE_WRITE,&viewer);
  PetscObjectSetName((PetscObject)p.next,"prev");
  VecView(p.next,    viewer);
  PetscViewerDestroy(&viewer);
          
  DestroyGrid(&p);
  PetscFinalize();
  return;
}
