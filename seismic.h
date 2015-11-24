
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


