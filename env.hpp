
#include <math.h>

#ifndef LIGHT_HPP
#define LIGHT_HPP

/*********************************************************************************************
=	LightModel Environment Class			Kipp Horel, Feb 5th 2009
=
=	This class holds a 3 Dimensional voxelized model of approximate Light conditions, represented as floats between 0 and 1.
=
=	Constructor takes 2 arguments: LightModel(# of Divisions (int) , Span (float) ); 
=	the represented space will be centered along X and Z axes, so unless an ODD number of Divisions is given, 
=	the origin will sit on the fence between 2 voxels.
=
=	void initExposureMap( )			- initializes each voxel with a random value between 0.9 and 1, which adds stochasticity and creates non-zero gradients.
=	Float getContinuousExposure( V3f )	- returns exposure value in the voxel which contains or is closest to the given world position.
=	void shadow( V3f, bool )		- when bool=true, reduces light exposure in voxels near the given world position. if bool=false, it reverses this effect. 
=	V3f getContinuousEnvDir( V3f pos, float rad )	
=						- returns 3D vector representing direction and strength of the relative gradient at the given world position.
=
*********************************************************************************************/

class LightModel
{
  int VOXEL_DENSITY;// # of voxels in each dimension, total # of voxels = DENSITY^3.
  float SPAN;// Amount of world space to be represented by the voxel set. size of an individual voxel = SPAN / DENSITY
  float WORLDtoVOX;// this a conversion factor for turning world coordinates into voxel indexes.
  float*** exposure;// Private pointer which will give this class access to it's dynamically created set of voxels.
  float SPREAD;// A higher spread value will produce a wider shadow, and lower = narrower. Must be > 0.
  int DEPTH;// max propagation depth in voxels
  float C1, C2;// propagation constants, C1/(C2+distance) is the shadowing factor. 
  float RAN, VERT;// Constants controlling the random factor and Vertical bias during light field initialization.

  
public:
  bool world_violated;// True = tree grow beyond the SPAN.
  LightModel(int d, float s)
  {							
    VOXEL_DENSITY=d;		
    SPAN=s;
    // this a conversion factor for turning world coordinates into voxel indexes.
    WORLDtoVOX=((float)(VOXEL_DENSITY)/(SPAN));

    exposure=new float** [d];// link main ptr to array of ptr-ptr's.	(one axis down, 2 to go)
    for(int i=0;i<d;i++)
      {
	*(exposure+i)=new float* [d];// link each ptr-ptr to an array of ptr's. (now we have a 2D grid of pointers) 	
	for(int j=0;j<d;j++) exposure[i][j]=new float[d];// link each ptr to an array of floats.	   
      }// (now we have a dynamically allocated 3D array of floats)
  }

  void setParam(float sp, int dp, float ca, float cb, float r, float v)
  {
    SPREAD=sp;// setting the most important model parameters is done with this function.
    DEPTH=dp;// these parameters are described in more depth above.
    C1=ca;
    C2=cb;
    RAN=r;
    VERT=v;
    world_violated = false;
  }

  void initExposureMap()
  {// voxels can be initialized with varying levels of randomization and depth based bias.
    for(int k=0;k<VOXEL_DENSITY;k++)
      {
	for(int i=0;i<VOXEL_DENSITY;i++)
	  {
	    for(int j=0;j+1<VOXEL_DENSITY;j++)
	      //Start from 0 => the "floor" level MUST be initialized in order to be
	      //multiplied to itself (see *= operator below). The bug was "j=1".
	      {
		exposure[i][j][k] = ( 1.0f - (RAN+VERT) ) + ran(RAN) + (VERT*(float)j)/(float)VOXEL_DENSITY;
		if(isnan(exposure[i][j][k]) || fabs(exposure[i][j][k]) > 1000000){
		    std::cout << "E(i,j,k) = " << "E(" << i << "," << j << "," << k << \
		      ")=" << exposure[i][j][k] << std::endl;
		}
	      }
	    //Even remove this unnecessary "lower quality" step. Let the "floor" be a regular level.
	    // exposure[i][0][k] *= 0.7f ;//the "floor", or lower boundary has lower light quality.
	    // if(isnan(exposure[i][0][k]) || fabs(exposure[i][0][k]) > 1000000){
	    //   std::cout << "E(i,0,k) = " << "E(" << i << ",0," << k << \
	    // 	")=" << exposure[i][0][k] << std::endl;
	    // }
	    exposure[i][VOXEL_DENSITY-1][k] = 0.0f ;//the "roof", or upper boundary is completely shaded.
    	    if(isnan(exposure[i][VOXEL_DENSITY-1][k]) ||		\
	       fabs(exposure[i][VOXEL_DENSITY-1][k]) > 1000000){
	      std::cout << "E(i,j,k) = " << "E(" << i << "," << VOXEL_DENSITY-1 << "," << k << \
		")=" << exposure[i][VOXEL_DENSITY-1][k] << std::endl;
	    }
	  }
      }
  }

  void getVoxelIndex(V3f pos, int &x, int &y, int &z)
  {
    // returns the index of the closest exposure voxel. starting point for the shadow
    // algorithm.
    // y is indexed only in the positive domain.  the X and Z dimensions both need an
    // offset, to center the space around the origin.
    y = (int)(pos.y*WORLDtoVOX);
    x = (int)(pos.x*WORLDtoVOX + ((float)VOXEL_DENSITY)/2.0f);
    z = (int)(pos.z*WORLDtoVOX + ((float)VOXEL_DENSITY)/2.0f);

    x = max(min(x,VOXEL_DENSITY-1),0);				
    y = max(min(y,VOXEL_DENSITY-1),0);
    z = max(min(z,VOXEL_DENSITY-1),0);
  }

  void getCornerVoxelIndex(V3f pos, int &x, int &y, int &z)
  {
    //when interpolating, a 2x2x2 cubic set of neighboring voxels is considered.
    //this function is much the same as above, except modifed to return indexes to
    //the top-left-front member of this set. by adding 1 to each dimension in
    //combination,
    //all 8 exposure values will combined to find a weighted average.
    y = (int)(pos.y*WORLDtoVOX);
    x = (int)(pos.x*WORLDtoVOX + ((float)VOXEL_DENSITY-1)/2.0f);
    z = (int)(pos.z*WORLDtoVOX + ((float)VOXEL_DENSITY-1)/2.0f);

    x = max(min(x,VOXEL_DENSITY-2),0);				
    y = max(min(y,VOXEL_DENSITY-2),0);				
    z = max(min(z,VOXEL_DENSITY-2),0);				
  }

  float getContinuousExposure(V3f pos,bool verb=false)
  {// returns exposure value @ world coordinates.	
    int x,y,z;
    // this conversion from world coordinates to voxel indexes enforces array
    // boundaries, and points us to the set of 8 surrounding data points for exposure
    // interpolation.
    getCornerVoxelIndex(pos,x,y,z);
    if(verb){
      std::cout << "Vox: " << x << "," << y << "," << z << std::endl;
    }
    // we convert to voxel space, and remove the Integer portion in each dimension,
    // leaving only a remainder between 0-1. This gives a set of U V W weightings
    // which in combination describe the proximity to each of the 8 voxels.
    float U(pos.x*WORLDtoVOX);		
    U-=floor(U);
    float V(pos.y*WORLDtoVOX);
    V-=floor(V);
    float W(pos.z*WORLDtoVOX);		
    W-=floor(W);
    if(verb){
      std::cout << "U-V-W: " << U << "," << V << "," << W << std::endl;
      std::cout << "pos: " << "(" << pos.x << "," << pos.y << "," << pos.z << ")" << std::endl;
      std::cout << "(x,y,z)=" << "(" << x << "," << y << "," << z << ")\n";
    }
    
    //INTERPOLATION:
    float interpolated=0;
    //if U=1.0, V=1.0, W=1.0, 1.0x1.0x1.0, adds 0.000*exposure[x][y][z] (far corner
    //no effect)
    //if U=0.1, V=0.1, W=0.1, 0.9x0.9x0.9, adds 0.729*exposure[x][y][z] (very close
    //full effect)
    interpolated+=(1.0f-U)*(1.0f-V)*(1.0f-W)*exposure[x][y][z];
    interpolated+=(U)*(1.0f-V)*(1.0f-W)*exposure[x+1][y][z];	
    interpolated+=(1.0f-U)*(V)*(1.0f-W)*exposure[x][y+1][z];
    interpolated+=(U)*(V)*(1.0f-W)*exposure[x+1][y+1][z];
    interpolated+=(1.0f-U)*(1.0f-V)*(W)*exposure[x][y][z+1];
    interpolated+=(U)*(1.0f-V)*(W)*exposure[x+1][y][z+1];		
    interpolated+=(1.0f-U)*(V)*(W)*exposure[x][y+1][z+1];
    interpolated+=(U)*(V)*(W)*exposure[x+1][y+1][z+1];
    //since these 8 combined weights sum to 1.0, by weighted combinination of all 8
    //members, we can interpolate smoothly across the 2x2x2 set.
    if(verb){
      std::cout << "xyz: " << exposure[x][y][z] << std::endl;
      std::cout << "x+1yz: " << exposure[x+1][y][z] << std::endl;
      std::cout << "xy+1z: " << exposure[x][y+1][z] << std::endl;
      std::cout << "xyz+1: " << exposure[x][y][z+1] << std::endl;
      std::cout << "x+1y+1z: " << exposure[x+1][y+1][z] << std::endl;
      std::cout << "x+1yz+1: " << exposure[x+1][y][z+1] << std::endl;
      std::cout << "xy+1z+1: " << exposure[x][y+1][z+1] << std::endl;
      std::cout << "x+1y+1z+1: " << exposure[x+1][y+1][z+1] << std::endl;
      std::cout << "Interp: " << interpolated << std::endl;
    }

    if(fabs(interpolated) > 1000000 || isnan(interpolated)){
      std::cout << "WARNING: exposure interpolated: " << interpolated << std::endl;
      std::cout << "U-V-W: " << U << "," << V << "," << W << std::endl;
      std::cout << "pos: " << "(" << pos.x << "," << pos.y << "," << pos.z << ")" << std::endl;
      std::cout << "(x,y,z)=" << "(" << x << "," << y << "," << z << ")\n";
      std::cout << "xyz: " << exposure[x][y][z] << std::endl;
      std::cout << "x+1yz: " << exposure[x+1][y][z] << std::endl;
      std::cout << "xy+1z: " << exposure[x][y+1][z] << std::endl;
      std::cout << "xyz+1: " << exposure[x][y][z+1] << std::endl;
      std::cout << "x+1y+1z: " << exposure[x+1][y+1][z] << std::endl;
      std::cout << "x+1yz+1: " << exposure[x+1][y][z+1] << std::endl;
      std::cout << "xy+1z+1: " << exposure[x][y+1][z+1] << std::endl;
      std::cout << "x+1y+1z+1: " << exposure[x+1][y+1][z+1] << std::endl;
      std::cout << "Interp: " << interpolated << std::endl;
    }
    return interpolated;
  }


  V3f getContinuousEnvDir(V3f pos, float rad)
  {// Measuring the gradient of the interpolated field
    V3f result;
    //enforce boundary conditions
    pos.x=max(min(pos.x,SPAN*0.9999f),SPAN*-0.9999f);
    pos.y=max(min(pos.y,SPAN*1.9999f),0.0001f);
    pos.z=max(min(pos.z,SPAN*0.9999f),SPAN*-0.9999f);
    // Each dimension is probed a given distance in each direction. The locations for
    // sampling of the environment are determined here, and boundary conditions are
    // enforced.
    float xa(max(pos.x-rad,SPAN*-0.9999f)),  xb(min(pos.x+rad,SPAN*0.9999f));
    float ya(max(pos.y-rad,0.0001f)), yb(min(pos.y+rad,SPAN*1.9999f));
    float za(max(pos.z-rad,SPAN*-0.9999f)), zb(min(pos.z+rad,SPAN*0.9999f));
    // the gradient approximation is calculated by balancing two opposing
    // measurements against each other along each axis.
    result.x+=getContinuousExposure(V3f(xb,pos.y,pos.z));
    result.x-=getContinuousExposure(V3f(xa,pos.y,pos.z));
    result.y+=getContinuousExposure(V3f(pos.x,yb,pos.z));
    result.y-=getContinuousExposure(V3f(pos.x,ya,pos.z));
    result.z+=getContinuousExposure(V3f(pos.x,pos.y,zb));
    result.z-=getContinuousExposure(V3f(pos.x,pos.y,za));
    // if(isinf(result.x) || isinf(result.y) || isinf(result.z)){
    //   std::cout << "za = " << za << ", zb = " << zb << std::endl;
    //   getContinuousExposure(V3f(pos.x,pos.y,za),true);
    //   std::cout << "pos: " << pos.x << "," << pos.y << "," << pos.z << std::endl;
    // }

    return result;
  }


  void shadow3D(V3f pos, bool sub)//This function adds or subtracts shadow from the environment.
  {
    float tmp1;
    float dist;
    int x, y, z, xb, yb, zb;
    getVoxelIndex(pos,x,y,z);// get our voxel indexes
    int max(DEPTH), min(1);
    //DEPTH variation should start with 1 NOT 0. Since the conditions for the factor
    // (1- C1*pow(C2,-dist)) to be nonnegative are:
    // 1. C1 > 0
    // 2. C2 >= 1
    // 3. dist >= 1 (this is made 0 when DEPTH = 0)
    //Additionally, DEPTH > 0 means there is not self-shading.

    if(sub) for(float i=min; i<max; i++)//if sub=true, we will be DECREASING light values.
	      {
		int n(i*SPREAD);
		for(float j=-n;j<=n;j++)			
		  for(float k=-n;k<=n;k++)			 
		    {
		      dist = sqrt(i*i+j*j+k*k);//this is NOT depth (int q) as stated in Palubicki et al. 2009
		      xb=x+j; yb=y-i; zb = z+k;// xb, yb, zb are now the indexes of the voxel being updated;
		      if(xb>=0&&xb<VOXEL_DENSITY&&yb>=0&&yb<VOXEL_DENSITY&&zb>=0&&zb<VOXEL_DENSITY)
			{// if it is within boundaries, modify the exposure value by
			 // a factor determined by distance from x,y,z
			  tmp1 = exposure[xb][yb][zb];
			  //I am not sure of the formula, see below for the increase case
			  //This formula contradicts to the one in Palubicki et al. 2009.
			  exposure[xb][yb][zb]*=( 1 - (C1*pow(C2,-dist)) );
			  //Avoid out-of-range (0.0,1.0) values
			  exposure[xb][yb][zb] = fminf(exposure[xb][yb][zb],(float)1.0);
			  if( exposure[xb][yb][zb] < 0.0 && tmp1 > 0.0){
			    std::cout << "--- DECREASE Negative Exposures ---" << std::endl;
			    std::cout << "Exp new: " << exposure[xb][yb][zb] << ", old: " << tmp1 << std::endl;
			    std::cout << "Pars: " << C1 << "," << C2 << "," << dist << std::endl;
			    std::cout << "i = " << i << ", j = " << j << ", k = " << k << std::endl;
			  }
			}
		    } 
	      }
    else for(float i=min; i<max; i++)// this alternative loop is almost identical, 
	   {// except that it INCREASES light values, negating the effect of the above process.
	     int n(i*SPREAD);
	     for(float j=-n;j<=n;j++)
	       for(float k=-n;k<=n;k++)			 
		 {
		   dist = sqrt(i*i+j*j+k*k);//this is NOT depth (int) as stated in Palubicki et al. 2009
		   xb=x+j; yb=y-i; zb = z+k;							
		   if(xb>=0&&xb<VOXEL_DENSITY&&yb>=0&&yb<VOXEL_DENSITY&&zb>=0&&zb<VOXEL_DENSITY)
		     {
		       tmp1 = exposure[xb][yb][zb];
		       //I am not sure of the formula, this might give the INF values!
		       //This formula contradicts to the one in Palubicki et al. 2009.
		       exposure[xb][yb][zb]/=( 1 - (C1*pow(C2,-dist)) );
		       //Avoid out-of-range (0.0,1.0) values
		       exposure[xb][yb][zb] = fminf(exposure[xb][yb][zb],(float)1.0);
		       if( exposure[xb][yb][zb] < 0.0 && tmp1 > 0.0){
			 std::cout << "--- INCREASE Negative Exposures ---" << std::endl;
		       	 std::cout << "Exp new: " << exposure[xb][yb][zb] << ", old: " << tmp1 << std::endl;
		       	 std::cout << "Pars: " << C1 << "," << C2 << "," << dist << std::endl;
			 std::cout << "i = " << i << ", j = " << j << ", k = " << k << std::endl;			 
		       }
		     } 							
		 }  
	   }
  }

};

#endif // LIGHT_HPP
