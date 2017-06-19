/************************************************************
SELF ORGANIZING TREE MODELS FOR IMAGE SYNTHESIS
	
Wojciech Palubicki, Kipp Horel, Steven Longay, Adam Runions,
Brendan Lane, Radomir Mech, Przemyslaw Prusinkiewicz

DEMO PROGRAM #1, Self-organizing growth in light environment 
by Kipp Horel, April 29th, 2009.

Copyrights (C) 2009 by the authors.

This original file is further modified for the needs of the
Bayes Forest by Ilya Potapov.
************************************************************/

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
#include "lpfgall.h"
#pragma clang diagnostic pop
#include <algorithm>
#include <vector>
using std::max;
using std::min;
#include <cstdio>
#include <ctime>
#include <climits>
#include "randist.cpp"
#include "quaternion.hpp"	
#include "hlu.hpp"
#include "env.hpp"
#include "parameters.inc"
#include "user.h" // The user specified input, overwrites the parameters.inc
#include "cyl.hpp" //Cylinder and branch definitions: geometry and topology
#include "methods.cpp" //Various technical routines
#include "mtg.cpp" //Multi-Graph file format
#include "scatter.cpp" //Scatter generation
using namespace std;

//Technical Parameters ---------------------
vector<CylData> cyls;
int cyl_resize_count = 0;
int nCyl = 0;
int nDel = 0;
int fatal_error = 0;
std::ofstream screenlog;
clock_t iter_t;
clock_t tot_t;
int nCyl_prev = 0;

//Segment Data------------------------------
struct SegmentData
{
  bool lateral;//we mark the base of each branch with the lateral state
  float vigor;//this resource is distributed from base to tip
  float accumulatedLight;//based on these factors, which are combined totals from tip to base.
  int count;//the number of distal segments
  float demand;//total resource demand (set in DOWN phase, referenced in UP phase) helps us divide resources
  float width;		
  float length;
  int cyl; //reference to the cylinder data

  V3f pos, heading, up, upGC;//current position & heading, up vector, and General Cylinder Corrected Up Vector.

  SegmentData(const V3f& _h=V3f(0,1,0), const V3f& _u=V3f(0,0,1), const V3f&
    _uGC=V3f(0,0,1), const V3f& _p=V3f(0,0,0), const bool& _l=false, const float&
    _ln=NODELENGTH, const int& _c=0)
    : heading(_h), up(_u), upGC(_uGC), pos(_p), lateral(_l), width(INITWIDTH), count(1),
      length(_ln), cyl(_c)
  {}
};

//MODULES----------------------------------

module Segment(SegmentData);			//Segment modules 	make up the visible bulk of the tree
module Seed(bool,float);			//The seed Module 	is at the base, creating the first Segment
module A(V3f,float);				//Apex modules 		are placed at the tip of a new branch, adding segments when resources allow.
module L(V3f,float);				//Lateral modules 	sprout new branches, when resources allow.
module P(float);				//Prune Cap modules	hold the place where a branch has died. deletes anything distal to itself.

// GLOBAL DECLARATIONS--------------------------------------
#ifdef VEDGE
//bool a = printf("VEDGE is defined!\n");
#define VOXDEN (round((float)SPAN/(float)VEDGE))
#endif
#ifdef VEDGEASNODE
//bool d = printf("VEDGEASNODE is defined!\n");
#define VOXDEN (round((float)SPAN/((float)VEDGEASNODE*(float)NODELENGTH)))
#endif

LightModel env(VOXDEN,SPAN);	//Creating an instance of the LightModel class
int PHASE;					//two phases of computation 1:forward 2:backward
Qf qid;						//Identity Quaternion, used in the between function

float age(0), lambda, oneMinusLambda;		//as tree age increases the preferences for lateral or successor will change over time
int iter = 0;

//AXIOM-----------------------------------------
derivation length: (1+2*(MAXITER-1));
Axiom:	  CurrentTexture(0) TextureVCoeff(0.25) SetColor(31)  Seed(true,1.0f);  

//FUNCTIONS-------------------------------------

Qf between(V3f from, V3f to, float strength)		//Returns a quaternion that will perform a partial rotation from one vector to another (strength: 0 - 1) 
{
  Qf full = rotateTo(from,to);			//calculate the full reorientation quaternion, 
  return slrp(qid,full,strength);			//then use spherical linear interpolation with identity to get the partial result. 
}

// Initialize computation -----------------

Start:
{
  age=0;
  if((int)(DEPTH) <= 0 || (float)(CONSTA) < 0.0 || (float)(CONSTB) <= 1.0){
    std::cout << "ERROR: \n";
    std::cout << "Please set: (int)DEPTH > 0 & (float)CONSTA >= 0.0 & (float)CONSTB > 1.0\n";
    fatal_error = 1;
    Stop();
    return;
  }
  env.setParam(SPREAD,DEPTH,CONSTA,CONSTB,RAN,VERT);//configure the light model parameters
  env.initExposureMap();//initialize the exposure values
  qid=axisAngle(V3f(0,1,0),0.0f);//set an identity value for later use
  PHASE=1;
  //CylData allocation
  cyls.resize((++cyl_resize_count)*NCYLINC);
  //Rand generator initialization
  if((int)RNDSEED){
    //std::cout << "Manual rand seed\n";
    sran((int)RNDSEED);
  }
  else{
    //std::cout << "Automatic rand seed\n";
    sran((int)time(NULL));
  }
  screenlog.open("screen.log",std::ofstream::out);
  time_t timestamp;
  time(&timestamp);
  screenlog << ctime(&timestamp) << std::endl;
  tot_t = clock();

}
/*********/StartEach:
//The only function of our starteach is to alternate between phases.
{
  iter_t = clock();
  if (PHASE==1)
    {
      UseGroup(1);	
      Forward();	// the distal phase moves left to right
      PHASE=2;
      iter++;
      screenlog << "Iter " << iter << " ==> " << std::endl;
    }
  else	
    {
      UseGroup(2);	
      Backward();   // the proximal phase moves right to left.
      PHASE=1;
      screenlog << "Iter " << iter << " <== " << std::endl;
      
      //age is incremented at the start of the proximal phase and drives a change in
      //lambda over time.
      if(age<DELAY+TRANSITION) age++;
      if(age>DELAY)
	lambda = ((DELAY+TRANSITION-age)*INITLAM + (age-DELAY)*TARGLAM)/TRANSITION;
      else
	lambda = INITLAM;
      oneMinusLambda = 1-lambda;//one minus lambda must always complement lambda.
    }
  fatal_error = env.world_violated;
  if( ((nCyl - nDel) > (int)MAXLIVECYL) || (nCyl > (int)MAXNCYL) ){
    std::cout << "Too many (live) cyls! Iter: " << iter << std::endl;
      fatal_error = 1;
  }
  if(fatal_error){
    Stop();//stop the simulation
    return;//skip everything in case of the fatal error
  }
}

/*********/	EndEach:
{
  iter_t = clock() - iter_t;
#ifdef TSTOP //TSTOP is defined in minutes
  if( ((float)iter_t/(60*CLOCKS_PER_SEC)) > TSTOP ){
    std::cout << "WARNING: long simulation (" <<  \
      (float)iter_t/(60*CLOCKS_PER_SEC) << " min)! Stop. " << std::endl;
    Stop();
    return;
  }
#endif
  screenlog << "\tIteration took " << ((float)iter_t)/CLOCKS_PER_SEC << " sec" << \
    std::endl;
  screenlog << "\t#cyl = " << nCyl << ", #live = " << nCyl-nDel << ", #added = " \
	    << nCyl - nCyl_prev <<	std::endl;
  nCyl_prev = nCyl;
  if(fatal_error){
    Stop();//stop the simulation
    return;//skip everything in case of the fatal error
  }

  // *************************************************************************** 
  // We want only the final form without animation, so no DisplayFrame() and the
  //anim.a file is corrected correspondingly (IP comment, modification for Bayes
  //Forest).
  // ***************************************************************************
  
  //if (PHASE==2) DisplayFrame();		//since we are using 2 phases per frame we
					//must use the display on request system.
}

/*********/ End:
{
  tot_t = clock() - tot_t;
  screenlog << "Total simulation time: " << ((float)tot_t)/CLOCKS_PER_SEC << " sec" \
	    << std::endl;
  screenlog.close();//close the screen log file
  if(fatal_error){
    cerr << "Error: world is violated!" << endl;
    remove("scatter.dat");
    remove("out.mtg");
    remove("info.dat");
  }
  //Printf("Age is %g\n",age);
  // for(int i=0; i<nCyl; i++)
  //   std::cout << i << ":\n" << cyls[i] << std::endl;
  std::cout << "%%% Simulation finished. Iterations: " << iter << " %%%" << \
    std::endl;
  std::cout << "Simulation time: " << ((float)tot_t)/CLOCKS_PER_SEC << " sec" << std::endl;
  std::cout << "Total # cyl's: " << nCyl << std::endl;
  std::cout <<  "  #deleted cyl's: " << nDel << std::endl;
  std::cout << "Diameter at base: " << 2*cyls[0].radius << std::endl;
  // std::cout << "Breast height(1.3m) diameter: " << breast_height_diam(cyls,nCyl,V3f(0.0,0.0,0.0)) \
  // 	    << std::endl;
  std::cout << "Height: " << tree_height(cyls,nCyl,(float)0.0) << std::endl;
  if(!fatal_error){//Skip writing output when error occurred
    // MTG output
    if(MTGOUT){
      std::cout << "Forming MTG output...";
      mtg_write(cyls,nCyl);//MTG-output
      std::cout << "Done.\n";
    }
    // Scatter generation
    std::vector<Branch> br;
    if(SCAOUT){
      std::cout << "Extracting branches...";
      extract_branches(cyls,nCyl,br);
      std::cout << "Done (" << br.size() << " branches).\n";
      // for(int i = 0; i < br.size(); i++){
      //   std::cout << "Br. #" << i+1 << "\n\t";
      //   std::cout << br[i];
      // }
      std::cout << "Forming scatters...";
      scatter_output(cyls,nCyl,br);
      std::cout << "Done.\n";
    }
    // Overal info output
    if(INFOUT){
      std::cout << "Overall info...";
      std::ofstream info;
      info.open("info.dat",std::ofstream::out);
      info << "Iter\tnCyl\tDiamBase\tnBr" << std::endl;
      if(SCAOUT)
	info << iter << "\t" << nCyl << "\t" << (2*cyls[0].radius) << "\t" \
	     << br.size() << std::endl;
      else
	info << iter << "\t" << nCyl << "\t" << (2*cyls[0].radius) << "\t" \
	     << "NA" << std::endl;
      info.close();
      std::cout << "Done.\n";
    }
  }
  // Log, whether fatal error occurred
  std::cout << "Log file...";
  std::ofstream logfile;
  logfile.open("out.log",std::ofstream::out);
  logfile << fatal_error << std::endl;
  logfile.close();
  std::cout << "Done.\n";
}


// ******* PROXIMAL PHASE *******
// In the Proximal phase of signalling, Light exposure is sampled at the buds and
// accumulates proximally through fast information transfer.
group 2:
//****
// We must provide a version of the transfer routine for each unique context, in
// order for chain to be complete.
//****
A(p,v): { produce A(p,env.getContinuousExposure(p));}
L(p,v): { produce L(p,env.getContinuousExposure(p));}

Segment(pid) >> SB() Segment(cid) EB() Segment(id):
{
  // this segment has 2 living child branches
  // we accumulate and increment the distal segment count value,
  pid.count=id.count+cid.count+1;
  // and the total light collected by each of the child segments.
  pid.accumulatedLight = id.accumulatedLight + cid.accumulatedLight;
  // these factor into a combined demand for resources, which will be used to help divide available resources in the next phase.
  pid.demand = id.count*lambda*pow(id.accumulatedLight/id.count,NONLINEAR) \
    + cid.count*((1-lambda))*pow(cid.accumulatedLight/cid.count,NONLINEAR);
  // width of branches accumulates according to the pipe model.
  pid.width=pow(pow(id.width,PIPEEXP)+pow(cid.width,PIPEEXP),1.0f/PIPEEXP);
  //Set radius of the CYLS
  cyls[pid.cyl].setRadius(pid.width/2);
  produce Segment(pid);
}

Segment(pid) >> SB() Segment(cid) EB() P(w):// this next production occurs if the main branch has been dropped away. 
{
  // it is similar in spirit to the last one, but modified to fit the situation.
  pid.count=cid.count+1;
  pid.accumulatedLight =cid.accumulatedLight;				 
  pid.demand = cid.count*(1-lambda)*pow(cid.accumulatedLight/cid.count,NONLINEAR);
  pid.width=pow(pow(w,PIPEEXP)+pow(cid.width,PIPEEXP),1.0f/PIPEEXP);
  //Set radius of the CYLS
  cyls[pid.cyl].setRadius(pid.width/2);
  produce Segment(pid);
}

Segment(pid) >> SB() P(w) EB() Segment(id):// here it is the lateral that has come off.	
{
  pid.count=id.count+1;
  pid.accumulatedLight = id.accumulatedLight;				
  pid.demand = id.count*lambda*pow(id.accumulatedLight/id.count,NONLINEAR);
  pid.width=pow(pow(w,PIPEEXP)+pow(id.width,PIPEEXP),1.0f/PIPEEXP);
  //Set radius of the CYLS
  cyls[pid.cyl].setRadius(pid.width/2);
  produce Segment(pid);
}

Segment(pid) >> SB() L(p,v) EB() Segment(id):// these ones deal with the segments still bearing buds.
{
  pid.count=id.count+2;
  pid.accumulatedLight = id.accumulatedLight + v;			
  pid.demand = id.count*lambda*pow(id.accumulatedLight/id.count,NONLINEAR)+((1-lambda))*pow(v,NONLINEAR);
  pid.width=id.width;
  //Set radius of the CYLS
  cyls[pid.cyl].setRadius(pid.width/2);
  produce Segment(pid);
}

Segment(pid) >> SB() L(lp,lv) EB() A(ap,av):// a primary apex
{
  pid.count=3;
  pid.accumulatedLight = av + lv;					
  pid.demand = lambda*pow(av,NONLINEAR)+((1-lambda))*pow(lv,NONLINEAR);
  produce Segment(pid);
}

Segment(pid) >> Segment(id):// just a link in the chain
{
  pid.count = id.count+1;
  pid.accumulatedLight = id.accumulatedLight;
  pid.demand = id.count*lambda*pow(id.accumulatedLight/id.count,NONLINEAR);
  pid.width=id.width;
  //CYLS info update
  cyls[pid.cyl].setRadius(pid.width/2);
  produce Segment(pid);
}

Segment(id) :// dead end
{
  id.count=1;
  id.accumulatedLight =0;
  id.demand = 1;
  produce Segment(id);
}

Seed(b,w) >> Segment(id):// base
{
  w=id.width;
  produce Seed(b,w);
}

// ******* DISTAL PHASE *******
// In the distal phase, resources are distributed from the base out to the buds,
// where new growth occurs.
group 1:
// Any branch may be dropped if not enough resources are allocated to it.

Seed(start,w):	//The seed production only activates in the very first step, to create the initial segment.
{
  if(start){//Create first cyl=segment
    cyls[nCyl].setAgeOrder(iter,0);
    cyls[nCyl].setLength((float)NODELENGTH);//See the default init of SegmentData
    cyls[nCyl].setRadius(INITWIDTH/2);
    cyls[nCyl].setAxis(V3f(0.0,0.0,0.0),V3f(0.0,1.0,0.0));
    cyls[nCyl].setParent(-1);
    cyls[nCyl].setExtension(-1);
    increase_cyls(nCyl,cyls,cyl_resize_count,NCYLINC);//nCyl++;
    produce Seed(false,w) \
      Segment( SegmentData(V3f(0,1,0),V3f(0,0,1),V3f(0,0,1),V3f(0,0,0),false,NODELENGTH,nCyl-1) ) \
      SB() L(V3f(0,1,0),0.0) EB() A(V3f(0,1,0),0.0);
  }
}

P(w) << A(p,v): { produce; }// being distal to a P module marks any module for deletion.
P(w) << L(p,v): { produce; }
P(w) << Segment(id):
{
  env.shadow3D(id.pos,false);//for segments the environment shadow should be removed
  //Remove the segment/cyl and all its offsets,i.e. children and extensions (pruning)
  delete_cyl(cyls,id.cyl);
  nDel++;
  produce;
}

//This situation contextually identifies the base segment, located next to the seed.
Seed(b,w) << Segment(id):
{
  // accumulated light is converted into resources with a multiplicative factor.
  id.vigor = ALPHA0+id.accumulatedLight*ALPHA1;
  produce Segment(id);// and the information is saved to the string.
}

// All subsequent segments will take a weighted share of their parent's resources,
Segment(pid) << Segment(id):
{
  // by comparing the parents total resource demand value to their own contribution.
  if(id.lateral)
    id.vigor = \
      pid.vigor * id.count*((1-lambda)*pow(id.accumulatedLight/id.count,NONLINEAR)/pid.demand);
  else
    id.vigor = \
      pid.vigor * id.count*(lambda*pow(id.accumulatedLight/id.count,NONLINEAR)/pid.demand);

  if((id.vigor/(float)id.count) >DROPTHRESH&&id.pos.y>0.0f)
    produce Segment(id);
  else 	{
    // If the resources allocated do not maintain a certain level with respect
    // to the distal count, the segment will die, and be replaced with a P module
    // marking all subsequent segments for deletion.
    env.shadow3D(id.pos,false);
    //P module goes instead of the Segment(id) => we need to remove id.cyl from cyls
    //std::cout << "Delete: " << id.cyl << "=" << cyls[id.cyl] << std::endl;
    delete_cyl(cyls,id.cyl);
    nDel++;
    produce P(id.width);
  } 
}

void grow(float v, V3f h, V3f u, V3f uGC, V3f p, int pcyl, bool l)
{//This loop handles PROLEPTIC growth for both Lateral and Apex modules.
  Qf q;
  int num = v;//number of segments to grow
  //num = 1;//DEBUG
  //length is tweaked so that total shoot length = vigor*NODELENGTH
  float length = (v/(float)num)*NODELENGTH;
  V3f trop_vec = V3f(0.0,TROPISM,0.0);
  V3f env_dir = env.getContinuousEnvDir(p,SAMPLERAD);
  V3f dir(REORIENT*env_dir+trop_vec);
  float strength=min(dir.Length(),(float)REDIRSTREN);
  dir.Normalize();
  int order = cyls[pcyl].order;
  if(l)
    order++;
  // p is the start position of the parent
  // the start position of the offset is, therefore, p + len*h, where
  p += cyls[pcyl].length*cyls[pcyl].axis;

  for(int i=0;i<num;i++)// for 0-n segments grown this season
    {
      //Growth does not take place outside the defined world space
      if(p.x<-(float)SPAN/2 || p.x>(float)SPAN/2 ||
	 p.z<-(float)SPAN/2 || p.z>(float)SPAN/2 ||
	 p.y<0 || p.y>(float)SPAN){
	// std::cout << "ERROR: world is violated (non-fatal): " << "(" << p.x \
	// 	  << "," << p.y << "," << p.z << std::endl;
	break;
      }
      q=between(h,dir,strength);//a partial reorientation quaternion 
      h=q.rotateVector(h);//is applied to the current heading
      uGC=q.rotateVector(uGC);//is applied to the up direction of the cylinder
      u=(axisAngle(h,PHYLROT)*q).rotateVector(u);//along with a PHYLROT (137.52 default) degree spiral phyllotaxis
      env.shadow3D(p,true);//the environment is modified, the START position makes shadow
      //check if the environment was violated, i.e. tree grew outside the given world.
      fatal_error = env.world_violated;
      if(fatal_error)
	return;
      //**** Cyl addition ****
      cyls[nCyl].setAgeOrder(iter,order);
      cyls[nCyl].setLength(length);
      cyls[nCyl].setAxis(p,h);
      cyls[nCyl].setRadius(INITWIDTH/2);
      if(i == 0){
	cyls[nCyl].setParent(pcyl);
	if(l)
	  cyls[pcyl].setChild(nCyl);
	else
	  cyls[pcyl].setExtension(nCyl);
      }
      else{
	cyls[nCyl].setParent(nCyl-1);
	cyls[nCyl-1].setExtension(nCyl);
      }
      increase_cyls(nCyl,cyls,cyl_resize_count,NCYLINC);//nCyl++;
      //************
      //a new segment is put into the string, along with a Lateral bud.
      nproduce Segment(SegmentData(h,u,uGC,p,l,length,nCyl-1)) SB L(p+length*h,0.0) EB;
      p += length*h;//position is modified only after the new segment introduced
      l=false;
    }
  nproduce A(p,v);//an Apex module is appended to the tip of the new shoot.
}

Segment(pid) << L(p,v)://LATERAL BUDS
{
  float vig = min(pid.vigor * (1-lambda)*(pow(v,(float)NONLINEAR)/pid.demand),(float)MAXGROW);
  //std::cout << "L:" << vig << "," << MAXGROW << std::endl;
  if(vig>1.0f)//If resources allow, create a new branch.
    {
      Qf q(axisAngle(pid.up,BRANCHANGLE));// Quaternion for a lateral divergence
      V3f h(q.rotateVector(pid.heading));// is applied to parent heading
      V3f uGC(q.rotateVector(pid.upGC));
      // then that new heading is used to apply a NONPLANTWIST(90) degree twist (NON-PLANARITY)
      V3f u((axisAngle(h,NONPLANTWIST)).rotateVector(pid.up));
      // before letting the standard grow function take over.
      grow(vig,h,u,uGC,pid.pos,pid.cyl,true);
      produce;				
    }
  else{produce L(p,v);}
}

Segment(pid) << A(p,v):
{
  if(SYMPOD)
    produce;//sympodial branching, just remove Apical buds
  // apexes are a bit simpler because we can let our grow() loop do all the work.
  float vig = min(pid.vigor * (lambda*pow(v,(float)NONLINEAR)/pid.demand),(float)MAXGROW);
  //std::cout << "A:" << vig << "," << MAXGROW << std::endl;
  grow(vig,pid.heading,pid.up,pid.upGC,pid.pos,pid.cyl,false); 
  produce;
}

interpretation:		

Seed(b,w):
{
  produce SetWidth(w) ;	//Using the Generalized Cylinder system, starting at the base,	
}

Segment(id)://the tree is displayed by feeding orientation and length into a turtle system.
{//below a certain width threshold a Leaf surface can be appended to each segment.
  V3f lab_pos;
  if(id.width<LEAFINESS){
    produce \
      SetHead(id.heading.x,id.heading.y,id.heading.z,id.upGC.x,id.upGC.y,id.upGC.z) \
      SetWidth(id.width) ContourSides(10) StartGC() G(id.length) EndGC() \
      SB  SetHead(id.heading.x,id.heading.y,id.heading.z,id.up.x,id.up.y,id.up.z) \
      Surface(leaf,0.3f) EB;
  }
  else {
    if(LABEL){
      std::string label;
      label.append(std::to_string(id.cyl));
      //End point: id.length*id.heading + id.pos
      lab_pos = id.pos + (id.length/2)*id.heading;
      nproduce MoveTo(lab_pos.x+id.width,lab_pos.y,lab_pos.z+id.width) \
	SetWidth(id.width) Label(label.c_str());
    }
    produce MoveTo(id.pos.x,id.pos.y,id.pos.z)\
    SetHead(id.heading.x,id.heading.y,id.heading.z,id.upGC.x,id.upGC.y,id.upGC.z) \
      SetWidth(id.width) ContourSides(10) StartGC() G(id.length) EndGC(); \
  }
}

