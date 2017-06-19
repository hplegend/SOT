#ifndef FSTREAM_H
#include <fstream>
#define FSTREAM_H 1
#endif


int find_branch_by_order(std::vector<Branch> & br, int order, \
			 std::vector<int> & br_ind)
{
  br_ind.resize(0);
  for(int i = 0; i < br.size(); i++){
    if(br[i].order == order)
      br_ind.push_back(i);
  }
  
  return 0;
}

int parent_branch(std::vector<CylData> & cyls, std::vector<Branch> & br, int cc)
{//Find the parent branch of a branch cotaining cyl CC as the base cyl.
  if(br.size() == 0)
    return -1;
  
  //std::cout << "Got here, branch " << br_ind << std::endl;
  //Parent cyl of the 1st cyl of the branch

  int par_cyl = cyls[cc].parent;
  if(par_cyl == -1)
    return -1;

  for(int i = 0; i < br.size(); i++){
    for(int j = 0; j < br[i].cyl_ind.size(); j++){
      if(par_cyl == br[i].cyl_ind[j])
	return i;
    }
  }
  
  return -2;//should not be possible
}

int scatter_output(std::vector<CylData> & cyls, int nCyl, std::vector<Branch> & br)
{
  std::ofstream scat;
  scat.open("scatter.dat",std::ofstream::out);
  std::vector<int> curr_brs;
  int order = 0;
  float len,rad,ang,gamma,zeta,tot_len,some;
  int n,m,k;
  V3f rax,newZ,newY,newX,newV;
  V3f newVprojXY,newVprojXZ;

  // Find trunk (branch of 0-order)
  find_branch_by_order(br,order,curr_brs);
  //std::cout << "Order: " << order << ", # br: " << curr_brs.size() << std::endl;

  //*****************************************************************
  //                          HEADER
  //*****************************************************************
  //Write out the header information on Branch and Segment based data-sets
  scat << "# Branch: bra az ltot rini lapar" << std::endl;
  scat << "# Segment: rad len gamma zeta" << std::endl;

  while(curr_brs.size() > 0){
    if(order > MAXORDER){
      break;
    }
    scat << "# order " << order << std::endl;
    //**********************************************
    //                    BRANCH
    //**********************************************
    if(order > 0){//Define for non-trunk branches
      for(int i = 0; i < curr_brs.size(); i++){
	n = br[curr_brs[i]].cyl_ind[0];//1st cyl of the br
	m = parent_branch(cyls,br,n);
	//******* BRA: branching angle of a branch
	some = cyls[n].axis*cyls[cyls[n].parent].axis;
	ang = (180.0/M_PI)*acos(some);
	if(isnan(ang)){
	  std::cout << "Bra: checking NaN instance..." << std::endl;
	  // Maybe precision is wrong?
	  if( (some - 1) < 1.0e-04 && (some - 1) > 0){//Effectively 1.0
	    std::cout << "Resolved to 1.0" << std::endl;
	    scat << (180.0/M_PI)*acos(0.99999999999999) << " ";
	  }
	  else if ( (some + 1)  > -1.0e-04 && (some + 1) < 0 ){//Effectively -1.0
	    std::cout << "Resolved to -1.0" << std::endl;
	    scat << (180.0/M_PI)*acos(-0.99999999999999) << " ";
	  }
	  else{//Finally, either out of range [-1,1] or some is NaN
	    std::cerr << "Outside of the range [-1,1] or NaN input to acos: " << std::endl;
	    std::cerr << cyls[n].axis << "; " << cyls[cyls[n].parent].axis << \
	      ": " << cyls[n].axis*cyls[cyls[n].parent].axis << std::endl;
	    scat << "NaN ";
	  }
	}
	else{
	  scat << ang << " ";
	}
	//******* AZ: azimuth angle (around parent segment/cylinder)
	//Rotation axis with respect to X, since this direction does not change when
	//transfer the structure to Matlab
	newX = V3f(1.0,0.0,0.0);//global X-axis
	rax = cyls[cyls[n].parent].axis % newX;//rot. axis
	ang = acos(cyls[cyls[n].parent].axis * newX);//ang to rotate
	//global Y-axis coincides with LPFG's (-Z)-axis
	//global Z-axis coincides with LPFG's Y-axis
	//(see also lpfg_to_normal_orientation(V3f&) in methods.cpp)
	newY = V3f(0.0,0.0,-1.0);//global Y-axis = (-Z)-lpfg
	newZ = V3f(0.0,1.0,0.0);//global Z-axis = Y-lpfg
	//The code below is exactly the same (accounting on the inconsistency in the
	//coordinate systems of LPFG and Matlab) as in GEN_SCATTER() function
	//extracting the QSM properties (see Bayes Forest Matlab interface).
	newY = rotate_v3f(newY,rax,-ang);
	newZ = rotate_v3f(newZ,rax,-ang);
	//Approximation, since the real angle would be from projections of
	//cyls[n].axis onto newY-newZ-plane
	ang = (180/M_PI)*atan2(cyls[n].axis*newZ,cyls[n].axis*newY);
	if(isnan(ang)){
	  std::cerr << "Az: instance of NaN: " << std::endl;
	  std::cerr << cyls[n].axis << "; " << cyls[cyls[n].parent].axis <<\
	    "; " << rax << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << ang << " ";
	}
	//******* LTOT: total length of a branch
	tot_len = 0.0;
	for(int j = 0; j < br[curr_brs[i]].cyl_ind.size(); j++){
	  tot_len += cyls[br[curr_brs[i]].cyl_ind[j]].length;
	}
	if(isnan(tot_len)){
	  std::cerr << "ltot: instance of NaN " << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << tot_len << " ";
	}
	//******* RINI: Base/inital radius of a branch
	rad = cyls[n].radius;
	if(isnan(rad)){
	  std::cerr << "rini: instance of NaN " << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << rad << " ";
	}
	//******* LAPAR: Distance from the parent's base (length along the parent)
	m = parent_branch(cyls,br,n);//parent branch of the branch i
	k = cyls[n].parent;//parent cyl of n
	len = 0.0;
	for(int j = 0; j < br[m].cyl_ind.size(); j++){
	  len += cyls[br[m].cyl_ind[j]].length;
	  if(br[m].cyl_ind[j] == k)
	    break;
	}
	if(isnan(len)){
	  std::cerr << "lapar: instance of NaN " << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << len << " ";
	}

	//******** Trailing newline (NOTE: put new data above this)
	scat << std::endl;
      }
      
      scat << std::endl << std::endl;
    }
    
    //**********************************************
    //                    SEGMENT
    //**********************************************

    for(int i = 0; i < curr_brs.size(); i++){//Over branches
      len = 0.0;
      for(int j = 0; j < br[curr_brs[i]].cyl_ind.size(); j++){//Over segments
	//******* RAD: Radius along a branch
	rad = cyls[br[curr_brs[i]].cyl_ind[j]].radius;
	if(isnan(rad)){
	  std::cerr << "RAD: instance of NaN" << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << rad << " ";
	}
	//******* LEN: length along a branch
	if(isnan(len)){
	  std::cerr << "len: instance of NaN" << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << len << " ";
	}
	len += cyls[br[curr_brs[i]].cyl_ind[j]].length;
	//******* GAMMA and ZETA: horizontal and vertical variations along a branch
	// ** TWO-ROTATION METHOD **
	// 1. Rotate around (Z) axis for X-axis to coincide with the j-1'th cyl(parent)
	// projection onto XY-plane. Rotate the Y-axis correspondingly and
	// do not change Z-axis. Calculate gamma in the newX-newY plane for j'th cyl(child).
	// 2. Rotate around (Y) axis for newX to coincide with the j-1'th cyl. Rotate
	// the Z-axis correspondingly and do not change Y-axis. Calculate zeta in the
	// newX-newZ plane for j'th cyl(child).
	// NOTE: X = X(lpfg), Y = -Z(lpfg), and Z = Y(lpfg)

	if(j == 0){
	  // for the 1st segment make gamma and zeta zeros
	  gamma = 0.0;
	  zeta = 0.0;
	}
	if(j > 0){
	  //For the segments other than the 1st take adjacent pairs
	  //1. Rotation around Z and XY projection angle gamma
	  modify_coord_xy(cyls[br[curr_brs[i]].cyl_ind[j-1]].axis,newX,newY);
	  gamma = (180/M_PI) * \
	    projection_angle(cyls[br[curr_brs[i]].cyl_ind[j]].axis,newX,newY);
	
	  //2. Rotation around Y and XZ projection angle zeta
	  modify_coord_xz(cyls[br[curr_brs[i]].cyl_ind[j-1]].axis,\
			  newX,newZ,(const V3f)newY);
	  zeta = (180/M_PI) * \
	    projection_angle(cyls[br[curr_brs[i]].cyl_ind[j]].axis,newX,newZ);
	}

	if(isnan(gamma)){
	  std::cerr << "gamma: instance of NaN" << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << gamma << " ";
	}
	if(isnan(zeta)){
	  std::cerr << "zeta: instance of NaN" << std::endl;
	  scat << "NaN ";
	}
	else{
	  scat << zeta << " ";
	}


	//******** Trailing newline (NOTE: put new data above this)
	scat << std::endl;
      }
    }
    scat << std::endl << std::endl;
    
    //+++++++++++ FIND BRANCHES OF THE NEXT ORDER ++++++++++++++
    order++;
    find_branch_by_order(br,order,curr_brs);
    //std::cout << "Order: " << order << ", # br: " << curr_brs.size() << std::endl;
  }

  scat.close();
  return 0;
}

int extract_branches(std::vector<CylData> & cyls, int nCyl, std::vector<Branch> & br)
{//Extract branches based on the analytical growth topology (i.e. resulting from
 //lateral and terminal buds expansion/growth) or empirical growth topology 
 //(the thickest is the extension, others children).
  if(nCyl == 0)
    return 1;

  br.resize(0);
  
  int cc = 0;
  std::vector<int> base_cyls(1,cc);//base cyls = bases of the branches

  // Debuggers
  // std::vector<int> num_bra_by_order((int)MAXORDER+1,0);
  // std::vector<int> num_cyl_by_order((int)MAXORDER+1,0);
  // std::vector<int> cyls_of_br;
  // int j;

  Branch br_tmp;//temporal container for the branch 
  int n_br = 0;//number of branches

  //Variables for the thickest pathway algorithm
  int k = -1;
  float rad = -1.0;
  int init_order = 0;

  //Iterate until the base cyl's array is empty
  while(base_cyls.size() > 0){
    cc = base_cyls[0];//Take the first cyl in the base cyl's array
    if( cc == 0 )//Explicitly set the first cyl order.
      cyls[cc].order = init_order;
    if(cyls[cc].order > (int)MAXORDER)//interrupt due to the max order limitation
      break;
    br_tmp.cyl_ind.resize(0);//temporary branch empty
    br_tmp.cyl_ind.push_back(cc);//base(=first) cyl of the next branch put to br_tmp

    if(!((int)EMPTOP)){
      //*** Extraction algorithm I: analytical topology
      while(cyls[cc].extension > -1){//go through extensions of the cyl's
	//Meanwhile, take care of the new base cyl's coming as child'n of the new branch
	for(int i=0;i<cyls[cc].children.size();i++)
	  base_cyls.push_back(cyls[cc].children[i]);
	//identify the extension of the curent cyl and put it into the branch
	cc = cyls[cc].extension;
	br_tmp.cyl_ind.push_back(cc);
      }
      //Reminder of the loop: new base cyl's as child'n of the last cyl added to the new branch
      for(int i=0;i<cyls[cc].children.size();i++)
	base_cyls.push_back(cyls[cc].children[i]);
      //***
    }
    else{
      //std::cout << "===== Empirical topology algorithm: " << EMPTOP << std::endl;
      //*** Extraction algorithm II: empirical topology (thickest pathway)
      //as if the final tree is observed and child'n+ext indicate only physical
      //connections: like it is done in the Quantitative Structure Models.
      //Here, we use an approximation to the quite complex QSM algorithm. Namely, we
      //go the path from a parent to the thickest offset (children+ext's) or an
      //offset with the largest radius.
      //ORDER: keep the order information in the base cyl's .order field, since we
      //determine the order of the branch by the base cyl. NOTE: the further scatter
      //formation is done branch-wise, thus no interest in setting all the cyls'
      //order information.
      //1. Determine the real extension based on the radius
      while(cyls[cc].extension > -1 || cyls[cc].children.size() > 0){//offsets exist
	//k is an index 0 -- ext, 1,2,... -- 1st,2nd etc. child
	//std::cout << cc << ": " << cyls[cc].extension << "," << cyls[cc].children.size() << std::endl;
	if(cyls[cc].extension > -1){
	  rad = cyls[cyls[cc].extension].radius;
	  k = 0;
	}
	for(int i=0; i<cyls[cc].children.size(); i++){
	  if(rad < cyls[cyls[cc].children[i]].radius){
	    rad = cyls[cyls[cc].children[i]].radius;
	    k = i+1;
	  }
	}
	//std::cout << "k is " << k << std::endl;
	//2. Copy the real extension to the new branch and form base cyl's from the
	//real children 
	if(k == 0){
	  //Nominal children are real children
	  for(int i=0;i<cyls[cc].children.size();i++){
	    base_cyls.push_back(cyls[cc].children[i]);
	    //Fix the order of the children = base cyl's
	    cyls[cyls[cc].children[i]].order = cyls[br_tmp.cyl_ind[0]].order + 1;
	  }
	  //Nominal extension is a real extension
	  cc = cyls[cc].extension;
	}
	else if(k > 0){
	  //Nominal extension is a children, if exists
	  if(cyls[cc].extension > -1){
	    base_cyls.push_back(cyls[cc].extension);
	    //1.2.7 SOT release fixes the below line
	    //from "cyls[cyls[cc].extension].order + 1"
	    // to  "cyls[br_tmp.cyl_ind[0]].order + 1"
	    cyls[cyls[cc].extension].order = cyls[br_tmp.cyl_ind[0]].order + 1;
	  }
	  //All nominal children are real children except for the (k-1)'th one.
	  for(int i=0; i<cyls[cc].children.size(); i++)
	    if(i != k-1){
	      base_cyls.push_back(cyls[cc].children[i]);
	      //Fix the order of the children = base cyl's
	      cyls[cyls[cc].children[i]].order = cyls[br_tmp.cyl_ind[0]].order + 1;
	    }
	  //Nominal (k-1)'th child is a real extension
	  cc = cyls[cc].children[k-1];
	}
	else{// k == -1; should not be reachable in practice!
	  std::cout << "WARNING: k == -1 code reached!" << std::endl;
	  rad = -1.0;
	  break;
	}
	//Reset k and rad
	k = -1;
	rad = -1.0;
	//Put the new extension to the new branch
	br_tmp.cyl_ind.push_back(cc);
      }
      //***
    }
    
    //Erase the first base cyl as it was processed
    base_cyls.erase(base_cyls.begin());

    //Copy to the branch array
    //determine order and parent by the 1st cyl
    br_tmp.order = cyls[br_tmp.cyl_ind[0]].order;
    br_tmp.parent_br = parent_branch(cyls,br,br_tmp.cyl_ind[0]);
    br.push_back(br_tmp);//copy the branch.
    n_br++;//number of branches increased

    // Debug
    // num_bra_by_order[br_tmp.order] += 1;
    // num_cyl_by_order[br_tmp.order] += br_tmp.cyl_ind.size();
    // j = (int)cyls_of_br.size();
    // cyls_of_br.resize(j + br_tmp.cyl_ind.size());
    // for(int i = 0; i < br_tmp.cyl_ind.size(); i++){//look in the current branch
    //   for (int m = 0; m < j; m++){//look in the previously collected cyl's of branches
    // 	if(br_tmp.cyl_ind[i] == cyls_of_br[m]){//Error: already was collected
    // 	  std::cout << "Cyl " << br_tmp.cyl_ind[i] << " was seen as " << \
    // 	    cyls_of_br[m] << std::endl;
    // 	}
    //   }
    //   cyls_of_br[j + i] = br_tmp.cyl_ind[i];
    // }
    //std::cout << br_tmp;
  }

  // std::ofstream temp;
  // temp.open("cyls.dat",std::ofstream::out);
  // j = 0;
  // for(int i = 0; i < br.size(); i++){
  //   if(br[i].order == 2){
  //     std::cout << ++j << " Branch " << i << ": #cyls: " << br[i].cyl_ind.size() << std::endl;
  //     std::cout << br[i];
  //     temp << br[i].cyl_ind.size() << std::endl;
  //   }
  // }
  // for(int i = 0; i <= (int)MAXORDER; i++){
  //   std::cout << "Order: " << i << ", # branches: " << num_bra_by_order[i] << std::endl;
  //   std::cout << "Order: " << i << ", # cyls: " << num_cyl_by_order[i] << std::endl;
  // }

  // temp.close();
  return 0;
}
