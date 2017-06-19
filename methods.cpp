#define TOL (1e-05)
#define POINTHROW 100 //how many time to throw Poisson rand var

//Functions' declarations
int delete_cyls(std::vector<CylData> &, std::vector<int> &, int);
float tree_height(std::vector<CylData> & , int , float );
float breast_height_diam(std::vector<CylData> & , int , V3f );
V3f lpfg_to_normal_orientation(V3f& );
int modify_coord_xy(V3f, V3f &, V3f &);
int modify_coord_xz(V3f , V3f & , V3f & , const V3f );
float projection_angle(V3f , V3f , V3f );
V3f rotate_v3f(V3f &,V3f& ,float);

//Some inline handy functions
inline float rad_to_deg(float rad){return ((180.0/M_PI)*rad);}
inline float deg_to_rad(float deg){return ((M_PI/180.0)*deg);}
//Set a random direction branching from the vector v, for rotations purposes
inline V3f rnd_dir(V3f v){
  V3f out(v.x+ran(2.0)-1,v.y+ran(2.0)-1,v.z+ran(2.0)-1);
  return (out.Normalize());
}




int get_all_offsets(int base_cyl,std::vector<CylData> & cyls,	\
		    std::vector<int> &offsets)
{//Get all offset of the base cyl (initiator of a branch), including all children's
 //branches. Deleted cyls do not have topology information (i.e. ext's and child'n,
 //so they will not be in the offsets.
  // if((cyls[base_cyl].parent > -1) &&				\
  //    (cyls[cyls[base_cyl].parent].extension == base_cyl)){
  //   std::cerr << "ERROR(all_offsets): not a base cyl" << std::endl;
  //   return 1;
  // }
  if(cyls[base_cyl].is_deleted){
    std::cerr << "ERROR(all_offsets): deleted cyl" << std::endl;
    return 2;
  }
  int cc;
  std::vector<int> bcyls(1,base_cyl);
  while(bcyls.size()){
    cc = bcyls[0];//Rename the first cyl in BCYLS
    offsets.push_back(cc);//Add it to the offsets
    //Get all ext's of the CC
    while(cyls[cc].extension > -1){
      offsets.push_back(cyls[cc].extension);
      //Add base cyls from the children array of CC
      for(int i=0;i<cyls[cc].children.size();i++)
	bcyls.push_back(cyls[cc].children[i]);
      //Change CC
      cc = cyls[cc].extension;
    }
    //For the last extension
    for(int i=0;i<cyls[cc].children.size();i++)
      bcyls.push_back(cyls[cc].children[i]);

    bcyls.erase(bcyls.begin());
  }
  return 0;
}

int delete_cyls(std::vector<CylData> & cyls, std::vector<int> & cyls_to_remove)
{//Remove CYLS from the tree. Removing is done via a special flag in
 //CylData. Additionally, the parent's information of the removed cyl is cleared so
 //that the removed cyl is no longer among extension/children of the parent
  if(cyls_to_remove.size() == 0)
    return 1;
  int i = 0;
  while( 1 ){
    cyls[cyls_to_remove[i]].is_deleted = true;
    //Adjust topological information of the parent
    if(cyls[cyls_to_remove[i]].parent > -1){
      if( cyls[cyls[cyls_to_remove[i]].parent].extension == cyls_to_remove[i] )
	cyls[cyls[cyls_to_remove[i]].parent].extension = -1;//among extensions
      else{//among children
	for(int j = 0; j < cyls[cyls[cyls_to_remove[i]].parent].children.size(); j++){
	  if(cyls[cyls[cyls_to_remove[i]].parent].children[j] == cyls_to_remove[i]){
	    cyls[cyls[cyls_to_remove[i]].parent].children.\
	      erase(cyls[cyls[cyls_to_remove[i]].parent].children.begin()+j);
	    break;
	  }
	}
      }
    }
    //Adjust topological information of the offsets.
    if(cyls[cyls_to_remove[i]].extension > -1){
      cyls[cyls[cyls_to_remove[i]].extension].parent = -1;
    }
    for(int j = 0; j < cyls[cyls_to_remove[i]].children.size(); j++){
      cyls[cyls[cyls_to_remove[i]].children[j]].parent = -1;
    }
    //Adjust topological information of the cyl itself
    cyls[cyls_to_remove[i]].extension = -1;
    cyls[cyls_to_remove[i]].children.resize(0);
    i++;//increment index of the cyls to be removed
    if(i == cyls_to_remove.size())
      break;
  }

  return 0;
}

int delete_cyl(std::vector<CylData> & cyls, int cyl_to_remove)
{//As delete_cyls but for a single cyl to remove
  std::vector<int> offsets(1,cyl_to_remove);
  return delete_cyls(cyls,offsets);
}

int increase_cyls(int & nCyl,std::vector<CylData> & cyls, \
		  int & cyl_resize_count,int chunk_size)
{
  nCyl++;
  if(nCyl >= (cyl_resize_count*chunk_size)){
    //std::cout << "nCyl = " << nCyl << ". Resizing: ";
    cyls.resize((++cyl_resize_count)*chunk_size);
    //std::cout << cyls.size() << std::endl;
  }
  return 0;
}
int num_deleted(int nCyl, std::vector<CylData>& cyls)
{
  int out = 0;
  for(int i=0; i<nCyl; i++)
    if(cyls[i].is_deleted)
      out++;

  return out;
}

// ********************* TREE HEIGHT AND DIAMETER *********************
float tree_height(std::vector<CylData> & cyls, int nCyl, float lowest)
{//Calculate the tree height
  float H = 0.0;
  for(int i=0;i<nCyl;i++){
    if(cyls[i].is_deleted)
      continue;
    if(H < cyls[i].end.y)
      H = cyls[i].end.y;
  }

  return (H-lowest);
}

float breast_height_diam(std::vector<CylData> & cyls, int nCyl, V3f sta)
{//Calculate the breast height (1.3 m) diameter of the tree.
  float D = 0.0;
  for(int i=0; i<nCyl; i++){
    if(cyls[i].is_deleted || cyls[i].order != 0)
      continue;
    if( (fabs(cyls[i].start.y-sta.y) < 1.3) && fabs(cyls[i].end.y-sta.y) >= 1.3 ){
      D = 2 * cyls[i].radius;
      break;
    }
  }
  
  return D;
}

// ******************* 3D ORIENTATION ****************************
V3f lpfg_to_normal_orientation(V3f& vector)
{
  V3f out = V3f(vector.x,-vector.z,vector.y);
  return out;
}

int modify_coord_xy(V3f vector, V3f & newX, V3f & newY)
{//Modifies the original X,Y,Z coord system such that
  // 1. newX = vector's projection onto XY-plane
  // 2. Y-axis = newY is rotated to correspond to the previous change
  // 3. Z is left unchanged.
  // There is a special procedure served when vector coincides with +/-Z direction.
  // The function account on the special LPFG orientation, namely,
  // X = Xlpfg, Y = -Zlpfg, Z = Ylpfg.
  // NOTE: series of modifications are applied to vector, but it is NOT returned.

  //For reference
  //V3f X = V3f(1.0,0.0,0.0);
  V3f Y = V3f(0.0,0.0,-1.0);
  V3f Z = V3f(0.0,1.0,0.0);
  //projection of vector onto XY
  float tmp = vector.y;
  vector.y = 0.0;
  if(vector.Length() < TOL){//Vertical vector = Z/-Z
    //std::cout << "WARNING: VERTICAL VECTOR!" << std::endl;
    newX = vector;
    Z = newX % Y;
    newY = Y;
  }
  else{//Normal procedure
    vector.Normalize();//normalize the non-zero vector
    //new X-axis
    newX = vector;
    //new Y-axis
    newY = Z % newX;//Y = cross product of Z and X
  }
  // ERRORS CHECK
  //Orthogonal check
  if((fabs(newX*newY) > TOL) || (fabs(newX*Z) > TOL)\
     || (fabs(Z*newY) > TOL)){
    std::cout << "Error(modify_coord_xy): Basis is not orthogonal: " \
	      << newX << ", " << newY << ", " << Z << std::endl;
    return 1;
  }
  //Length check
  if((newX.Length() > (1.0 + TOL)) || (newY.Length() > (1.0 + TOL))){
    std::cout << "Error(modify_coord_xy): Basis length is not one." << std::endl;
    return 2;
  }
  return 0;
}

int modify_coord_xz(V3f vector, V3f & newX, V3f & newZ, const V3f newY)
{//Modifies the coordinate system obtained by MODIFY_COORD_XY() function such that
  // 1. newX = vector
  // 2. Z-axis is rotated to correspond to the previous change
  // 3. Y(newY) is left unchanged (see modify_coord_xy)
  // The function account on the special LPFG orientation, namely,
  // X = Xlpfg, Y = -Zlpfg, Z = Ylpfg.
  // NOTE: series of modifications are applied to vector, but it is NOT returned.

  //For reference
  // V3f X = V3f(1.0,0.0,0.0);
  // V3f Y = V3f(0.0,0.0,-1.0);
  // V3f Z = V3f(0.0,1.0,0.0);
  //new X-axis
  newX = vector;
  //new Z-axis
  newZ = newX % newY;
  // ERRORS CHECK
  //Orthogonal check
  if((fabs(newX*newY) > TOL) || (fabs(newX*newZ) > TOL)\
     || (fabs(newZ*newY) > TOL)){
    std::cout << "Error(modify_coord_xz): Basis is not orthogonal: " \
	      << newX << ", " << newY << ", " << newZ << std::endl;
    return 1;
  }
  //Length check
  if((newX.Length() > (1.0 + TOL)) || (newZ.Length() > (1.0 + TOL))){
    std::cout << "Error(modify_coord_xz): Basis length is not one." << std::endl;
  }
  return 0;
}

float projection_angle(V3f vector, V3f ortho1, V3f ortho2)
{//Calculates the angle of the projection of the vector onto the plane determined by
 //ortho1 and ortho2 with the ortho1 direction.

  //Plane's normal
  V3f normal = ortho1 % ortho2;
  //Find the projection by 2 consecutive cross products
  vector = vector % normal;
  vector = normal % vector;
  if(vector.Length() > TOL){//Normalize non-zero vector
    vector.Normalize();
  }
  //Find angle
  return (atan2f(vector * ortho2, vector * ortho1));
}

V3f rotate_v3f(V3f &vector,V3f& rot_axis,float angle)
{
  // Rotate vector around rot_axis for the specified angle(in rad). This is done
  // only when vector operations on <V3f> type are defined (they are defined by
  // default).

  //Cosine and sine of the angle
  float C = cos((double)angle);
  float S = sin((double)angle);
  //Define the rotation matrix by rows (3 in total)
  rot_axis.Normalize();//normalize the rotation matrix
  float X = rot_axis.x;
  float Y = rot_axis.y;
  float Z = rot_axis.z;
  V3f row1(X*X+(1-X*X)*C,			\
	   X*Y*(1-C)-Z*S,\
	   X*Z*(1-C)+Y*S);
  V3f row2(X*Y*(1-C)+Z*S,\
	   Y*Y+(1-Y*Y)*C,\
	   Y*Z*(1-C)-X*S);
  V3f row3(X*Z*(1-C)-Y*S,\
	   Y*Z*(1-C)+X*S,\
	   Z*Z+(1-Z*Z)*C);
  //vector rotated is a result of matrix multiplication of Rotation matrix and the
  //vector that needs to be rotated. Thus, we can define the resulting vector with
  //components obtained by dot products of the corresponding matrix and old vector
  //components.
  V3f out;
  out.x = row1 * vector;
  out.y = row2 * vector;
  out.z = row3 * vector;
  
  return out;
}

