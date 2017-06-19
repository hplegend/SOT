#include <vector>

//Bit shift operator for V3f type
inline std::ostream& operator<< (std::ostream& out, V3f v){
  out << "(" << v.x << "," << v.y << "," << v.z << ")";
  return out;
}

//Declare the CylData class
class CylData {
public:
  CylData(): age(0),order(0),length(0.0),radius(0.0),start(0.0,0.0,0.0),\
	     end(0.0,0.0,0.0),axis(0.0,0.0,0.0),\
	     parent(-1),extension(-1),children(){}
  void setAgeOrder(int A, int O){//Sets age and order info
    age = A;
    order = O;
    is_deleted = false;
    is_dead = false;
  }
  void setLength(float L){//Sets length
    length = L;
  }
  void setRadius(float R){//Sets radius
    radius = R;
  }
  void setAxis(V3f sta, V3f Ax){//Sets axis of the cyl
    start = sta;
    axis = Ax;
    end = start + length*axis;
  }
  void setParent(int p){//Sets the parent for the cyl
    parent = p;
  }
  void setChild(int new_child){//Adds a child of the cyl
    children.push_back(new_child);//
  }
  void setExtension(int ext){//Set the extension of the cyl
    extension = ext;
  }
  int age;// cyl age
  int order;// order of the branch cyl belongs to
  float length;// length of the cyl
  float radius;// radius of the cyl
  V3f start;// starting point
  V3f end;// ending point
  V3f axis;// axis, namely, <V3f>start - <V3f>end (a bit redundant)
  int parent;// cyl's parent, i.e. a cyl it connects to directly
  int extension;//extension of the cyl
  std::vector<int> children;//children of the cyl
  bool is_deleted;//is this cyl deleted?
  bool is_dead;//is dead, i.e. included in the final structure
  float Wf0;// initial foliage mass
  float Wf;// foliage mass of the cyl
  float hwR;// radius of the heartwood of the cyl
};
std::ostream& operator<< (std::ostream& out, CylData cyl){
  out << cyl.age << " ";
  out << cyl.order << " ";
  out << cyl.length << " ";
  out << cyl.radius << " ";
  out << cyl.start << " ";
  out << cyl.end << " ";
  out << cyl.axis << " ";
  out << cyl.Wf0 << " ";
  out << cyl.Wf << " ";
  out << cyl.hwR << " ";
  out << cyl.parent << " ";
  out << cyl.extension << " ";
  out << "(";
  for(int i=0;i<cyl.children.size();i++){
    out << cyl.children[i];
    if(i < cyl.children.size()-1)
      out << ",";
  }
  out << ") ";
  out << cyl.is_deleted;
  out << std::endl;
  return out;
}

//Declare the Branch class
class Branch{
 public:
  std::vector<int> cyl_ind;
  int order;
  int parent_br;
};

std::ostream& operator<< (std::ostream& out, Branch br){
  out << "Order: " << br.order << ",  ";
  out << "Parent: " << br.parent_br << ", ";
  out << "Cyls: (";
  for(int i=0;i<br.cyl_ind.size();i++){
    out << br.cyl_ind[i];
    if(i < br.cyl_ind.size()-1)
      out << ",";
  }
  out << ") ";
  out << std::endl;
  return out;
}
