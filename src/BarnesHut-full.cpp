/*
 Lonestar BarnesHut: Simulation of the gravitational forces in a
 galactic cluster using the Barnes-Hut n-body algorithm

 Author: Martin Burtscher
 Center for Grid and Distributed Computing
 The University of Texas at Austin

 Copyright (C) 2007, 2008 The University of Texas at Austin

 Licensed under the Eclipse Public License, Version 1.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.eclipse.org/legal/epl-v10.html

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 File: BarnesHut.cpp
 Modified: Feb. 19, 2008 by Martin Burtscher (initial C++ version)
 Modified: May 26, 2009 by Nicholas Chen (uses TBB constructs)
 */

//#define DEBUG

//#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sys/time.h>
//#include <iostream>

#include "asap.h"
//#include "tbb/blocked_range.h"
//#include "tbb/parallel_for.h"
//#include "tbb/parallel_invoke.h"
#include "../../include/blocked_range.h"
//#include "../../include/parallel_for.h"
//#include "../../include/parallel_invoke.h"

namespace tbb {
    template<typename Func0, typename Func1, typename Func2, typename Func3,
             typename Func4, typename Func5, typename Func6, typename Func7>
    void parallel_invoke[[asap::param("P0, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10, P11, P12, P13, P14, P15")]]
                        (const Func0 &F0 [[asap::arg("P0, P1")]],   const Func1 &F1 [[asap::arg("P2, P3")]], 
                         const Func2 &F2 [[asap::arg("P4, P5")]],   const Func3 &F3 [[asap::arg("P6, P7")]],
                         const Func4 &F4 [[asap::arg("P8, P9")]],   const Func5 &F5 [[asap::arg("P10, P11")]], 
                         const Func6 &F6 [[asap::arg("P12, P13")]], const Func7 &F7 [[asap::arg("P14, P15")]]);

} // end namespace


//#include "tbb/task_scheduler_init.h"
//#include "tbb/task_group.h"
#line 47 "BarnesHut.cpp"

static double dtime; // length of one time step
static double eps; // potential softening parameter
static double tol; // tolerance for stopping recursion, should be less than 0.57 for 3D case to bound error

static double dthf, epssq, itolsq;

static int step;

using namespace tbb;
using namespace std;

enum {
  CELL, BODY
};

class /*PARAM(R)*/ OctTreeLeafNode;

static OctTreeLeafNode **bodies ARG_(Global, Global, Global) ; // the n bodies


class PARAM(PN) OctTreeNode {
public:
  int type ARG_(PN); // CELL or BODY
  double mass ARG_(PN);
  double posx ARG_(PN);
  double posy ARG_(PN);
  double posz ARG_(PN);

  OctTreeNode(int type, double mass, double posx, double posy, double posz) : type(type), mass(mass), posx(posx), posy(posy), posz(posz) {} // expected-warning{{Inferred Effect Summary for OctTreeNode: [reads(rpl([rLOCAL],[]))]}}

  OctTreeNode() : type(0), mass(0), posx(0), posy(0), posz(0) {}

#ifdef DEBUG
  friend std::ostream& operator<<(std::ostream &stream, const OctTreeNode &node);
#endif
}; // end class OctTreeNode

#ifdef DEBUG
std::ostream &operator<<(std::ostream &stream, const OctTreeNode &node) {
  stream << "{posx=" << node.posx << ", posy=" << node.posy << ", posz="
    << node.posz << "}\n" << flush;
  return stream;
}
#endif

// The internal nodes are cells that summarize their children's properties
class
REGION_(
Rc0, Rc1, Rc2, Rc3, Rc4, Rc5, Rc6, Rc7,
Rp0, Rp1, Rp2, Rp3, Rp4, Rp5, Rp6, Rp7)
PARAM(PiN)
BASEARG(OctTreeNode, PiN)
OctTreeInternalNode: public OctTreeNode {
public:

  OctTreeInternalNode(double posx, double posy, double posz) : OctTreeNode(CELL,
      0, posx, posy, posz), child0(0), child1(0), child2(0), child3(0),
  child4(0), child5(0), child6(0), child7(0) { }  // expected-warning{{Inferred Effect Summary for OctTreeInternalNode: [reads(rpl([rLOCAL],[]))]}}

  //static void NewNode PARAM(R) (OctTreeInternalNode * in ARG_(R), double px, double py, double pz);
  void copyFrom PARAM(Rn) WRITES(PiN) READS(Rn) (OctTreeInternalNode *node ARG_(Rn)) { // expected-warning{{Inferred Effect Summary for copyFrom: [reads(rpl([p51_Rn],[])),reads(rpl([rLOCAL],[])),writes(rpl([p50_PiN],[]))]}}
    type = node->type;
    mass = node->mass;
    posx = node->posx;
    posy = node->posy;
    posz = node->posz;
  }
  //static void RecycleTree() {
  //  freelist = head;
  //}

public:
  //void Insert(OctTreeLeafNode * const b, const double r); // builds the tree
  void InsertAll PARAM(Rb) WRITES(Rb:*, PiN:*) (OctTreeLeafNode **b ARG_(Rb, Rb), int n, double r);
  // recursively summarizes info about subtrees
  void ComputeCenterOfMass WRITES(*) (int &curr ARG_(Global));
  OctTreeNode** GetChildRef ARG_(*, *) (int i);
  OctTreeNode* GetChild ARG_(*, *) READS(*) (int i);

  //OctTreeLeafNode **GetPartition(int i, OctTreeLeafNode **partition0,
  //  OctTreeLeafNode **partition1, OctTreeLeafNode **partition2,
  //  OctTreeLeafNode **partition3, OctTreeLeafNode **partition4,
  //  OctTreeLeafNode **partition5, OctTreeLeafNode **partition6,
  //  OctTreeLeafNode **partition7);

  OctTreeLeafNode **GetPartition
  ARG_(Rb:*, Rb:*)
  PARAM(Rb)
  (OctTreeLeafNode **b ARG_(Rb, Rb), int i,
   OctTreeLeafNode **partition0 ARG_(Rb:Rp0, Rb:Rp0),
   OctTreeLeafNode **partition1 ARG_(Rb:Rp1, Rb:Rp1),
   OctTreeLeafNode **partition2 ARG_(Rb:Rp2, Rb:Rp2),
   OctTreeLeafNode **partition3 ARG_(Rb:Rp3, Rb:Rp3),
   OctTreeLeafNode **partition4 ARG_(Rb:Rp4, Rb:Rp4),
   OctTreeLeafNode **partition5 ARG_(Rb:Rp5, Rb:Rp5),
   OctTreeLeafNode **partition6 ARG_(Rb:Rp6, Rb:Rp6),
   OctTreeLeafNode **partition7 ARG_(Rb:Rp7, Rb:Rp7));

  static void ChildIDToPos PARAM(Rx, Ry, Rz) WRITES(Rx, Ry, Rz)
    (int childID, double radius, double *x ARG_(Rx), double *y ARG_(Ry), double *z ARG_(Rz));

  OctTreeNode *child0 ARG_(PiN:Rc0, PiN:Rc0);
  OctTreeNode *child1 ARG_(PiN:Rc1, PiN:Rc1);
  OctTreeNode *child2 ARG_(PiN:Rc2, PiN:Rc2);
  OctTreeNode *child3 ARG_(PiN:Rc3, PiN:Rc3);
  OctTreeNode *child4 ARG_(PiN:Rc4, PiN:Rc4);
  OctTreeNode *child5 ARG_(PiN:Rc5, PiN:Rc5);
  OctTreeNode *child6 ARG_(PiN:Rc6, PiN:Rc6);
  OctTreeNode *child7 ARG_(PiN:Rc7, PiN:Rc7);

private:
  //OctTreeInternalNode *link; // links all internal tree nodes so they can be recycled
  // free list for recycling
  static OctTreeInternalNode *head, *freelist;

  int ChildID PARAM(Rb) READS(PiN, Rb) (OctTreeLeafNode *const b ARG_(Rb));

  void Partition
  PARAM(Rb, Rps) READS(PiN) WRITES(Rps, Rb:*)
  (OctTreeLeafNode **b ARG_(Rb, Rb), int n, int *partitionSize ARG_(Rps),
  OctTreeLeafNode **partition0 ARG_(Rb:Rp0, Rb:Rp0),
  OctTreeLeafNode **partition1 ARG_(Rb:Rp1, Rb:Rp1),
  OctTreeLeafNode **partition2 ARG_(Rb:Rp2, Rb:Rp2),
  OctTreeLeafNode **partition3 ARG_(Rb:Rp3, Rb:Rp3),
  OctTreeLeafNode **partition4 ARG_(Rb:Rp4, Rb:Rp4),
  OctTreeLeafNode **partition5 ARG_(Rb:Rp5, Rb:Rp5),
  OctTreeLeafNode **partition6 ARG_(Rb:Rp6, Rb:Rp6),
  OctTreeLeafNode **partition7 ARG_(Rb:Rp7, Rb:Rp7)); 

  void InsertChildren
  PARAM(Pp0, Pp1, Pp2, Pp3, Pp4, Pp5, Pp6, Pp7, Pps)
  READS(Pps)
  WRITES(PiN:*, Pp0:*, Pp1:*, Pp2:*, Pp3:*, Pp4:*, Pp5:*, Pp6:*, Pp7:*)
  (double r, int *partitionSize ARG_(Pps),
  OctTreeLeafNode **partition0 ARG_(Pp0, Pp0),
  OctTreeLeafNode **partition1 ARG_(Pp1, Pp1),
  OctTreeLeafNode **partition2 ARG_(Pp2, Pp2),
  OctTreeLeafNode **partition3 ARG_(Pp3, Pp3),
  OctTreeLeafNode **partition4 ARG_(Pp4, Pp4),
  OctTreeLeafNode **partition5 ARG_(Pp5, Pp5),
  OctTreeLeafNode **partition6 ARG_(Pp6, Pp6),
  OctTreeLeafNode **partition7 ARG_(Pp7, Pp7));
}; // end class OctTreeInternalNode

// the tree leaves are the bodies
class PARAM(R) BASEARG(OctTreeNode, R)
OctTreeLeafNode: public OctTreeNode {
public:
  WRITES(R) OctTreeLeafNode() {  // expected-warning{{Inferred Effect Summary for OctTreeLeafNode: [writes(rpl([p69_R],[]))]:}}
    type = BODY;
    mass = 0.0;
    posx = 0.0;
    posy = 0.0;
    posz = 0.0;
    velx = 0.0;
    vely = 0.0;
    velz = 0.0;
    accx = 0.0;
    accy = 0.0;
    accz = 0.0;
  }

  OctTreeLeafNode(int type, double mass, double posx, double posy, double posz,
      double velx, double vely, double velz, double accx, double accy, double
      accz) : OctTreeNode(type, mass, posx, posy, posz), velx(velx), vely(vely),
  velz(velz), accx(accx), accy(accy), accz(accz) {} // expected-warning{{Inferred Effect Summary for OctTreeLeafNode: [reads(rpl([rLOCAL],[]))]:}}

  //OctTreeLeafNode(const OctTreeLeafNode& node) {
  //}

  void copyFrom PARAM(Rn) READS(Rn) WRITES(R) (OctTreeLeafNode *node ARG_(Rn)) { // expected-warning{{Inferred Effect Summary for copyFrom: [reads(rpl([p70_Rn],[])),reads(rpl([rLOCAL],[])),writes(rpl([p69_R],[]))]}}
    type = node->type;
    mass = node->mass;
    posx = node->posx;
    posy = node->posy;
    posz = node->posz;
    velx = node->velx;
    vely = node->vely;
    velz = node->velz;
    accx = node->accx;
    accy = node->accy;
    accz = node->accz;
  }

  void setVelocity WRITES(R) (const double x, const double y, const double z) { // expected-warning{{Inferred Effect Summary for setVelocity: [reads(rpl([rLOCAL],[])),writes(rpl([p69_R],[]))]:}}
    velx = x;
    vely = y;
    velz = z;
  }

  // advances a body's velocity and position by one time step
  void Advance READS(Global) WRITES(R) ();

  // computes the acceleration and velocity of a body
  void ComputeForce READS(*) WRITES(R) (const OctTreeInternalNode * const root, const double size);

private:
  // recursively walks the tree to compute the force on a body
  void RecurseForce READS(*) WRITES(R) (const OctTreeNode * const n, double dsq);

  double velx ARG_(R);
  double vely ARG_(R);
  double velz ARG_(R);
  double accx ARG_(R);
  double accy ARG_(R);
  double accz ARG_(R);
}; // end class OctTreeLeafNode

OctTreeInternalNode *OctTreeInternalNode::head = 0;
OctTreeInternalNode *OctTreeInternalNode::freelist = 0;

/*
void OctTreeInternalNode::NewNode PARAM(R) (OctTreeInternalNode * in ARG_(R), double px, double py, double pz) {
#ifdef DEBUG
  cout << "NewNode(px = " << px << ", py = " << py << ", pz = " << pz << ")" << endl;
#endif
  //OctTreeInternalNode *in ARG_(Local, R) = new OctTreeInternalNode;
  in->type = CELL;
  in->mass = 0.0;
  in->posx = px;
  in->posy = py;
  in->posz = pz;
  in->child0 = 0;
  in->child1 = 0;
  in->child2 = 0;
  in->child3 = 0;
  in->child4 = 0;
  in->child5 = 0;
  in->child6 = 0;
  in->child7 = 0;
}
*/

int OctTreeInternalNode::ChildID /*PARAM(Rb) READS(PiN, Rb)*/ (OctTreeLeafNode *const b ARG_(Rb)) { // expected-warning{{Inferred Effect Summary for ChildID: [reads(rpl([p50_PiN],[])),reads(rpl([p57_Rb],[])),writes(rpl([rLOCAL],[]))]}}
  int i = 0;

  if (posx < b->posx) {
    i = 1;
  }
  if (posy < b->posy) {
    i += 2;
  }
  if (posz < b->posz) {
    i += 4;
  }
#ifdef DEBUG
  cout << "ChildID(*this=" << *this << ", *b=" << *b << ") = " << i << endl;
#endif
  return i;
}

void OctTreeInternalNode::ChildIDToPos /*PARAM(Rx, Ry, Rz)*/ WRITES(Rx, Ry, Rz)
  (int childID, double radius, double *x ARG_(Rx), double *y ARG_(Ry), double *z ARG_(Rz)) { // expected-warning{{Inferred Effect Summary for ChildIDToPos: [writes(rpl([p54_Rx],[])),writes(rpl([p55_Ry],[])),writes(rpl([p56_Rz],[])),writes(rpl([rLOCAL],[]))]}}
  *x = *y = *z = 0;
  if (childID >= 4) {
    *z = radius;
    childID -= 4;
  }
  if (childID >= 2) {
    *y = radius;
    childID -= 2;
  }
  if (childID >= 1) {
    *x = radius;
    // childID -= 1;
  }
}

OctTreeNode** OctTreeInternalNode::GetChildRef ARG_(*, *) (int i) { // expected-warning{{Inferred Effect Summary for GetChildRef: [writes(rpl([rLOCAL],[]))]:}}
  OctTreeNode **childi ARG_(Local, *, *) = 0;
  switch (i) {
    case 0:
      childi = &child0;
      break;
    case 1:
      childi = &child1;
      break;
    case 2:
      childi = &child2;
      break;
    case 3:
      childi = &child3;
      break;
    case 4:
      childi = &child4;
      break;
    case 5:
      childi = &child5;
      break;
    case 6:
      childi = &child6;
      break;
    case 7:
      childi = &child7;
      break;
  }
  return childi;
}

OctTreeNode* OctTreeInternalNode::GetChild ARG_(*, *) READS(*) (int i) { // expected-warning{{Inferred Effect Summary for GetChild: [reads(rpl([rSTAR],[])),writes(rpl([rLOCAL],[]))]}}
  return *GetChildRef(i);
}

/*
void OctTreeInternalNode::Insert(OctTreeLeafNode * const b, const double r) // builds the tree
{
#ifdef DEBUG
  cout << "Insert(*this = {type = " << type << ", mass = " << mass << ", posx = " << posx << ", posy = " << posy << ", posz = " << posz << "}, " << "*b = {type = " << b->type << ", mass = " << b->mass << ", posx = " << b->posx << ", posy = " << b->posy << ", posz = " << b->posz << "}" << ", r = " << r << ")" << endl;
#endif
  int i = 0;
  double x = 0.0, y = 0.0, z = 0.0;

  if (posx < b->posx) {
    i = 1;
    x = r;
  }
  if (posy < b->posy) {
    i += 2;
    y = r;
  }
  if (posz < b->posz) {
    i += 4;
    z = r;
  }

  if (GetChild(i) == 0) {
    *GetChildRef(i) = b;
  } else if (GetChild(i)->type == CELL) {
    ((OctTreeInternalNode *) (GetChild(i)))->Insert(b, 0.5 * r);
  } else {
    const double rh = 0.5 * r;
    OctTreeInternalNode * const cell = NewNode(posx - rh + x, posy - rh + y, posz - rh + z);
    cell->Insert(b, rh);
    cell->Insert((OctTreeLeafNode *) (GetChild(i)), rh);
    *GetChildRef(i) = cell;
  }
}*/

OctTreeLeafNode **OctTreeInternalNode::GetPartition
ARG_(Rb:*, Rb:*)
/*PARAM(Rb)*/
(OctTreeLeafNode **b ARG_(Rb, Rb), int i,
 OctTreeLeafNode **partition0 ARG_(Rb:Rp0, Rb:Rp0),
 OctTreeLeafNode **partition1 ARG_(Rb:Rp1, Rb:Rp1),
 OctTreeLeafNode **partition2 ARG_(Rb:Rp2, Rb:Rp2),
 OctTreeLeafNode **partition3 ARG_(Rb:Rp3, Rb:Rp3),
 OctTreeLeafNode **partition4 ARG_(Rb:Rp4, Rb:Rp4),
 OctTreeLeafNode **partition5 ARG_(Rb:Rp5, Rb:Rp5),
 OctTreeLeafNode **partition6 ARG_(Rb:Rp6, Rb:Rp6),
 OctTreeLeafNode **partition7 ARG_(Rb:Rp7, Rb:Rp7)) { // expected-warning{{Inferred Effect Summary for GetPartition: [writes(rpl([rLOCAL],[]))]}}
  OctTreeLeafNode **partitioni ARG_(Local, Rb:*, Rb:*) = 0;
  switch (i) {
    case 0:
      partitioni = partition0;
      break;
    case 1:
      partitioni = partition1;
      break;
    case 2:
      partitioni = partition2;
      break;
    case 3:
      partitioni = partition3;
      break;
    case 4:
      partitioni = partition4;
      break;
    case 5:
      partitioni = partition5;
      break;
    case 6:
      partitioni = partition6;
      break;
    case 7:
      partitioni = partition7;
      break;
  }
  return partitioni;
}

void OctTreeInternalNode::Partition
/*PARAM(Rb, Rps) READS(PiN) WRITES(Rps, Rb:*) */
(OctTreeLeafNode **b ARG_(Rb, Rb), int n, int *partitionSize ARG_(Rps),
OctTreeLeafNode **partition0 ARG_(Rb:Rp0, Rb:Rp0),
OctTreeLeafNode **partition1 ARG_(Rb:Rp1, Rb:Rp1),
OctTreeLeafNode **partition2 ARG_(Rb:Rp2, Rb:Rp2),
OctTreeLeafNode **partition3 ARG_(Rb:Rp3, Rb:Rp3),
OctTreeLeafNode **partition4 ARG_(Rb:Rp4, Rb:Rp4),
OctTreeLeafNode **partition5 ARG_(Rb:Rp5, Rb:Rp5),
OctTreeLeafNode **partition6 ARG_(Rb:Rp6, Rb:Rp6),
OctTreeLeafNode **partition7 ARG_(Rb:Rp7, Rb:Rp7)) { // expected-warning{{Inferred Effect Summary for Partition: [reads(rpl([p50_PiN],[])),writes(rpl([p58_Rb,rSTAR],[])),writes(rpl([p59_Rps],[])),writes(rpl([rLOCAL],[]))]}}
  for (int i = 0; i < 8; ++i) {
    partitionSize[i] = 0;
  }
  for (int i = 0; i < n; ++i) {
    int partitionID = ChildID(b[i]);
    OctTreeLeafNode **partitioni ARG_(Local, Rb:*, Rb:*) = GetPartition(b, partitionID,
        partition0, partition1, partition2, partition3,
        partition4, partition5, partition6, partition7);
    //partitioni[partitionSize[partitionID]] = b[i];
    partitioni[partitionSize[partitionID]] = new OctTreeLeafNode();
    partitioni[partitionSize[partitionID]]->copyFrom(b[i]);
    ++partitionSize[partitionID];
  }
}

class PARAM(Rp, Rc) ChildrenInserter {
private:
  double r ARG_(Rc);
  double posx ARG_(Rc);
  double posy ARG_(Rc);
  double posz ARG_(Rc);
  int childID ARG_(Rc);
  OctTreeLeafNode **partition ARG_(Rp, Rp, Rp);
  int partitionSize ARG_(Rp);
  OctTreeNode **child ARG_(Rc, Rc, Rc);

public:
  ChildrenInserter(double r, double posx, double posy, double posz, int childID, OctTreeLeafNode **partition ARG_(Rp, Rp), int partitionSize, OctTreeNode **child ARG_(Rc, Rc))
                  : r(r), posx(posx), posy(posy), posz(posz), childID(childID), partition(partition), partitionSize(partitionSize), child(child) {}  // expected-warning{{Inferred Effect Summary for ChildrenInserter: [reads(rpl([rLOCAL],[]))]}}

  void operator() WRITES(Rp:*, Rc:*) () const { // expected-warning{{Inferred Effect Summary for operator(): [writes(rpl([p71_Rp,rSTAR],[])),writes(rpl([p72_Rc,rSTAR],[])),writes(rpl([rLOCAL],[]))]}}
    if (partitionSize > 1) {
      double x, y, z;
      OctTreeInternalNode::ChildIDToPos(childID, r, &x, &y, &z);
      double rh = 0.5 * r;
      OctTreeInternalNode *cell ARG_(Rc) = new OctTreeInternalNode(posx - rh + x, posy - rh + y, posz - rh + z);
      //OctTreeInternalNode *cell ARG_(Rc) = OctTreeInternalNode::NewNode(posx - rh + x, posy - rh + y, posz - rh + z);
      *child = cell;
      cell->InsertAll(partition, partitionSize, rh);
    } else if (partitionSize == 1) {
      //*child = partition[0];
      OctTreeLeafNode *leaf ARG_(Rc) = new OctTreeLeafNode();
      leaf->copyFrom(partition[0]);
      *child = leaf;
      /*
      *child = new OctTreeLeafNode();
      ((OctTreeLeafNode *) (*child))->copyFrom(partition[0]);*/
    }
  }
}; // end class ChildrenInserter

void OctTreeInternalNode::InsertChildren
//PARAM(Pp0, Pp1, Pp2, Pp3, Pp4, Pp5, Pp6, Pp7, Pps)
//READS(Pps)
//WRITES(PiN:*, Pp0:*, Pp1:*, Pp2:*, Pp3:*, Pp4:*, Pp5:*, Pp6:*, Pp7:*)
(double r, int *partitionSize ARG_(Pps),
OctTreeLeafNode **partition0 ARG_(Pp0, Pp0),
OctTreeLeafNode **partition1 ARG_(Pp1, Pp1),
OctTreeLeafNode **partition2 ARG_(Pp2, Pp2),
OctTreeLeafNode **partition3 ARG_(Pp3, Pp3),
OctTreeLeafNode **partition4 ARG_(Pp4, Pp4),
OctTreeLeafNode **partition5 ARG_(Pp5, Pp5),
OctTreeLeafNode **partition6 ARG_(Pp6, Pp6),
OctTreeLeafNode **partition7 ARG_(Pp7, Pp7)) { // expected-warning{{Inferred Effect Summary for InsertChildren: [reads(rpl([p50_PiN],[])),reads(rpl([p68_Pps],[])),writes(rpl([p50_PiN,r0_Rc0,rSTAR],[])),writes(rpl([p50_PiN,r1_Rc1,rSTAR],[])),writes(rpl([p50_PiN,r2_Rc2,rSTAR],[])),writes(rpl([p50_PiN,r3_Rc3,rSTAR],[])),writes(rpl([p50_PiN,r4_Rc4,rSTAR],[])),writes(rpl([p50_PiN,r5_Rc5,rSTAR],[])),writes(rpl([p50_PiN,r6_Rc6,rSTAR],[])),writes(rpl([p50_PiN,r7_Rc7,rSTAR],[])),writes(rpl([p60_Pp0,rSTAR],[])),writes(rpl([p61_Pp1,rSTAR],[])),writes(rpl([p62_Pp2,rSTAR],[])),writes(rpl([p63_Pp3,rSTAR],[])),writes(rpl([p64_Pp4,rSTAR],[])),writes(rpl([p65_Pp5,rSTAR],[])),writes(rpl([p66_Pp6,rSTAR],[])),writes(rpl([p67_Pp7,rSTAR],[])),writes(rpl([rLOCAL],[]))]}}
  ChildrenInserter inserter0 ARG_(Pp0, PiN:Rc0) (r, posx, posy, posz, 0, partition0, partitionSize[0], &child0);
  ChildrenInserter inserter1 ARG_(Pp1, PiN:Rc1) (r, posx, posy, posz, 1, partition1, partitionSize[1], &child1);
  ChildrenInserter inserter2 ARG_(Pp2, PiN:Rc2) (r, posx, posy, posz, 2, partition2, partitionSize[2], &child2);
  ChildrenInserter inserter3 ARG_(Pp3, PiN:Rc3) (r, posx, posy, posz, 3, partition3, partitionSize[3], &child3);
  ChildrenInserter inserter4 ARG_(Pp4, PiN:Rc4) (r, posx, posy, posz, 4, partition4, partitionSize[4], &child4);
  ChildrenInserter inserter5 ARG_(Pp5, PiN:Rc5) (r, posx, posy, posz, 5, partition5, partitionSize[5], &child5);
  ChildrenInserter inserter6 ARG_(Pp6, PiN:Rc6) (r, posx, posy, posz, 6, partition6, partitionSize[6], &child6);
  ChildrenInserter inserter7 ARG_(Pp7, PiN:Rc7) (r, posx, posy, posz, 7, partition7, partitionSize[7], &child7);
  parallel_invoke(inserter0, inserter1, inserter2, inserter3, inserter4, inserter5, inserter6, inserter7);
}

void OctTreeInternalNode::InsertAll /*PARAM(Rb) WRITES(Rb:*, PiN:*) */(OctTreeLeafNode **b ARG_(Rb, Rb), int n, double r) { // expected-warning{{Inferred Effect Summary for InsertAll: [reads(rpl([p50_PiN],[])),writes(rpl([p50_PiN,r0_Rc0,rSTAR],[])),writes(rpl([p50_PiN,r1_Rc1,rSTAR],[])),writes(rpl([p50_PiN,r2_Rc2,rSTAR],[])),writes(rpl([p50_PiN,r3_Rc3,rSTAR],[])),writes(rpl([p50_PiN,r4_Rc4,rSTAR],[])),writes(rpl([p50_PiN,r5_Rc5,rSTAR],[])),writes(rpl([p50_PiN,r6_Rc6,rSTAR],[])),writes(rpl([p50_PiN,r7_Rc7,rSTAR],[])),writes(rpl([p52_Rb,rSTAR],[])),writes(rpl([rLOCAL],[]))]}}
#ifdef DEBUG
  cout << "InsertAll(*this = " << *this << ", n = " << n << ", r = " << r << ")" << endl;
#endif
  if (n == 0) {
    return;
  }
  //OctTreeLeafNode **partition = new OctTreeLeafNode**[8];
  OctTreeLeafNode **partition0 ARG_(Local, Rb:Rp0, Rb:Rp0);
  OctTreeLeafNode **partition1 ARG_(Local, Rb:Rp1, Rb:Rp1);
  OctTreeLeafNode **partition2 ARG_(Local, Rb:Rp2, Rb:Rp2);
  OctTreeLeafNode **partition3 ARG_(Local, Rb:Rp3, Rb:Rp3);
  OctTreeLeafNode **partition4 ARG_(Local, Rb:Rp4, Rb:Rp4);
  OctTreeLeafNode **partition5 ARG_(Local, Rb:Rp5, Rb:Rp5);
  OctTreeLeafNode **partition6 ARG_(Local, Rb:Rp6, Rb:Rp6);
  OctTreeLeafNode **partition7 ARG_(Local, Rb:Rp7, Rb:Rp7);
  int partitionSize ARG_(Local) [8];
  partition0 = new OctTreeLeafNode*[n];
  partition1 = new OctTreeLeafNode*[n];
  partition2 = new OctTreeLeafNode*[n];
  partition3 = new OctTreeLeafNode*[n];
  partition4 = new OctTreeLeafNode*[n];
  partition5 = new OctTreeLeafNode*[n];
  partition6 = new OctTreeLeafNode*[n];
  partition7 = new OctTreeLeafNode*[n];
  Partition(b, n, partitionSize, partition0, partition1, partition2,
      partition3, partition4, partition5, partition6, partition7);
  InsertChildren(r, partitionSize, partition0, partition1, partition2,
      partition3, partition4, partition5, partition6, partition7);
  //for (int i = 0; i < 8; ++i) {
  //  delete partition[i];
  //}
  //delete partition;
  //for (int i = 0; i < 8; ++i) {
  //  delete [] GetPartition(b, i, partition0, partition1, partition2,
  //      partition3, partition4, partition5, partition6,
  //      partition7);
  //}
}

// Recursively summarizes info about subtrees
void OctTreeInternalNode::ComputeCenterOfMass WRITES(*) (int &curr ARG_(Global)) { // expected-warning{{Inferred Effect Summary for ComputeCenterOfMass: [writes(rpl([rSTAR],[]))]}}
  double m, px = 0.0, py = 0.0, pz = 0.0;
  OctTreeNode *ch ARG_(*, *);

  int j = 0;
  mass = 0.0;
  for (int i = 0; i < 8; i++) {
    //ch = child[i];
    ch = *GetChildRef(i);
    if (ch != 0) {
      //child[i] = 0;
      //child[j++] = ch;
      // move non-NULL children to the front (needed to make other code faster)
      *GetChildRef(i) = 0;
      *GetChildRef(j++) = ch;

      if (ch->type == BODY) {
        // sort bodies in tree order (approximation of putting nearby nodes
        // together for locality)
        bodies[curr++] = (OctTreeLeafNode *) ch;
      } else {
        ((OctTreeInternalNode *) ch)->ComputeCenterOfMass(curr);
      }
      m = ch->mass;
      mass += m;
      px += ch->posx * m;
      py += ch->posy * m;
      pz += ch->posz * m;
    }
  }

  m = 1.0 / mass;
  posx = px * m;
  posy = py * m;
  posz = pz * m;
}

// Advances a body's velocity and position by one time step
void OctTreeLeafNode::Advance READS(Global) WRITES(R) () { // expected-warning{{Inferred Effect Summary for Advance: [reads(rpl([rGLOBAL],[])),writes(rpl([p69_R],[])),writes(rpl([rLOCAL],[]))]}}
  double dvelx, dvely, dvelz;
  double velhx, velhy, velhz;

  dvelx = accx * dthf;
  dvely = accy * dthf;
  dvelz = accz * dthf;

  velhx = velx + dvelx;
  velhy = vely + dvely;
  velhz = velz + dvelz;

  posx += velhx * dtime;
  posy += velhy * dtime;
  posz += velhz * dtime;

  velx = velhx + dvelx;
  vely = velhy + dvely;
  velz = velhz + dvelz;
}

// Computes the acceleration and velocity of a body
void OctTreeLeafNode::ComputeForce READS(*) WRITES(R) (const OctTreeInternalNode * const root, const double size) { // expected-warning{{Inferred Effect Summary for ComputeForce: [reads(rpl([rSTAR],[])),writes(rpl([p69_R],[])),writes(rpl([rLOCAL],[]))]}}
  double ax, ay, az;

  ax = accx;
  ay = accy;
  az = accz;

  accx = 0.0;
  accy = 0.0;
  accz = 0.0;

  RecurseForce(root, size * size * itolsq);

  if (step > 0) {
    velx += (accx - ax) * dthf;
    vely += (accy - ay) * dthf;
    velz += (accz - az) * dthf;
  }
}

// Recursively walks the tree to compute the force on a body
void OctTreeLeafNode::RecurseForce READS(*) WRITES(R) (const OctTreeNode * const n ARG_(Global), double dsq) { // expected-warning{{Inferred Effect Summary for RecurseForce: [reads(rpl([rSTAR],[])),writes(rpl([p69_R],[])),writes(rpl([rLOCAL],[]))]}}
  double drx, dry, drz, drsq, nphi, scale, idr;

  drx = n->posx - posx;
  dry = n->posy - posy;
  drz = n->posz - posz;
  drsq = drx * drx + dry * dry + drz * drz;
  if (drsq < dsq) {
    if (n->type == CELL) {
      OctTreeInternalNode *in ARG_(Global) = (OctTreeInternalNode *) n;
      dsq *= 0.25;
      if (in->GetChild(0) != 0) {
        RecurseForce(in->GetChild(0), dsq);
        if (in->GetChild(1) != 0) {
          RecurseForce(in->GetChild(1), dsq);
          if (in->GetChild(2) != 0) {
            RecurseForce(in->GetChild(2), dsq);
            if (in->GetChild(3) != 0) {
              RecurseForce(in->GetChild(3), dsq);
              if (in->GetChild(4) != 0) {
                RecurseForce(in->GetChild(4), dsq);
                if (in->GetChild(5) != 0) {
                  RecurseForce(in->GetChild(5), dsq);
                  if (in->GetChild(6) != 0) {
                    RecurseForce(in->GetChild(6), dsq);
                    if (in->GetChild(7) != 0) {
                      RecurseForce(in->GetChild(7), dsq);
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else { // n is a body
      if (n != this) {
        drsq += epssq;
        idr = 1 / sqrt(drsq);
        nphi = n->mass * idr;
        scale = nphi * idr * idr;
        accx += drx * scale;
        accy += dry * scale;
        accz += drz * scale;
      }
    }
  } else { // node is far enough away, don't recurse any deeper
    drsq += epssq;
    idr = 1 / sqrt(drsq);
    nphi = n->mass * idr;
    scale = nphi * idr * idr;
    accx += drx * scale;
    accy += dry * scale;
    accz += drz * scale;
  }
}

static int nbodies; // number of bodies in system
static int timesteps; // number of time steps to run
static int grainSize; // number of parallel tasks

static inline int ReadInput WRITES(Global) (char *filename) { // expected-warning{{Inferred Effect Summary for ReadInput: [writes(rpl([rGLOBAL],[])),writes(rpl([rLOCAL],[]))]}}
  double vx, vy, vz;
  FILE *f ARG_(Local, *);

  f = fopen(filename, "r+t");
  if (f == 0) {
    fprintf(stderr, "file not found: %s\n", filename);
    //exit(-1);
    return 1;
  }

  fscanf(f, "%d", &nbodies);
  fscanf(f, "%d", &timesteps);
  fscanf(f, "%lf", &dtime);
  fscanf(f, "%lf", &eps);
  fscanf(f, "%lf", &tol);

  dthf = 0.5 * dtime;
  epssq = eps * eps;
  itolsq = 1.0 / (tol * tol);

  if (bodies == 0) {
    if (grainSize != 0)
      fprintf(stderr,
      "configuration: %d bodies, %d time steps, %d grain size (manually specified)\n", nbodies, timesteps,
          grainSize);
    else
      fprintf(stderr,
      "configuration: %d bodies, %d time steps and automatically determined grain size\n", nbodies, timesteps);

    bodies = new OctTreeLeafNode*[nbodies];
    for (int i = 0; i < nbodies; i++)
      bodies[i] = new OctTreeLeafNode();
  }

  for (int i = 0; i < nbodies; i++) {
    fscanf(f, "%lE", &(bodies[i]->mass));
    fscanf(f, "%lE", &(bodies[i]->posx));
    fscanf(f, "%lE", &(bodies[i]->posy));
    fscanf(f, "%lE", &(bodies[i]->posz));
    fscanf(f, "%lE", &vx);
    fscanf(f, "%lE", &vy);
    fscanf(f, "%lE", &vz);
    bodies[i]->setVelocity(vx, vy, vz);
  }

  fclose(f);
  return 0;
}

static inline void ComputeCenterAndDiameter WRITES(Global)
  (const int n, double &diameter ARG_(Global), double &centerx ARG_(Global),
   double &centery ARG_(Global), double &centerz ARG_(Global)) { // expected-warning{{Inferred Effect Summary for ComputeCenterAndDiameter: [writes(rpl([rGLOBAL],[])),writes(rpl([rLOCAL],[]))]}}
  double minx, miny, minz;
  double maxx, maxy, maxz;
  double posx, posy, posz;

  minx = 1.0E90;
  miny = 1.0E90;
  minz = 1.0E90;
  maxx = -1.0E90;
  maxy = -1.0E90;
  maxz = -1.0E90;

  for (int i = 0; i < n; i++) {
    posx = bodies[i]->posx;
    posy = bodies[i]->posy;
    posz = bodies[i]->posz;

    if (minx > posx)
      minx = posx;
    if (miny > posy)
      miny = posy;
    if (minz > posz)
      minz = posz;

    if (maxx < posx)
      maxx = posx;
    if (maxy < posy)
      maxy = posy;
    if (maxz < posz)
      maxz = posz;
  }

  diameter = maxx - minx;
  if (diameter < (maxy - miny))
    diameter = (maxy - miny);
  if (diameter < (maxz - minz))
    diameter = (maxz - minz);

  centerx = (maxx + minx) * 0.5;
  centery = (maxy + miny) * 0.5;
  centerz = (maxz + minz) * 0.5;
}

static inline int min(long a, long b) { // expected-warning{{Inferred Effect Summary for min: [writes(rpl([rLOCAL],[]))]}}
  if (a < b)
    a = b;
  return a;
}

static void PrintDouble(double d) { // expected-warning{{Inferred Effect Summary for PrintDouble: [writes(rpl([rLOCAL],[]))]}}
  int i;
  char str ARG_(Local) [16];

  sprintf(str, "%.4lE", d);

  i = 0;
  while ((i < 16) && (str[i] != 0)) {
    if ((str[i] == 'E') && (str[i + 1] == '-') && (str[i + 2] == '0') && (str[i + 3] == '0')) {
      printf("E00");
      i += 3;
    } else if (str[i] != '+') {
      printf("%c", str[i]);
    }
    i++;
  }
}

static OctTreeInternalNode *root ARG_(Global);
static double gDiameter;

class ParallelForProcessor {
public:
  void operator() READS(*) WRITES(Global) (const blocked_range<int>& range) const { // expected-warning{{Inferred Effect Summary for operator(): [reads(rpl([rSTAR],[])),writes(rpl([rGLOBAL],[])),writes(rpl([rLOCAL],[]))]}}
    for (int i = range.begin(); i != range.end(); i++) {
      bodies[i]->ComputeForce(root, gDiameter);
    }
  }
};

int main WRITES(*) (int argc, char **argv ARG_(Global, Global)) { // expected-warning{{Inferred Effect Summary for main: [writes(rpl([rGLOBAL,r0_Rc0,rSTAR],[])),writes(rpl([rGLOBAL,r1_Rc1,rSTAR],[])),writes(rpl([rGLOBAL,r2_Rc2,rSTAR],[])),writes(rpl([rGLOBAL,r3_Rc3,rSTAR],[])),writes(rpl([rGLOBAL,r4_Rc4,rSTAR],[])),writes(rpl([rGLOBAL,r5_Rc5,rSTAR],[])),writes(rpl([rGLOBAL,r6_Rc6,rSTAR],[])),writes(rpl([rGLOBAL,r7_Rc7,rSTAR],[])),writes(rpl([rLOCAL],[])),writes(rpl([rSTAR],[]))]}}
  //task_scheduler_init init;
  ParallelForProcessor parallelProcessor;

  fprintf(stderr, "\n");
  fprintf(stderr, "Lonestar benchmark suite\n");
  fprintf(stderr,
  "Copyright (C) 2007, 2008 The University of Texas at Austin\n");
  fprintf(stderr, "http://iss.ices.utexas.edu/lonestar/\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "application: BarnesHut v1.0\n");

  bodies = 0;

  if (argc == 3) {
    fprintf(stderr, "The support for grain size is turned off.\n");
    //grainSize = atoi(argv[2]);
  }

  timeval starttime, endtime;
  long runtime, lasttime, mintime;
  int run;

  runtime = 0;
  lasttime = -1;
  mintime = -1;
  run = 0;

  while (((run < 3) || (abs(lasttime - runtime) * 64 > min(lasttime, runtime))) && (run < 7)) {
    if (ReadInput(argv[1]) != 0) {
      return 1;
    }

    lasttime = runtime;
    gettimeofday(&starttime, 0);

    for (step = 0; step < timesteps; step++) { // time-step the system
      double diameter ARG_(Global), centerx ARG_(Global), centery ARG_(Global), centerz ARG_(Global);
      ComputeCenterAndDiameter(nbodies, diameter, centerx, centery, centerz);

      // create the tree's root
      OctTreeInternalNode *local_root ARG_(Global) = new OctTreeInternalNode(centerx, centery, centerz); 
      //OctTreeInternalNode *local_root ARG_(Local) = OctTreeInternalNode::NewNode(centerx, centery, centerz);

      const double radius = diameter * 0.5;
//      for (int i = 0; i < nbodies; ++i) {
//        local_root->Insert(bodies[i], radius); // grow the tree by inserting each body
//      }
      local_root->InsertAll(bodies, nbodies, radius);

      int curr ARG_(Global) = 0;
      // summarize subtree info in each internal node (plus restructure tree and
      // sort bodies for performance reasons)
      local_root->ComputeCenterOfMass(curr);

      root = local_root;
      gDiameter = diameter;

      /* Commenting out because we need index parameterized arrays to support this
      if (grainSize != 0) {
        blocked_range<int> grainedRange(0, nbodies, grainSize);
        parallel_for(grainedRange, parallelProcessor);
      } else {
        blocked_range<int> range(0, nbodies);
        parallel_for(range, parallelProcessor);
      }*/

      //OctTreeInternalNode::RecycleTree(); // recycle the tree

      // the iterations are independent: they can be executed in any order and in parallel
      for (int i = 0; i < nbodies; i++) {
        // advance the position and velocity of each body
        bodies[i]->Advance();
      }
    } // end of time step

    gettimeofday(&endtime, 0);
    runtime = (long) (endtime.tv_sec * 1000.0 + endtime.tv_usec / 1000.0 - starttime.tv_sec * 1000.0
        - starttime.tv_usec / 1000.0 + 0.5);

    if ((runtime < mintime) || (run == 0))
      mintime = runtime;
    run++;
  }

  fprintf(stderr, "runtime: %ld ms\n\n", mintime);

  for (int i = 0; i < nbodies; i++) { // print result
    PrintDouble(bodies[i]->posx);
    printf(" ");
    PrintDouble(bodies[i]->posy);
    printf(" ");
    PrintDouble(bodies[i]->posz);
    printf("\n");
  }

  return 0;
}

