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

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sys/time.h>
#include <iostream>

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_invoke.h"
//#include "../include/blocked_range.h"
//#include "../include/parallel_for.h"
//#include "../include/parallel_invoke.h"
//#include "tbb/task_scheduler_init.h"
//#include "tbb/task_group.h"

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

class OctTreeLeafNode;

static OctTreeLeafNode **bodies; // the n bodies


class OctTreeNode {
public:
	int type; // CELL or BODY
	double mass;
	double posx;
	double posy;
	double posz;
	friend std::ostream& operator<<(std::ostream &stream, const OctTreeNode &node);
};

std::ostream &operator<<(std::ostream &stream, const OctTreeNode &node) {
	stream << "{posx=" << node.posx << ", posy=" << node.posy << ", posz="
		<< node.posz << "}" << flush;
	return stream;
}

class OctTreeInternalNode: public OctTreeNode { // the internal nodes are cells that summarize their children's properties
public:
	static OctTreeInternalNode *NewNode(const double px, const double py, const double pz);

	static void RecycleTree() {
		freelist = head;
	}

public:
	void Insert(OctTreeLeafNode * const b, const double r); // builds the tree
	void InsertAll(OctTreeLeafNode ** const b, const int n, const double r);
	void ComputeCenterOfMass(int &curr); // recursively summarizes info about subtrees
	OctTreeNode** GetChildRef(int i);
	OctTreeNode* GetChild(int i);
	OctTreeLeafNode ***GetPartitionRef(int i, OctTreeLeafNode **partition0,
			OctTreeLeafNode **partition1, OctTreeLeafNode **partition2,
			OctTreeLeafNode **partition3, OctTreeLeafNode **partition4,
			OctTreeLeafNode **partition5, OctTreeLeafNode **partition6,
			OctTreeLeafNode **partition7);
	OctTreeLeafNode **GetPartition(int i, OctTreeLeafNode **partition0,
			OctTreeLeafNode **partition1, OctTreeLeafNode **partition2,
			OctTreeLeafNode **partition3, OctTreeLeafNode **partition4,
			OctTreeLeafNode **partition5, OctTreeLeafNode **partition6,
			OctTreeLeafNode **partition7);
	static void ChildIDToPos(int childID, double radius, double &x, double &y, double &z);
	//OctTreeNode *child[8];
	OctTreeNode *child0, *child1, *child2, *child3, *child4, *child5,
		    *child6, *child7;

private:
	//OctTreeInternalNode *link; // links all internal tree nodes so they can be recycled
	static OctTreeInternalNode *head, *freelist; // free list for recycling

	int ChildID(OctTreeLeafNode * const b);
	void Partition(OctTreeLeafNode **b, const int n, int *partitionSize,
			OctTreeLeafNode **partition0, OctTreeLeafNode **partition1,
			OctTreeLeafNode **partition2, OctTreeLeafNode **partition3,
			OctTreeLeafNode **partition4, OctTreeLeafNode **partition5,
			OctTreeLeafNode **partition6, OctTreeLeafNode **partition7); 
	void InsertChildren(double r, int *partitionSize, OctTreeLeafNode
			**partition0, OctTreeLeafNode **partition1, OctTreeLeafNode
			**partition2, OctTreeLeafNode **partition3, OctTreeLeafNode
			**partition4, OctTreeLeafNode **partition5, OctTreeLeafNode
			**partition6, OctTreeLeafNode **partition7); 
};

class OctTreeLeafNode: public OctTreeNode { // the tree leaves are the bodies
public:
	OctTreeLeafNode();
	~OctTreeLeafNode() {
	}

	void setVelocity(const double x, const double y, const double z) {
		velx = x;
		vely = y;
		velz = z;
	}

	// advances a body's velocity and position by one time step
	void Advance();

	// computes the acceleration and velocity of a body
	void ComputeForce(const OctTreeInternalNode * const root, const double size);

private:
	void RecurseForce(const OctTreeNode * const n, double dsq); // recursively walks the tree to compute the force on a body

	double velx;
	double vely;
	double velz;
	double accx;
	double accy;
	double accz;
};

OctTreeInternalNode *OctTreeInternalNode::head = NULL;
OctTreeInternalNode *OctTreeInternalNode::freelist = NULL;

OctTreeInternalNode *OctTreeInternalNode::NewNode(const double px, const double py, const double pz) {
#ifdef DEBUG
	cout << "NewNode(px = " << px << ", py = " << py << ", pz = " << pz << ")" << endl;
#endif
	OctTreeInternalNode *in;

	in = new OctTreeInternalNode();
	//Accesses to freelist should be atomic when NewNode is called
	//concurrently.
	//if (freelist == NULL) {
	//	in = new OctTreeInternalNode();
	//	in->link = head;
	//	head = in;
	//} else { // get node from freelist
	//	in = freelist;
	//	freelist = freelist->link;
	//}

	in->type = CELL;
	in->mass = 0.0;
	in->posx = px;
	in->posy = py;
	in->posz = pz;
	//for (int i = 0; i < 8; i++)
	//	in->child[i] = NULL;
	in->child0 = NULL;
	in->child1 = NULL;
	in->child2 = NULL;
	in->child3 = NULL;
	in->child4 = NULL;
	in->child5 = NULL;
	in->child6 = NULL;
	in->child7 = NULL;

	return in;
}

int OctTreeInternalNode::ChildID(OctTreeLeafNode * const b)
{
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

void OctTreeInternalNode::ChildIDToPos(int childID, double radius, double &x, double &y, double &z) {
	x = y = z = 0;
	if (childID >= 4) {
		z = radius;
		childID -= 4;
	}
	if (childID >= 2) {
		y = radius;
		childID -= 2;
	}
	if (childID >= 1) {
		x = radius;
		childID -= 1;
	}
}

OctTreeNode** OctTreeInternalNode::GetChildRef(int i) {
	OctTreeNode** childi = 0;
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

OctTreeNode* OctTreeInternalNode::GetChild(int i) {
	return *GetChildRef(i);
}

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

	if (GetChild(i) == NULL) {
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
}

//void OctTreeInternalNode::InsertAll(OctTreeLeafNode ** const b, const int n, const double r)
//{
//#ifdef DEBUG
//	cout << "InsertAll(*this = " << *this << ", n = " << n << ", r = " << r << ")" << endl;
//#endif
//	if (n == 0) {
//		return;
//	}
//	OctTreeLeafNode ***partition = new OctTreeLeafNode**[8];
//	int partitionSize[8];
//	for (int i = 0; i < 8; ++i) {
//		partition[i] = new OctTreeLeafNode*[n];
//		partitionSize[i] = 0;
//	}
//	for (int i = 0; i < n; ++i) {
//		int partitionID = ChildID(b[i]);
//		partition[partitionID][partitionSize[partitionID]] = b[i];
//		++partitionSize[partitionID];
//	}
//	// https://software.intel.com/en-us/node/506118
//	tbb::task_group g;
//	for (int i = 0; i < 8; ++i) {
//		if (partitionSize[i] > 1) {
//			double x, y, z;
//			ChildIDToPos(i, r, x, y, z);
//			const double rh = 0.5 * r;
//			OctTreeInternalNode * const cell = NewNode(posx - rh + x, posy - rh + y, posz - rh + z);
//			child[i] = cell;
//			//cell->InsertAll(partition[i], partitionSize[i], rh);
//			g.run([=]{cell->InsertAll(partition[i], partitionSize[i], rh);});
//		} else if (partitionSize[i] == 1) {
//			child[i] = partition[i][0];
//		}
//	}
//	g.wait();
//	for (int i = 0; i < 8; ++i) {
//		delete partition[i];
//	}
//	delete partition;
//}

OctTreeLeafNode ***OctTreeInternalNode::GetPartitionRef(int i, OctTreeLeafNode
		**partition0, OctTreeLeafNode **partition1, OctTreeLeafNode
		**partition2, OctTreeLeafNode **partition3, OctTreeLeafNode
		**partition4, OctTreeLeafNode **partition5, OctTreeLeafNode
		**partition6, OctTreeLeafNode **partition7) {
	OctTreeLeafNode ***partitioni = 0;
	switch (i) {
		case 0:
			partitioni = &partition0;
			break;
		case 1:
			partitioni = &partition1;
			break;
		case 2:
			partitioni = &partition2;
			break;
		case 3:
			partitioni = &partition3;
			break;
		case 4:
			partitioni = &partition4;
			break;
		case 5:
			partitioni = &partition5;
			break;
		case 6:
			partitioni = &partition6;
			break;
		case 7:
			partitioni = &partition7;
			break;
	}
	return partitioni;
}

OctTreeLeafNode **OctTreeInternalNode::GetPartition(int i, OctTreeLeafNode **partition0,
		OctTreeLeafNode **partition1, OctTreeLeafNode **partition2,
		OctTreeLeafNode **partition3, OctTreeLeafNode **partition4,
		OctTreeLeafNode **partition5, OctTreeLeafNode **partition6,
		OctTreeLeafNode **partition7) {
	return *GetPartitionRef(i, partition0, partition1, partition2,
			partition3, partition4, partition5, partition6,
			partition7);
}

void OctTreeInternalNode::Partition(OctTreeLeafNode **b, const int n, int
		*partitionSize, OctTreeLeafNode **partition0, OctTreeLeafNode
		**partition1, OctTreeLeafNode **partition2, OctTreeLeafNode
		**partition3, OctTreeLeafNode **partition4, OctTreeLeafNode
		**partition5, OctTreeLeafNode **partition6, OctTreeLeafNode
		**partition7) {
	for (int i = 0; i < 8; ++i) {
		//partition[i] = new OctTreeLeafNode*[n];
		partitionSize[i] = 0;
	}
	for (int i = 0; i < n; ++i) {
		int partitionID = ChildID(b[i]);
		//partition[partitionID][partitionSize[partitionID]] = b[i];
		OctTreeLeafNode **partitioni = GetPartition(partitionID,
				partition0, partition1, partition2, partition3,
				partition4, partition5, partition6, partition7);
		partitioni[partitionSize[partitionID]] = b[i];
		++partitionSize[partitionID];
	}
}

class ChildrenInserter {
private:
	double r, posx, posy, posz;
	int childID;
	OctTreeLeafNode **partition;
	int partitionSize;
	OctTreeNode **child;

public:
	ChildrenInserter(double r, double posx, double posy, double posz, int childID, OctTreeLeafNode **partition, int partitionSize, OctTreeNode **child): r(r), posx(posx), posy(posy), posz(posz), childID(childID), partition(partition), partitionSize(partitionSize), child(child) {}

	void operator()() const {
		if (partitionSize > 1) {
			double x, y, z;
			OctTreeInternalNode::ChildIDToPos(childID, r, x, y, z);
			double rh = 0.5 * r;
			OctTreeInternalNode *cell = OctTreeInternalNode::NewNode(posx - rh + x, posy - rh + y, posz - rh + z);
			*child = cell;
			cell->InsertAll(partition, partitionSize, rh);
		} else if (partitionSize == 1) {
			*child = partition[0];
		}
	}
};

void OctTreeInternalNode::InsertChildren(double r, int *partitionSize, OctTreeLeafNode
		**partition0, OctTreeLeafNode **partition1, OctTreeLeafNode
		**partition2, OctTreeLeafNode **partition3, OctTreeLeafNode
		**partition4, OctTreeLeafNode **partition5, OctTreeLeafNode
		**partition6, OctTreeLeafNode **partition7) {
	ChildrenInserter inserter0 = ChildrenInserter(r, posx, posy, posz, 0, partition0, partitionSize[0], GetChildRef(0));
	ChildrenInserter inserter1 = ChildrenInserter(r, posx, posy, posz, 1, partition1, partitionSize[1], GetChildRef(1));
	ChildrenInserter inserter2 = ChildrenInserter(r, posx, posy, posz, 2, partition2, partitionSize[2], GetChildRef(2));
	ChildrenInserter inserter3 = ChildrenInserter(r, posx, posy, posz, 3, partition3, partitionSize[3], GetChildRef(3));
	ChildrenInserter inserter4 = ChildrenInserter(r, posx, posy, posz, 4, partition4, partitionSize[4], GetChildRef(4));
	ChildrenInserter inserter5 = ChildrenInserter(r, posx, posy, posz, 5, partition5, partitionSize[5], GetChildRef(5));
	ChildrenInserter inserter6 = ChildrenInserter(r, posx, posy, posz, 6, partition6, partitionSize[6], GetChildRef(6));
	ChildrenInserter inserter7 = ChildrenInserter(r, posx, posy, posz, 7, partition7, partitionSize[7], GetChildRef(7));
	parallel_invoke(inserter0, inserter1, inserter2, inserter3, inserter4, inserter5, inserter6, inserter7);
}

//void OctTreeInternalNode::InsertChildren(double r, int *partitionSize, OctTreeLeafNode
//		**partition0, OctTreeLeafNode **partition1, OctTreeLeafNode
//		**partition2, OctTreeLeafNode **partition3, OctTreeLeafNode
//		**partition4, OctTreeLeafNode **partition5, OctTreeLeafNode
//		**partition6, OctTreeLeafNode **partition7) {
//	// https://software.intel.com/en-us/node/506118
//	tbb::task_group g;
//	for (int i = 0; i < 8; ++i) {
//		OctTreeLeafNode **partitioni = GetPartition(i, partition0,
//				partition1, partition2, partition3, partition4,
//				partition5, partition6, partition7);
//		OctTreeNode **childi = GetChildRef(i);
//		if (partitionSize[i] > 1) {
//			double x, y, z;
//			ChildIDToPos(i, r, x, y, z);
//			const double rh = 0.5 * r;
//			OctTreeInternalNode * const cell = NewNode(posx - rh + x, posy - rh + y, posz - rh + z);
//			*childi = cell;
//			//child[i] = cell;
//			//cell->InsertAll(partition[i], partitionSize[i], rh);
//			//g.run([=]{cell->InsertAll(partition[i], partitionSize[i], rh);});
//			g.run([=]{cell->InsertAll(partitioni, partitionSize[i], rh);});
//		} else if (partitionSize[i] == 1) {
//			//child[i] = partition[i][0];
//			*childi = partitioni[0];
//		}
//	}
//	g.wait();
//}

void OctTreeInternalNode::InsertAll(OctTreeLeafNode ** const b, const int n, const double r)
{
#ifdef DEBUG
	cout << "InsertAll(*this = " << *this << ", n = " << n << ", r = " << r << ")" << endl;
#endif
	if (n == 0) {
		return;
	}
	//OctTreeLeafNode **partition = new OctTreeLeafNode**[8];
	OctTreeLeafNode **partition0, **partition1, **partition2, **partition3,
			**partition4, **partition5, **partition6, **partition7;
	int partitionSize[8];
	partition0 = new OctTreeLeafNode*[n];
	partition1 = new OctTreeLeafNode*[n];
	partition2 = new OctTreeLeafNode*[n];
	partition3 = new OctTreeLeafNode*[n];
	partition4 = new OctTreeLeafNode*[n];
	partition5 = new OctTreeLeafNode*[n];
	partition6 = new OctTreeLeafNode*[n];
	partition7 = new OctTreeLeafNode*[n];
	Partition(b, n, partitionSize, partition0, partition1, partition2,
			partition3, partition4, partition5, partition6,
			partition7);
	InsertChildren(r, partitionSize, partition0, partition1, partition2,
			partition3, partition4, partition5, partition6,
			partition7);
	//for (int i = 0; i < 8; ++i) {
	//	delete partition[i];
	//}
	//delete partition;
	for (int i = 0; i < 8; ++i) {
		delete GetPartition(i, partition0, partition1, partition2,
				partition3, partition4, partition5, partition6,
				partition7);
	}
}

void OctTreeInternalNode::ComputeCenterOfMass(int &curr) // recursively summarizes info about subtrees
{
	double m, px = 0.0, py = 0.0, pz = 0.0;
	OctTreeNode *ch;

	int j = 0;
	mass = 0.0;
	for (int i = 0; i < 8; i++) {
		//ch = child[i];
		ch = *GetChildRef(i);
		if (ch != NULL) {
			//child[i] = NULL; // move non-NULL children to the front (needed to make other code faster)
			//child[j++] = ch;
			*GetChildRef(i) = NULL; // move non-NULL children to the front (needed to make other code faster)
			*GetChildRef(j++) = ch;

			if (ch->type == BODY) {
				bodies[curr++] = (OctTreeLeafNode *) ch; // sort bodies in tree order (approximation of putting nearby nodes together for locality)
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

OctTreeLeafNode::OctTreeLeafNode() {
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

void OctTreeLeafNode::Advance() // advances a body's velocity and position by one time step
{
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

void OctTreeLeafNode::ComputeForce(const OctTreeInternalNode * const root, const double size) // computes the acceleration and velocity of a body
{
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

void OctTreeLeafNode::RecurseForce(const OctTreeNode * const n, double dsq) // recursively walks the tree to compute the force on a body
{
	double drx, dry, drz, drsq, nphi, scale, idr;

	drx = n->posx - posx;
	dry = n->posy - posy;
	drz = n->posz - posz;
	drsq = drx * drx + dry * dry + drz * drz;
	if (drsq < dsq) {
		if (n->type == CELL) {
			OctTreeInternalNode *in = (OctTreeInternalNode *) n;
			dsq *= 0.25;
			if (in->GetChild(0) != NULL) {
				RecurseForce(in->GetChild(0), dsq);
				if (in->GetChild(1) != NULL) {
					RecurseForce(in->GetChild(1), dsq);
					if (in->GetChild(2) != NULL) {
						RecurseForce(in->GetChild(2), dsq);
						if (in->GetChild(3) != NULL) {
							RecurseForce(in->GetChild(3), dsq);
							if (in->GetChild(4) != NULL) {
								RecurseForce(in->GetChild(4), dsq);
								if (in->GetChild(5) != NULL) {
									RecurseForce(in->GetChild(5), dsq);
									if (in->GetChild(6) != NULL) {
										RecurseForce(in->GetChild(6), dsq);
										if (in->GetChild(7) != NULL) {
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


static inline void ReadInput(char *filename) {
	double vx, vy, vz;
	FILE *f;

	f = fopen(filename, "r+t");
	if (f == NULL) {
		fprintf(stderr, "file not found: %s\n", filename);
		exit(-1);
	}

	fscanf(f, "%d", &nbodies);
	fscanf(f, "%d", &timesteps);
	fscanf(f, "%lf", &dtime);
	fscanf(f, "%lf", &eps);
	fscanf(f, "%lf", &tol);

	dthf = 0.5 * dtime;
	epssq = eps * eps;
	itolsq = 1.0 / (tol * tol);

	if (bodies == NULL) {
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
}

static inline void ComputeCenterAndDiameter(const int n, double &diameter, double &centerx, double &centery,
		double &centerz) {
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

static inline int min(int a, int b) {
	if (a < b)
		a = b;
	return a;
}

static void PrintDouble(double d) {
	int i;
	char str[16];

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

static OctTreeInternalNode *root;
static double gDiameter;

class ParallelForProcessor {
public:
	void operator()(const blocked_range<int>& range) const {
		for (int i = range.begin(); i != range.end(); i++) {
			bodies[i]->ComputeForce(root, gDiameter);
		}
	}
};

int main(int argc, char *argv[]) {
	//task_scheduler_init init;
	ParallelForProcessor parallelProcessor;

	fprintf(stderr, "\n");
	fprintf(stderr, "Lonestar benchmark suite\n");
	fprintf(stderr,
	"Copyright (C) 2007, 2008 The University of Texas at Austin\n");
	fprintf(stderr, "http://iss.ices.utexas.edu/lonestar/\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "application: BarnesHut v1.0\n");

	bodies = NULL;

	if (argc == 3)
		grainSize = atoi(argv[2]);

	timeval starttime, endtime;
	long runtime, lasttime, mintime;
	int run;

	runtime = 0;
	lasttime = -1;
	mintime = -1;
	run = 0;

	while (((run < 3) || (abs(lasttime - runtime) * 64 > min(lasttime, runtime))) && (run < 7)) {
		ReadInput(argv[1]);

		lasttime = runtime;
		gettimeofday(&starttime, NULL);

		for (step = 0; step < timesteps; step++) { // time-step the system
			double diameter, centerx, centery, centerz;
			ComputeCenterAndDiameter(nbodies, diameter, centerx, centery, centerz);

			OctTreeInternalNode *local_root = OctTreeInternalNode::NewNode(centerx, centery, centerz); // create the tree's root

			const double radius = diameter * 0.5;
//			for (int i = 0; i < nbodies; ++i) {
//				local_root->Insert(bodies[i], radius); // grow the tree by inserting each body
//			}
			local_root->InsertAll(bodies, nbodies, radius);

			int curr = 0;
			local_root->ComputeCenterOfMass(curr); // summarize subtree info in each internal node (plus restructure tree and sort bodies for performance reasons)

			root = local_root;
			gDiameter = diameter;

			if (grainSize != 0)
				parallel_for(blocked_range<int> (0, nbodies, grainSize), parallelProcessor);
			else
				parallel_for(blocked_range<int> (0, nbodies), parallelProcessor);

			OctTreeInternalNode::RecycleTree(); // recycle the tree

			for (int i = 0; i < nbodies; i++) { // the iterations are independent: they can be executed in any order and in parallel
				bodies[i]->Advance(); // advance the position and velocity of each body
			}
		} // end of time step

		gettimeofday(&endtime, NULL);
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

