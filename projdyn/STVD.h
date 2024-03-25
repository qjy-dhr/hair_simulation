#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "projdyn_types.h"
#include <queue>

#define STVD_MODE_GRAPH_DIJKSTRA 0
#define STVD_MODE_ST_DIJKSTRA 1
#define STVD_DEFAULT_K 20

using namespace ProjDyn;

typedef std::pair<unsigned int, float> QueVert;

class QueVertComp
{
public:
	bool operator() (QueVert v1, QueVert v2) {
		return v1.second > v2.second;
	}
};

class STVD {

public:
	STVD();

	void init(Positions const& verts, Triangles const& tris = Triangles(0, 3), Tetrahedrons const& tets = Tetrahedrons(0, 4));

	STVD(Positions const& verts, Triangles const& tris = Triangles(0, 3), Tetrahedrons const& tets = Tetrahedrons(0, 4));

	/*
	Removes all source points from the mesh.
	*/
	void resetSources();

	/*
	Marks a vertex as a source, such that computeDistances() will yield the minimal geodesic distances
	to this vertex and all other marked vertices.
	*/
	void addSource(unsigned int vertexInd);

	/*
	Fills the vector distances with per-vertex values that correspond to the geodesic distance to the
	nearest marked vertex.
	*/
	void computeDistances(bool update = false, unsigned int mode = STVD_MODE_GRAPH_DIJKSTRA, Scalar1 maxDist = -1);

	/*
	Sets all distances to -1 (i.e. infinity), so even if computeDistances is run with update == true,
	the distances will be new.
	*/
	void resetDistances();

	/* 返回索引为vInd的向量到通过addSource设置的源的当前距离。需要通过computeDistances（）更新距离距离将是新的。*/
	double getDistance(unsigned int vInd);
	Vector1& getDistances();

	/* Set the parameter k which defines the size of the short term memory. Only
	   affects the process if the mode is set to STVD_MODE_ST_DIJKSTRA. */
	void setK(unsigned int k);

private:
	unsigned int m_numVerts;
	bool m_forTets;
	std::vector<unsigned int> m_sources;
	std::vector<std::vector<unsigned int>> m_neighbours;
	unsigned int m_k;
	unsigned int m_mode;
	bool m_isUpdated;
	Vector1 m_distances;
	Positions m_positions;
	Positions m_outerVertNormals;

	double updateVertDist(unsigned int v1, unsigned int v2, std::vector<int>& predecessors);
};