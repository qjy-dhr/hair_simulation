
#include <iostream>
#include <limits>

#include "STVD.h"

STVD::STVD()
{
}

Scalar1 triangleArea(Eigen::Vector3i tri, Positions& positions) {
	Scalar1 l1 = (positions.row(tri(0)) - positions.row(tri(1))).norm();
	Scalar1 l2 = (positions.row(tri(1)) - positions.row(tri(2))).norm();
	Scalar1 l3 = (positions.row(tri(2)) - positions.row(tri(0))).norm();
	Scalar1 p = (1.f / 2.f) * (l1 + l2 + l3);
	return std::sqrt(p * (p - l1) * (p - l2) * (p - l3));
}

Vector3 getTriangleNormal(Eigen::Vector3i& triangle, Positions& positions) {
	Vector3 normal = (positions.row(triangle(1)) - positions.row(triangle(0))).cross(positions.row(triangle(2)) - positions.row(triangle(0)));
	normal.normalize();
	return normal;
}


Positions getVertexNormals(Triangles& triangles, Positions& positions, int numOuterVertices)
{
	int numVerts = numOuterVertices;
	if (numVerts <= 0) numVerts = positions.rows();

	Positions vertNormals(numVerts, 3);
	vertNormals.setZero();

	for (int t = 0; t < triangles.rows(); t++) {
		Eigen::Vector3i tri = triangles.row(t);
		Vector3 tn = getTriangleNormal(tri, positions);
		Scalar1 area = triangleArea(tri, positions);
		for (int v = 0; v < 3; v++) {
			vertNormals.row(tri(v)) += tn * area;
		}
	}

	for (int v = 0; v < numVerts; v++) {
		vertNormals.row(v).normalize();
	}

	return vertNormals;
}

void STVD::init(Positions const& verts, Triangles const& tris, Tetrahedrons const& tets)
{
	m_numVerts = verts.rows();
	m_neighbours.resize(m_numVerts);
	m_mode = 0;
	m_k = STVD_DEFAULT_K;
	m_distances.resize(m_numVerts);
	m_distances.setConstant(-1);
	m_isUpdated = false;
	m_positions = verts;
	/* Establish neighbourhood structure for tets or triangles */
	if (tets.rows() > 0) {
		m_forTets = true;
		for (unsigned int tet = 0; tet < tets.rows(); tet++) {
			for (unsigned int locV1 = 0; locV1 < 4; locV1++) {
				int locV1Ind = tets(tet, locV1);
				for (unsigned int locV2 = 0; locV2 < 4; locV2++) {
					if (locV2 != locV1) {
						int locV2Ind = tets(tet, locV2);
						if (std::find(m_neighbours.at(locV1Ind).begin(), m_neighbours.at(locV1Ind).end(), locV2Ind) == m_neighbours.at(locV1Ind).end()) {
							m_neighbours.at(locV1Ind).push_back(locV2Ind);
						}
					}
				}
			}
		}
	}
	else {
		Triangles trisCopy = tris;
		Positions vertsCopy = verts;
		m_outerVertNormals = getVertexNormals(trisCopy, vertsCopy, -1);
		for (unsigned int tri = 0; tri < tris.rows(); tri++) {
			for (unsigned int locV1 = 0; locV1 < 3; locV1++) {
				int locV1Ind = tris(tri, locV1);
				for (unsigned int locV2 = 0; locV2 < 3; locV2++) {
					if (locV2 != locV1) {
						int locV2Ind = tris(tri, locV2);
						if (std::find(m_neighbours.at(locV1Ind).begin(), m_neighbours.at(locV1Ind).end(), locV2Ind) == m_neighbours.at(locV1Ind).end()) {
							m_neighbours.at(locV1Ind).push_back(locV2Ind);
						}
					}
				}
			}
		}
	}

}

STVD::STVD(Positions const& verts, Triangles const& tris, Tetrahedrons const& tets)
{
	init(verts, tris, tets);
}

void STVD::resetSources()
{
	m_sources.clear();
}

void STVD::addSource(unsigned int vertexInd)
{
	if (vertexInd < m_numVerts) {
		m_isUpdated = false;
		m_sources.push_back(vertexInd);
	}
}

void STVD::computeDistances(bool update, unsigned int mode, Scalar1 maxDist)
{
	m_mode = mode;

	/* Initialize distances, queue, predecessor list and "is final" status */
	if (!update) m_distances.setConstant(-1.);
	std::priority_queue<QueVert, std::vector<QueVert>, QueVertComp > queue;
	std::vector<int> predecessors;
	predecessors.resize(m_numVerts, -1);
	std::vector<bool> isFinal;
	isFinal.clear();
	isFinal.resize(m_numVerts, false);

	/* Add source vertices into queue */
	for (unsigned int v : m_sources) {
		m_distances[v] = 0;
		QueVert vd;
		vd.first = v;
		vd.second = 0.;
		queue.push(vd);
	}

	/* Main loop */
	while (!queue.empty()) {
		/* Extract vertex with shortest distance to sources */
		QueVert curV = queue.top();
		queue.pop();
		/* We check if vertex has already been finalized, since we don't check if vertices are
		already in the queue when inserting them. */
		if (!isFinal[curV.first]) {
			/* Vertex gets finalized. */
			isFinal[curV.first] = true;
			/* If this vertex is beyond maximal distance, do not add neighbours and stop distance evaluation */
			if (maxDist > 0 && m_distances[curV.first] > maxDist) {
				continue;
			}
			/* 对于所有未敲定的邻居。。。 */
			for (unsigned int nv : m_neighbours[curV.first]) {
				if (!isFinal[nv]) {
					/* ...如果新的距离更好，则更新它们的距离和前一个距离。 */
					double updatedDistance = updateVertDist(curV.first, nv, predecessors);
					if (m_distances[nv] < 0 || updatedDistance < m_distances[nv]) {
						/* Update distance, predecessors and push into queue. */
						m_distances[nv] = updatedDistance;
						predecessors[nv] = curV.first;
						QueVert vd;
						vd.first = nv;
						vd.second = updatedDistance;
						queue.push(vd);
					}
				}
			}
		}
	}

	m_isUpdated = true;
}

void STVD::resetDistances()
{
	m_distances.resize(m_numVerts);
	m_distances.setConstant(-1);
}

double STVD::getDistance(unsigned int vInd)
{
	if (!m_isUpdated) {
		std::cout << "Attempt to get out of date distances!" << std::endl;
		return -1;
	}
	if (vInd >= m_numVerts) {
		std::cout << "Attempt to get distance of illegal vertex!" << std::endl;
		return -1;
	}
	return m_distances(vInd);
}

Vector1& STVD::getDistances()
{
	if (!m_isUpdated) {
		std::cout << "Attempt to get out of date distances!" << std::endl;
	}
	return m_distances;
}

void STVD::setK(unsigned int k)
{
	m_k = k;
}

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

/* Returns the signed angle between v1 and v2 with respect to a normal n, i.e.
   the positive angle phi between them if v2 has to be rotated counterclockwise by phi to get to v1
   and the negative angle phi, if v2 has to be rotated clockwise by |phi| to get to v1*/
Scalar1 getSignedAngle(Vector3 v1, Vector3 v2, Vector3 n) {
	return std::acos(v1.dot(v2) / (v1.norm() * v2.norm())) * sgn(n.cross(v1).dot(v2));
}

double STVD::updateVertDist(unsigned int v1, unsigned int v2, std::vector<int>& predecessors)
{
	if (m_mode == STVD_MODE_GRAPH_DIJKSTRA) {
		return m_distances(v1) + (m_positions.row(v1) - m_positions.row(v2)).norm();
	}
	else if (m_mode == STVD_MODE_ST_DIJKSTRA) {
		// Implementation of the STVD update_dist function from [Campen et al. 2013]

		// This is the previous edge that was added to the sum, still in 3d
		Vector3 prevEdge = (m_positions.row(v2) - m_positions.row(v1));

		// Store actual predecessor temporarily and assume v1 would be predecessor
		// of v2
		int tmpPred = predecessors[v2];
		predecessors[v2] = v1;

		// The "front-most" vertex in the edge chain
		int prevV = v2;
		int curV = v1;

		Vector3 e_sum3d = prevEdge;
		Scalar1 eLength = e_sum3d.norm();

		// The following is only used for surface meshes
		// This will be the sum of the unfolded 2d vectors
		Vector2 e_sum2d;
		e_sum2d.setZero();
		e_sum2d(1) = eLength;
		// The angle between the last edge that has been unfolded and the x-axis
		Scalar1 curAngle = 0;

		// The currently best distance is the standard dijkstra update
		double bestDist = m_distances(v1) + eLength;

		// Now we check if unfolded edge sums can enhance this distance
		for (int i = 2; i <= m_k; i++) {
			if (predecessors[curV] < 0) break;

			// Get next edge in the chain
			int nextV = predecessors[curV];
			Vector3 nextEdge = (m_positions.row(curV) - m_positions.row(nextV));

			Scalar1 curDist = -1;
			// If both edges run along the surface, unfold them into the common plane
			// defined by the vertex normal
			if (!m_forTets) {
				Vector3 n = m_outerVertNormals.row(curV);
				Vector3 nextEdgeFlat = nextEdge - n * (n.dot(nextEdge));
				Vector3 prevEdgeFlat = prevEdge - n * (n.dot(prevEdge));
				Scalar1 angle = getSignedAngle(prevEdgeFlat, nextEdgeFlat, n);
				curAngle += angle;
				Scalar1 l = nextEdgeFlat.norm();
				Vector2 nextEdge2d(std::sin(curAngle) * l, std::cos(curAngle) * l);
				e_sum2d += nextEdge2d;
				curDist = m_distances[nextV] + e_sum2d.norm();
			}
			else {
				e_sum3d += nextEdge;
				curDist = m_distances[nextV] + e_sum3d.norm();
			}
			if (curDist >= 0 && curDist < bestDist) bestDist = curDist;

			prevV = curV;
			curV = nextV;
			prevEdge = nextEdge;
		}

		return bestDist;
	}
	return m_distances(v1) + (m_positions.row(v1) - m_positions.row(v2)).norm();
}
