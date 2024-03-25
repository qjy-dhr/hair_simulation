#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <memory>

namespace ProjDyn {

	class Constraint;
	typedef std::shared_ptr<Constraint> ConstraintPtr;


	typedef double Scalar1;										// A scalar type, double or float, as defined
	typedef size_t Index;										// Type for array and vector sizes and indices.

	//Dense
	template < int Rows, int Cols, int Options = (Eigen::ColMajor) >


	using MatrixT = Eigen::Matrix<Scalar1, Rows, Cols, Options>; // A typedef of the dense matrix of Eigen.Eigen 的稠密矩阵的类型定义。

	typedef MatrixT<Eigen::Dynamic, 3> Positions;				// 一个 n x 3 标量矩阵，通常用于位置或 m_velocities

	typedef MatrixT<Eigen::Dynamic, 1> Vector1;					// Vector of scalars		/ 标量矢量

	typedef MatrixT<3, 1> Vector3;								// Vector of 3 scalars		/ 3 个标量的矢量

	typedef MatrixT<2, 1> Vector2;								// Vector of 2 scalars		/ 2 个标量的矢量

	typedef MatrixT<Eigen::Dynamic, Eigen::Dynamic> Matrix;		// Arbitrary size Matrix

	typedef Eigen::Matrix<int, Eigen::Dynamic, 3> Triangles;	// / 每行三个索引的列表，用于三角形

	typedef Eigen::Matrix<int, 1, 3> Triangle;				// 	/ 单个三角形（即四个索引）。

	typedef Eigen::Matrix<int, Eigen::Dynamic, 4>	Tetrahedrons;			//每行四个索引的列表，用于四面体	。

	typedef Eigen::Matrix<int, 1, 4> Tetrahedron;				//  单个四面体（即四个索引）


	//Sparse
	template<int Options = Eigen::ColMajor>
	using SparseMatrixT = Eigen::SparseMatrix<Scalar1, Options>;	// A typedef of the sparse matrix of Eigen. 本征稀疏矩阵的类型定义。
	typedef SparseMatrixT<> SparseMatrix;						// A column-major sparse matrix.列主稀疏矩阵。
	typedef SparseMatrixT<Eigen::RowMajor> SparseMatrixRM;		// A row-major sparse matrix.行优先稀疏矩阵。
	typedef Eigen::Triplet<Scalar1> Triplet;						// A triplet, used in the sparse triplet representation for matrices.一个三元组，用于矩阵的稀疏三元组表示

	//Quaternions
	using Quaternion = Eigen::Quaterniond;
    template < int Rows, int Cols>
	using ArrayQ = Eigen::Array<Quaternion, Rows, Cols>;
    typedef ArrayQ<Eigen::Dynamic, 1> Orientations;

    //Constants
    constexpr static float gravity = -9.81;
    constexpr static float E = 1; //Young's modulus
    constexpr static float cr_radius = 3;
    constexpr static float poisson = 0.377f;
    constexpr static float cr_unit_weight = 0.8f;
    constexpr static float cr_density = 1.3f;

}
