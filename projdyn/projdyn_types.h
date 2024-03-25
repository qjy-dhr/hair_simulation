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


	using MatrixT = Eigen::Matrix<Scalar1, Rows, Cols, Options>; // A typedef of the dense matrix of Eigen.Eigen �ĳ��ܾ�������Ͷ��塣

	typedef MatrixT<Eigen::Dynamic, 3> Positions;				// һ�� n x 3 ��������ͨ������λ�û� m_velocities

	typedef MatrixT<Eigen::Dynamic, 1> Vector1;					// Vector of scalars		/ ����ʸ��

	typedef MatrixT<3, 1> Vector3;								// Vector of 3 scalars		/ 3 ��������ʸ��

	typedef MatrixT<2, 1> Vector2;								// Vector of 2 scalars		/ 2 ��������ʸ��

	typedef MatrixT<Eigen::Dynamic, Eigen::Dynamic> Matrix;		// Arbitrary size Matrix

	typedef Eigen::Matrix<int, Eigen::Dynamic, 3> Triangles;	// / ÿ�������������б�����������

	typedef Eigen::Matrix<int, 1, 3> Triangle;				// 	/ ���������Σ����ĸ���������

	typedef Eigen::Matrix<int, Eigen::Dynamic, 4>	Tetrahedrons;			//ÿ���ĸ��������б�����������	��

	typedef Eigen::Matrix<int, 1, 4> Tetrahedron;				//  ���������壨���ĸ�������


	//Sparse
	template<int Options = Eigen::ColMajor>
	using SparseMatrixT = Eigen::SparseMatrix<Scalar1, Options>;	// A typedef of the sparse matrix of Eigen. ����ϡ���������Ͷ��塣
	typedef SparseMatrixT<> SparseMatrix;						// A column-major sparse matrix.����ϡ�����
	typedef SparseMatrixT<Eigen::RowMajor> SparseMatrixRM;		// A row-major sparse matrix.������ϡ�����
	typedef Eigen::Triplet<Scalar1> Triplet;						// A triplet, used in the sparse triplet representation for matrices.һ����Ԫ�飬���ھ����ϡ����Ԫ���ʾ

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
