#pragma once
#pragma once

#define M_PI 3.14159265358979323846

#define	USE_QUARTIC_POL true

#include <math.h>
#include "projdyn_types.h"
#include "Edge.h"

#include "PDConstraint.h"

#include<Eigen/Geometry>
#include <memory>

#include <iostream>
#include <limits>
#include "STVD.h"

namespace ProjDyn {
    typedef Eigen::SimplicialLDLT<SparseMatrix> SparseSolver;
    // typedef Eigen::SimplicialLDLT<Matrix> m_subspaceSolver;
    typedef Eigen::LLT<Matrix> m_subspaceSolver;

    //模拟器类，用于处理以特征矩阵表示的对象，以及上面的约束来运行投影动力学模拟。
    class Simulator {
    public:
        Simulator()
        {
            //构造模拟器对象
            m_hasGrab = false;
            is_ready = false;
            is_empty = true;
            is_simulating = false;
            use_cosserat_rods = false;
            // use_subspace = true;
            use_subspace = true;
            cr_wind = false;
            m = 0;
        }

        void setMesh(Positions& pos) {
            // 传递几何体的位置
            m = pos.rows();
            m_positions_init = pos;
            m_positions = pos;
        }

        void createPositionSubspace(unsigned int numSamples, bool useSkinningSpace)
        {
            getSamples(numSamples);

            numSamples = samples.size();
            std::sort(samples.begin(), samples.end());
            samples.erase(std::unique(samples.begin(), samples.end()), samples.end());

            Scalar1 furthestDist = getSampleDiameter(samples);
            Scalar1 r = furthestDist * m_baseFunctionRadius;

            m_baseFunctionWeights = getRadialBaseFunctions(samples, true, r);

            /* if (useSkinningSpace) {
                 bool isFlat = false;
                 if (m_positions.col(2).norm() < 1e-10) isFlat = true;
                 m_baseFunctions = createSkinningSpace(m_positions, m_baseFunctionWeights, nullptr, 1U, nullptr, nullptr, isFlat);
             }
             else {*/
            m_baseFunctions = m_baseFunctionWeights;
            int ind = m_baseFunctions.cols();
            m_baseFunctions.conservativeResize(m_baseFunctions.rows(), m_baseFunctions.cols() + 3);
            m_baseFunctions.col(ind) = mom_positions.col(0);
            m_baseFunctions.col(ind + 1) = mom_positions.col(1);
            m_baseFunctions.col(ind + 2) = mom_positions.col(2);
            // }
        }

        Matrix getRadialBaseFunctions(std::vector<unsigned int>& samples, bool partitionOfOne, double r, double eps = -1., int numSmallSamples = -1, double smallSampleRadius = 1.)
        {
            if (eps < 0) eps = std::sqrt(-std::log(0.0001)) / r;

            Scalar1 pi = 3.14159265358979323846;
            int nSamples = samples.size();
            Matrix baseFunctions;
            baseFunctions.setZero(m_nVertices, nSamples);
            Scalar1 a = (1. / std::pow(r, 4.));
            Scalar1 b = -2. * (1. / (r * r));

            for (int i = 0; i < nSamples; i++) {
                if (numSmallSamples > 0 && i > nSamples - numSmallSamples) {
                    r = smallSampleRadius;
                    eps = std::sqrt(-std::log(0.0001)) / r;
                }

                unsigned int curSample = samples[i];
                clearSources();
                addSource(curSample);
                computeDistances(r, false);

                for (int v = 0; v < m_nVertices; v++) {
                    double curDist = getDistance(v);
                    double val = 0;
                    if (curDist < 0) {
                        val = 0;
                    }
                    else if (USE_QUARTIC_POL) {
                        if (curDist >= r) {
                            val = 0;
                        }
                        else {
                            val = a * std::pow(curDist, 4.) + b * (curDist * curDist) + 1;
                        }
                    }
                    else {
                        val = std::exp(-(curDist * eps * curDist * eps));
                        if (val < 0.0001) val = 0;
                    }
                    baseFunctions(v, i) = val;
                }

            }

            if (partitionOfOne) {
                for (int v = 0; v < m_nVertices; v++) {
                    Scalar1 sum = baseFunctions.row(v).sum();
                    if (sum < 1e-6) {
                        std::cout << "Warning: a vertex isn't properly covered by any of the radial basis functions!" << std::endl;
                        baseFunctions(v, indexOfMaxColCoeff(baseFunctions, v)) = 1.;
                    }
                    else {
                        baseFunctions.row(v) /= sum;
                    }
                }
            }

            return baseFunctions;
        }

        int indexOfMaxColCoeff(Matrix& A, unsigned int row)
        {
            int j = 0;
            double max = 0;
            for (int i = 0; i < A.cols(); i++) {
                double curNorm = std::abs(A(row, i));
                if (curNorm > max) {
                    max = curNorm;
                    j = i;
                }
            }
            return j;
        }

        void getSamples(unsigned int numSamples)
        {
            //std::vector<unsigned int> samples;

            clearSources();

            //unsigned int firstVert = rand() % m_nVertices;
            unsigned int firstVert = 0;
            samples.push_back(firstVert);
            addSource(firstVert);

            addSamples(numSamples - 1);

        }

        void clearSources()
        {
            m_stvd.resetSources();
            m_stvd.resetDistances();
        }

        void addSource(unsigned int index)
        {
            m_stvd.addSource(index);
            m_upToDate = false;
        }

        void addSamples(unsigned int numSamples) {
            if (numSamples <= 0) return;
            Matrix heatHistory(m_nVertices, numSamples);
            heatHistory.setZero();
            for (int i = 0; i < numSamples - 1; i++) {
                computeDistances();
                heatHistory.col(i) = m_distances;
                double maxDist = 0;
                unsigned int bestVert = -1;
                for (int v = 0; v < m_nVertices; v++) {
                    double curDist = getDistance(v);
                    if (maxDist < curDist) {
                        bestVert = v;
                        maxDist = curDist;
                    }
                }
                if (bestVert == -1) {
                    std::cout << "Error during sampling, returning with fewer samples... " << std::endl;
                    return;
                }
                if (std::find(samples.begin(), samples.end(), bestVert) != samples.end()) {
                    std::cout << "Duplicate vertex was selected in sampling, canceling." << std::endl;
                    return;
                }
                addSource(bestVert);
                samples.push_back(bestVert);
            }
        }

        double getDistance(unsigned int index)
        {
            if (!m_upToDate) {
                std::cout << "Sources have been added and computeDistances() needs to be called first (or computeDistances() has never been called)!" << std::endl;
                return 0;
            }
            return m_distances(index);
        }

        void computeDistances(Scalar1 maxDist = -1, bool accurate = false)
        {
            if (accurate) {
                m_stvd.computeDistances(true, STVD_MODE_ST_DIJKSTRA, maxDist);
            }
            else {
                m_stvd.computeDistances(true, STVD_MODE_GRAPH_DIJKSTRA, maxDist);
            }
            m_distances = m_stvd.getDistances();

            m_upToDate = true;
        }

        Scalar1 getSampleDiameter(std::vector<unsigned int>& samples) {
            clearSources();
            for (int v = 0; v < samples.size(); v++) addSource(samples[v]);
            computeDistances();
            double furthestDist = 0;
            for (unsigned int v = 0; v < m_distances.rows(); v++) {
                if (getDistance(v) > furthestDist) {
                    furthestDist = getDistance(v);
                }
            }

            return furthestDist;
        }

        void finalizeBaseFunctions() {
            m_baseFunctionsTransposed = m_baseFunctions.transpose();
            m_baseFunctionsSquared = m_baseFunctionsTransposed * m_baseFunctions;
        }

        //创建杆
        void setRods(const std::vector<Positions> rods) {

            assert(rods.size() >= 1); //至少有一个
            is_empty = false;
            is_ready = false;
            use_cosserat_rods = true;
            rod_indices.clear();
            cr_num_positions = 0;
            cr_positions_init = Vector1();    cr_positions_init_3d = Matrix(0,3);
            

            rod_size = rods.size();
            Index rod_index = 0;
            for (Index i = 0; i < rods.size(); i++) {
                Positions pos = rods.at(i);
                assert(pos.rows() >= 2); //至少有一段
                rod_indices.push_back(rod_index);

                //假定所有节段有相同的长度，杆的长度定义为第一节段的长度
                cr_segments_length.push_back((pos.row(1) - pos.row(0)).norm());
                rod_index += pos.rows();
                cr_num_positions += pos.rows();

                //更新位置
                Vector1 new_pos = pos2flatVec(pos);                                    
                Vector1 temp_pos(cr_positions_init.size() + new_pos.size());           
                temp_pos << cr_positions_init, new_pos;                                
                cr_positions_init = temp_pos;      

                Matrix temp_pos_3d(cr_positions_init_3d.rows() + pos.rows(), 3);
                temp_pos_3d << cr_positions_init_3d, pos;
                cr_positions_init_3d = temp_pos_3d;

               
                      
            }
            cr_positions = cr_positions_init;          
            
            cr_positions_3d = cr_positions_init_3d;


            

            //18=3*(8-2)
            cr_num_flat_pos = 3 * (cr_num_positions - rod_indices.size());
            //24=3*8
            cr_num_coord = 3 * cr_num_positions;
            cr_num_quaternions = cr_num_positions - rod_indices.size();
            cr_size = 3 * cr_num_positions + 4 * cr_num_quaternions;
            cr_velocities.setZero(cr_num_coord);
            cr_angular_velocities.setZero(cr_num_flat_pos);

            //3d
            cr_size_3d = cr_num_positions + cr_num_quaternions;
            cr_velocities_3d.setZero(cr_num_positions,3);
            cr_angular_velocities_3d.setZero(cr_num_quaternions,3);
            
            std::cout << "Successfully added " << rods.size() << " rods" << std::endl;

        }

        //Simple scene with only rods fixed to their initial placement
        void addDefaultConstraints() {
            cr_constraints.clear();
            addFixedPosConstraint();
            addSSConstraints();
            addBTConstraints();
        }

        //Rods stuck to a ball
        void addMeshConstraints() {
            cr_constraints.clear();
            addSSConstraints(100);
            addBTConstraints(10);
            addMovingPointConstraints(1e5, 1e4);
            addSphereCollisionConstraints(10, 2);
            addSelfCollisionConstraints(100, 2);
        }

        void resetPositions() {
            //将顶点位置重置为初始状态，并将m_velocities设置为0
            m_positions = m_positions_init;
            if (use_cosserat_rods) {
                cr_positions = cr_positions_init;
                cr_velocities.setZero(3 * cr_num_positions);
                cr_angular_velocities.setZero(cr_num_flat_pos);
                cr_orientations = quatFromPos(cr_positions);
                initializeSystem(m_default_constraints);
            }
        }

        bool isInitialized() {
            return is_ready && !is_empty;
        }

        void initializeLHS() {
            //Compute left-hand side of system (can be done earlier)
            cr_lhs.resize(cr_size, cr_size);
            cr_lhs.setZero();
            cr_lhs = cr_M_star / (h * h);

            for (Index i = 0; i < cr_constraints.size(); i++) {
                auto c = cr_constraints.at(i);
                auto Ai = c->getAiMatrix();
                cr_lhs += c->getSelectionMatrixWeighted() * Ai.transpose() * Ai * c->getSelectionMatrix();
            }
            if (use_subspace) {

                /* m_subspaceLHS_mom = m_baseFunctionsTransposed * m_lhsMatrix * m_baseFunctions * (m_timeStep * m_timeStep);
                 m_subspaceLHS_inner = m_baseFunctionsTransposed * conMat * m_baseFunctions;
                 m_lhsMatrixSampled = m_subspaceLHS_mom * (1. / (m_timeStep * m_timeStep)) + m_subspaceLHS_inner;

                 m_denseSolver.compute(m_lhsMatrixSampled);*/

                 //std::cout << cr_lhs << " ";

                // Matrix LHS = cr_lhs*m_baseFunctions.sparseView();
                 //sub_solver.compute(LHS);
                cr_solver.compute(cr_lhs);
            }
            else {
                cr_solver.compute(cr_lhs);
            }
        }

        bool initializeSystem(bool default_constraints) {
            // Setup simulation (collect constraints, set up and factorize global step system)
            //设置模拟（收集约束、设置和分解全局步长系统）
            if (is_simulating || is_empty) {
                return false;
            }

            m_default_constraints = default_constraints;

            if (use_cosserat_rods) {
                /****************
                 * CR Variables *
                 ****************/


                /*for (int i = 0; i < m_baseFunctions.rows(); i++) {
                    for (int j = 0; j < m_baseFunctions.size() / m_baseFunctions.rows(); j++) {
                        std::cout << m_baseFunctions(i, j) << " ";
                    }
                    std::cout << std::endl;
                }*/



                //Set default external forces
                cr_f_ext.setZero(cr_num_coord);
                for (Index i = 0; i < cr_num_coord; i++) {
                    cr_f_ext.coeffRef(i) = 0;
                    cr_f_ext.coeffRef(++i) = 0; //gravity;
                    cr_f_ext.coeffRef(++i) = 0;
                }
                //3d
                cr_f_ext_3d.setZero(cr_num_positions,3);
                for (Index i = 0; i < cr_num_positions; i++) {
                    for (int j = 0; j < 3; j++) {
                        cr_f_ext_3d.coeffRef(i,j) = 0;
                    }
                }

                cr_orientations = quatFromPos(cr_positions);
                cr_torques.setZero(cr_num_flat_pos);

                
                Vector1 masses_flat = Vector1::Ones(cr_num_coord) * cr_unit_weight;
                cr_masses = masses_flat.asDiagonal();
                cr_masses_inv = cr_masses.cwiseInverse();

                //3d ?????
                cr_torques_3d.setZero(cr_num_quaternions, 3);

                Vector1 masses_flat_3d = Vector1::Ones(cr_num_positions) * cr_unit_weight;
                cr_masses_3d = masses_flat_3d.asDiagonal();
                cr_masses_inv_3d = cr_masses_3d.cwiseInverse();


                //Build J matrix
                const float j_val = M_PI * cr_radius / 4;
                Vector1 J_flat_vec;
                J_flat_vec.setZero(3);
                J_flat_vec << j_val, j_val, j_val * 2;
                J_flat_vec *= cr_density;
                cr_J_vec = stackDiagonal(J_flat_vec, cr_num_positions - rod_indices.size());
                Vector1 J_flat_vec_inv;
                J_flat_vec_inv.setZero(3);
                J_flat_vec_inv << 1 / j_val, 1 / j_val, 1 / (j_val * 2);
                J_flat_vec_inv /= cr_density;
                cr_J_vec_inv = stackDiagonal(J_flat_vec_inv, cr_num_positions - rod_indices.size());

                Vector1 J_flat_quat;
                J_flat_quat.setZero(4);
                J_flat_quat << 0, j_val, j_val, j_val * 2;
                J_flat_quat *= cr_density;
                cr_J_quat = stackDiagonal(J_flat_quat, cr_num_quaternions);
                Vector1 J_flat_quat_inv;
                J_flat_quat_inv.setZero(4);
                J_flat_quat_inv << 0, 1 / j_val, 1 / j_val, 1 / (j_val * 2);
                J_flat_quat_inv /= cr_density;
                cr_J_quat_inv = stackDiagonal(J_flat_quat_inv, cr_num_quaternions);

                //按相应的段长度缩放
                for (Index i = 0; i < rod_indices.size(); i++) {
                    Index rod_index = rod_indices.at(i);
                    Index next_index = i == rod_indices.size() - 1 ? cr_num_positions : rod_indices.at(i + 1);
                    for (Index j = rod_index; j < next_index; j++) {
                        float seg_len = cr_segments_length.at(i);
                        cr_J_vec.row(j) *= seg_len;
                        cr_J_vec_inv.row(j) /= seg_len;
                        cr_J_quat.row(j) *= seg_len;
                        cr_J_quat_inv.row(j) /= seg_len;
                    }
                }

                //Build M_star matrix
                cr_M_star.resize(cr_size, cr_size);
                cr_M_star.setZero();
                for (Index i = 0; i < cr_masses.rows(); i++) {
                    cr_M_star.row(i) = cr_masses.row(i);
                }
                for (Index i = 0; i < cr_J_quat.rows(); i++) {
                    for (Index j = 0; j < cr_J_quat.cols(); j++) {
                        cr_M_star.coeffRef(i + cr_num_coord, j + cr_num_coord) = cr_J_quat.coeff(i, j);
                    }
                }



                /*if (m_default_constraints) {
                    addDefaultConstraints();
                } else {*/
                addMeshConstraints();
                // }

                std::cout << "Number of constraints: " << cr_constraints.size() << std::endl;
                //预先计算左侧矩阵
                initializeLHS();

            }

            is_ready = true;
            return is_ready;
        }

        bool step(int num_iterations) {
            if (!isInitialized()) return false;

            is_simulating = true;

            const Index num_constraints = cr_constraints.size();
            std::vector<Vector1> projections(num_constraints);

            //sx(t) = x(t) +hv(t) +h2M−1fext
            Vector1 s_x = cr_positions + h * cr_velocities + h * h * cr_masses_inv * cr_f_ext;
            Vector1 cross = vectorCross(cr_angular_velocities, cr_J_vec * cr_angular_velocities);

            Vector1 s_w = cr_angular_velocities + h * cr_J_vec_inv * (cr_torques - cross);
            Orientations s_w_u = pos2quat(s_w);
            Orientations s_u = pos2quat(quat2pos(cr_orientations) + h / 2 * quat2pos(cr_orientations * s_w_u));

            Vector1 s_t = posQuatConcat(s_x, s_u);
            cr_q_t = s_t;

            Index step = 0;

            if (use_subspace) {
                while (step < num_iterations) {
                    step++;
                    /****************
                     ** Local step **
                     ****************/
#pragma omp parallel for
                    for (int i = 0; i < num_constraints; i++) {
                        auto proj = cr_constraints.at(i)->projectOnConstraintSet(cr_q_t);
                        projections.at(i) = proj;
                    }

                    /*****************
                     ** Global step **
                     *****************/

                     //Compute right-hand side
                    cr_rhs.setZero();
                    cr_rhs = cr_M_star * s_t / (h * h);
                    for (Index i = 0; i < num_constraints; i++) {
                        auto c = cr_constraints.at(i);
                        cr_rhs += c->getLHS() * projections.at(i);
                    }

                    cr_q_t = cr_solver.solve(cr_rhs);

                    if (cr_solver.info() == Eigen::Success) {
                        //Update velocities and angular velocities 更新速度和角速度
                        Orientations new_quat;
                        Vector1 new_pos;
                        separatePosQuat(&cr_q_t, new_pos, new_quat, cr_num_coord);
                        cr_velocities = (new_pos - cr_positions) / h;
                        cr_positions = new_pos;
                        cr_angular_velocities = 2 / h * quat2pos(conjugateQuat(cr_orientations) * new_quat);
                        cr_orientations = new_quat;
                    }
                    else {
                        std::cout << "Unable to solve rods constraints" << std::endl;
                        return false;
                    }
                }
            }
            else {
                while (step < num_iterations) {
                    step++;
                    /****************
                     ** Local step **
                     ****************/
#pragma omp parallel for
                    for (int i = 0; i < num_constraints; i++) {
                        auto proj = cr_constraints.at(i)->projectOnConstraintSet(cr_q_t);
                        projections.at(i) = proj;
                    }

                    /*****************
                     ** Global step **
                     *****************/
                     //Compute right-hand side
                    cr_rhs.setZero();
                    cr_rhs = cr_M_star * s_t / (h * h);
                    for (Index i = 0; i < num_constraints; i++) {
                        auto c = cr_constraints.at(i);
                        cr_rhs += c->getLHS() * projections.at(i);
                    }

                    cr_q_t = cr_solver.solve(cr_rhs);

                    if (cr_solver.info() == Eigen::Success) {
                        //Update velocities and angular velocities 更新速度和角速度
                        Orientations new_quat;
                        Vector1 new_pos;
                        separatePosQuat(&cr_q_t, new_pos, new_quat, cr_num_coord);
                        cr_velocities = (new_pos - cr_positions) / h;
                        cr_positions = new_pos;
                        cr_angular_velocities = 2 / h * quat2pos(conjugateQuat(cr_orientations) * new_quat);
                        cr_orientations = new_quat;
                    }
                    else {
                        std::cout << "Unable to solve rods constraints" << std::endl;
                        return false;
                    }
                }
            }

            //Update external forces if enabled
            updateWind();

            is_simulating = false;
            return true;
        }

        void setGrab(const std::vector<Index>& grabVerts, const std::vector<Eigen::Vector3f>& grabPos) {
            // 如果执行了 alt-click 并且
             // 顶点被抓取。
             // 你可以使用抓取的顶点列表和新的位置来强制顶点移动
             // 在你的 step() 方法中。

            m_grabVerts = grabVerts;
            m_grabPos = grabPos;
            m_hasGrab = true;

            assert(grabVerts.size() >= 1);
            Vector3 p = m_positions.row(grabVerts.at(0));
            Eigen::Vector3f d = grabPos.at(0);
            Vector3 f(d.x(), d.y(), d.z());
            Vector3 diff = f - p;
#pragma omp parallel for
            for (int i = 0; i < m_positions.rows(); i++) {
                m_positions.row(i) += diff;
            }
        }

        void releaseGrab() {
            m_hasGrab = false;
        }

        void setTimestep(const float timestep) {
            h = timestep;
        }

        void toggleWind() {
            if (cr_wind) {
                cr_f_ext.setZero();
            }
            cr_wind = !cr_wind;
        }

        Positions* getPositions() {
            return &m_positions;
        }

        Positions* getRodsPositions() {
            upload_pos = vec2pos(cr_positions);
            // upload_pos = m_baseFunctions * upload_pos;

            /*尝试加入点
            upload_pos_end.resize(upload_pos.rows() / 3 * 5, 3);
             for (int i = 0; i < upload_pos.rows() / 3; i++) {
                 for (int j = 0; j < 3; j++) {
                     upload_pos_end(i * 3,j) = upload_pos(i * 3,j);
                     upload_pos_end(i * 3 + 1,j) = 0;
                     upload_pos_end(i * 3 + 2,j) = upload_pos(i * 3 + 1,j);
                     upload_pos_end(i * 3 + 3,j) = 0;
                     upload_pos_end(i * 3 + 4,j) = upload_pos(i * 3 + 2,j);
                 }

             }*/

            return &upload_pos;
        }

        std::vector<Index> getRodIndices() {
            return rod_indices;
        }

        Positions* getRodsTangents() {
            upload_tan.resize(cr_num_positions, 3);
            Vector3 t = Vector3::UnitY();
            Index j = 0;
            for (Index i = 0; i < cr_num_quaternions; i++) {
                Quaternion q = cr_orientations[i];
                upload_tan.row(j++) = q.normalized().toRotationMatrix() * t;
                if (i != 0 && std::find(rod_indices.begin(), rod_indices.end(), i) != rod_indices.end()) {
                    upload_tan.row(j++) = q.normalized().toRotationMatrix() * t;
                }
            }
            upload_tan.row(upload_tan.rows() - 1) = upload_tan.row(upload_tan.rows() - 2);

            //upload_tan = m_baseFunctions * upload_tan;

            return &upload_tan;
        }

        Positions* getRodsNormals() {
            upload_normals.resize(cr_num_positions, 3);
            Vector3 n = -Vector3::UnitZ();
            Index j = 0;
            for (Index i = 0; i < cr_num_quaternions; i++) {
                Quaternion q = cr_orientations[i];
                upload_normals.row(j++) = q.normalized().toRotationMatrix() * n;
                if (i != 0 && std::find(rod_indices.begin(), rod_indices.end(), i) != rod_indices.end()) {
                    upload_normals.row(j++) = q.normalized().toRotationMatrix() * n;
                }
            }
            upload_normals.row(upload_normals.rows() - 1) = upload_normals.row(upload_normals.rows() - 2);
            // upload_normals = m_baseFunctions * upload_normals;
            return &upload_normals;
        }


    public:
        bool m_hasGrab = false;
        std::vector<Index> m_grabVerts;
        std::vector<Eigen::Vector3f> m_grabPos;
        std::vector<Edge> edges;
        double h = 0.05; //Simulation step size

        /******************
         * Mesh variables *
         ******************/
        Index m;
        Positions m_positions;
        Positions m_positions_init;

        /***************************
         * Cosserat Rods variables *
         ***************************/
        std::vector<float> cr_segments_length;
        std::vector<Index> rod_indices;

        //位置变量（维度：cr_num_coord）。
        Index cr_size; //平面变量的完整大小
        Vector1 cr_positions_init;
        Vector1 cr_positions;

        //3d
        Index cr_size_3d;
        Matrix cr_positions_init_3d;
        Matrix cr_positions_3d;

        Vector1 cr_q_t; //保存所有变量（位置和四元数）的平面向量。
        Index cr_num_coord;
        Index cr_num_positions; //Number of positions
        Vector1 cr_velocities;    
        Vector1 cr_f_ext;         
        SparseMatrix cr_masses;
        SparseMatrix cr_masses_inv;

        //3d
        Matrix cr_q_t_3d;
        Matrix cr_velocities_3d;
        Matrix cr_f_ext_3d;
        SparseMatrix cr_masses_3d;
        SparseMatrix cr_masses_inv_3d;
        

        //角度变量（维：cr_num_flat_pos）
        Index cr_num_flat_pos;
        Vector1 cr_angular_velocities;  Matrix cr_angular_velocities_3d;
        Vector1 cr_torques;             Matrix cr_torques_3d;
        SparseMatrix cr_J_vec;


        SparseMatrix cr_J_vec_inv; //TODO: 使用对角矩阵类型
        SparseMatrix cr_J_quat;
        SparseMatrix cr_J_quat_inv;
        SparseMatrixRM cr_M_star;
        
        //3d
        SparseMatrix cr_J_vec_inv_3d; 
        SparseMatrix cr_J_quat_3d;
        SparseMatrix cr_J_quat_inv_3d;
        SparseMatrixRM cr_M_star_3d;

        Index cr_num_quaternions;
        Orientations cr_orientations;

        std::vector<CRConstraint*> cr_constraints;
        SparseSolver cr_solver;

        Vector1 cr_rhs;
        SparseMatrix cr_lhs;

        //状态变量。
        bool is_ready, is_empty;
        bool is_simulating;
        bool use_cosserat_rods;
        bool m_default_constraints = true;
        bool cr_wind;
        float time = 0.0f;
        Positions upload_pos, upload_tan, upload_normals;

        //子空间
        int rod_size;
        bool use_subspace;
        std::vector< unsigned int > samples;
        STVD m_stvd;
        Vector1 m_distances;
        bool m_upToDate = false;
        Matrix m_baseFunctionWeights;
        Matrix m_baseFunctions;   //N*4k
        Matrix m_baseFunctionsTransposed;
        Matrix m_baseFunctionsSquared;
        double m_baseFunctionRadius = 1.1; //这个数字越大，基本函数的支持就越大。
        int m_nVertices;
        Positions mom_positions;
        m_subspaceSolver sub_solver;
        m_subspaceSolver sub_solver_down;


        Positions upload_pos_end;

        //拉伸和剪切
        void addSSConstraints(Scalar1 weight = 1.0) {
            for (Index ind = 0; ind < rod_indices.size(); ind++) {
                Index rod_index = rod_indices.at(ind);
                Index next_index = ind == rod_indices.size() - 1 ? cr_num_positions : rod_indices.at(ind + 1);

                for (Index i = rod_index; i < next_index - 1; i++) {
                    auto new_c = new StretchShearConstraint(cr_size, weight, i * 3,
                        cr_num_coord + (i - ind) * 4, cr_segments_length.at(ind));
                    cr_constraints.push_back(new_c);
                }
            }
        }

        //弯曲和扭曲
        void addBTConstraints(Scalar1 weight = 1.0) {
            for (Index ind = 0; ind < rod_indices.size(); ind++) {
                Index rod_index = rod_indices.at(ind);
                Index next_index = ind == rod_indices.size() - 1 ? cr_num_positions : rod_indices.at(ind + 1);
                for (Index i = rod_index; i < next_index - 2; i++) {
                    auto new_c = new BendTwistConstraint(cr_size, weight, cr_num_coord + (i - ind) * 4, cr_segments_length.at(ind));
                    cr_constraints.push_back(new_c);
                }
            }
        }

        //定点
        void addFixedPosConstraint() {
            for (auto i : rod_indices) {
                Vector3 p;
                p << cr_positions.coeff(i * 3), cr_positions.coeff(i * 3 + 1), cr_positions.coeff(i * 3 + 2);
                auto fpc = new FixedPointConstraint(cr_size, 10000, i * 3, p);
                cr_constraints.push_back(fpc);
            }
        }
        //移动位置
        void addMovingPointConstraints(float moving_weight = 1e5, float normal_weight = 10) {
            if (m_positions.size() == 0) return;

            Vector3 center = computeCenter(m_positions);

            Index ind = 0;
            for (auto i : rod_indices) {
                //Find closest point to the rod on the mesh
                Index index = 0;
                float min_distance = std::numeric_limits<float>::max();
                Vector3 closest_point;
                for (Index j = 0; j < m_positions.rows(); j++) {
                    Vector3 p = m_positions.row(j);
                    Vector3 r(cr_positions.coeff(i * 3), cr_positions.coeff(i * 3 + 1), cr_positions.coeff(i * 3 + 2));
                    float dist = (p - r).norm();
                    if (dist < min_distance) {
                        index = j;
                        min_distance = dist;
                        closest_point = r;
                    }
                }

                auto mpc = new MovingPointConstraint(cr_size, moving_weight, i * 3, &m_positions, index);
                cr_constraints.push_back(mpc);

                Vector3 normal = (center - closest_point).normalized();
                auto nc = new NormalConstraint(cr_size, normal_weight, i * 3, normal);
                //cr_constraints.push_back(nc);

                //Quaternion q_normal = cr_orientations.coeff(i);
                Quaternion q_normal = Quaternion::FromTwoVectors(-Vector3::UnitY(), normal);
                auto ncq = new NormalConstraintQuat(cr_size, normal_weight, cr_num_coord + (i - ind) * 4, q_normal);
                cr_constraints.push_back(ncq);

                ind++;
            }
        }

        //球体碰撞
        void addSphereCollisionConstraints(Scalar1 weight = 1.0, Scalar1 forceFactor = 10.0) {
            if (m_positions.size() == 0) return;

            //Because the mesh is a ball every point is at the same distance to the center
            const Index default_index = 0;

            Vector3 center = computeCenter(m_positions);
            Vector3 random_point = m_positions.row(default_index);
            Scalar1 radius = (center - random_point).norm();

            for (Index ind = 0; ind < rod_indices.size(); ind++) {
                Index rod_index = rod_indices.at(ind);
                Index next_index = ind == rod_indices.size() - 1 ? cr_num_positions : rod_indices.at(ind + 1);


                for (Index i = rod_index; i < next_index; i++) {
                    auto scc = new SphereConstraint(cr_size, weight, radius * 1.2, i * 3, center, &m_positions, default_index, forceFactor);
                    cr_constraints.push_back(scc);
                }
            }
        }

        //自碰撞
        void addSelfCollisionConstraints(Scalar1 weight = 1.0, Scalar1 forceFactor = 1.0) {
            if (m_positions.size() == 0) return;

            const Index default_index = 0;

            Vector3 center = computeCenter(m_positions);
            Vector3 random_point = m_positions.row(default_index);
            Scalar1 radius = (center - random_point).norm();

            //For each index of each rod
            for (Index ind = 0; ind < rod_indices.size(); ind++) {
                Index rod_index = rod_indices.at(ind);
                Index next_index = ind == rod_indices.size() - 1 ? cr_num_positions : rod_indices.at(ind + 1);

                //TODO: Only setup these constraints for a subset of the rods (neighboring rods for example)
                for (Index i = rod_index; i < next_index - 1; i++) {
                    //Then for each other rod create a collision constraint
                    for (Index j = 0; j < rod_indices.size(); j++) {
                        if (j == ind) continue;
                        Index coll_curr_index = rod_indices.at(j);
                        Index coll_next_index = j == rod_indices.size() - 1 ? cr_num_positions : rod_indices.at(j + 1);
                        for (Index coll_index = coll_curr_index; coll_index < coll_next_index; coll_index++) {
                            auto scc = new SelfCollisionConstraint(cr_size, weight, cr_radius * 0.2, coll_index * 3, i * 3, (i + 1) * 3, forceFactor);
                            cr_constraints.push_back(scc);
                        }
                    }
                }
            }
        }

        void updateWind() {
            if (!isInitialized()) return;
            if (cr_wind) {
                time += h;
                float forceX = 100.0;
                float forceY = 30.0;
                float factorIntensity = 0.5;
                float factorY = 0.2;
                float factorX = 0.9;

                //Make the intensity vary over time
                forceX *= sin(time * factorIntensity);

                //Using the sphere vector representation we
                // can make the force rotate
                Vector3 vec(
                    forceY * cos(time * factorY),
                    forceX * sin(time * factorY) * cos(time * factorX),
                    forceX * sin(time * factorY) * sin(time * factorX)
                );

                //Set the external forces to the wind force
#pragma omp parallel for
                for (int i = 0; i < cr_num_coord; i += 3) {
                    cr_f_ext.coeffRef(i) = vec.x();
                    cr_f_ext.coeffRef(i + 1) = clamp(vec.y(), 0, forceX);
                    cr_f_ext.coeffRef(i + 2) = vec.z();
                }
            }
        }

        //从向量创建四元数的虚部。实部为0。
        static Orientations pos2quat(const Vector1 p) {
            //The size of the vector must a multiple of 3
            assert(p.size() % 3 == 0);
            Index q_size = p.size() / 3;
            Orientations u(q_size);
#pragma omp parallel for
            for (int i = 0; i < q_size; i++) {
                //For a normalized quaternion q = w + xi + yj + zk
                Scalar1 w, x, y, z;
                x = p.coeff(i * 3);
                y = p.coeff(i * 3 + 1);
                z = p.coeff(i * 3 + 2);
                //Recover the real part
                w = sqrt(1.0 - x * x - y * y - z * z);
                w = isnan(w) ? 0.0 : w;
                Quaternion quat(w, x, y, z);
                u[i] = quat;
            }
            return u;
        }

        //从四元数的虚部创建向量。
        static Vector1 quat2pos(const Orientations quat) {
            Vector1 pos;
            pos.setZero(quat.size() * 3);
#pragma omp parallel for
            for (int i = 0; i < quat.size(); i++) {
                Quaternion q = quat[i].normalized();
                pos[i * 3] = q.x();
                pos[i * 3 + 1] = q.y();
                pos[i * 3 + 2] = q.z();
            }
            return pos;
        }

        //将平面向量表示的位置转换为3xN矩阵。
        static Positions vec2pos(const Vector1 p) {
            assert(p.size() % 3 == 0);
            Index q_size = p.size() / 3;
            Positions q;
            q.setZero(q_size, 3);

            for (Index i = 0; i < q_size; i++) {
                Vector3 vec;
                vec << p[i * 3], p[i * 3 + 1], p[i * 3 + 2];
                q.row(i) = vec;
            }

            return q;
        }

        //将3xN矩阵表示的位置转换为平面向量     
        static Vector1 pos2flatVec(const Positions pos) {
            Vector1 v;
            v.setZero(pos.rows() * 3);

            for (Index i = 0; i < pos.rows(); i++) {
                v[i * 3] = pos.coeff(i, 0);
                v[i * 3 + 1] = pos.coeff(i, 1);
                v[i * 3 + 2] = pos.coeff(i, 2);
            }

            return v;
        }

        //从给定的位置创建四元数。这些四元数表示从切线基向量到杆段的旋转
        Orientations quatFromPos(const Vector1 pos) {
            assert(pos.size() % 3 == 0);
            Index num_pos = pos.size() / 3;
            Orientations quat;
            quat.resize(num_pos - rod_indices.size());
            const Vector3 tangent = Vector3::UnitY();

            Index j = 0;
            for (Index i = 0; i < num_pos - 1; i++) {
                Index index = i * 3;

                Vector3 v_i, v_next;
                v_i << pos[index], pos[index + 1], pos[index + 2];
                v_next << pos[index + 3], pos[index + 4], pos[index + 5];

                if (i != 0 && std::find(rod_indices.begin(), rod_indices.end(), i + 1) != rod_indices.end()) i++;

                Quaternion q = Quaternion::FromTwoVectors(tangent, (v_next - v_i).normalized());
                quat[j++] = q;
            }

            return quat;
        }

        //给定一个四元数列表，执行逐元素的共轭操作。
        static Orientations conjugateQuat(const Orientations quat) {
            Orientations new_quat;
            new_quat.resize(quat.size());
            for (Index i = 0; i < quat.size(); i++) {
                new_quat[i] = quat[i].conjugate();
            }
            return new_quat;
        }

        //对两个向量进行逐列的叉积运算。
        static Vector1 vectorCross(const Vector1 a, const Vector1 b) {
            Vector1 cross;
            assert(a.size() == b.size() && a.size() % 3 == 0);
            cross.setZero(a.rows());
            Vector3 a_v, b_v, cross_v;

            for (Index i = 0; i < cross.size(); i += 3) {
                a_v << a[i], a[i + 1], a[i + 2];
                b_v << b[i], b[i + 1], b[i + 2];
                cross_v = a_v.cross(b_v);
                cross[i] = cross_v[0];
                cross[i + 1] = cross_v[1];
                cross[i + 2] = cross_v[2];
            }
            return cross;
        }

        //将位置的平面向量和四元数列表连接成一个大的平面向量
        static Vector1 posQuatConcat(const Vector1 pos, const Orientations u) {
            Index m = pos.size();
            Index dim = m + 4 * u.rows();
            Vector1 fp(dim);

            for (Index i = 0; i < m; i++) {
                fp(i) = pos[i];
            }

            for (Index i = 0; i < u.rows(); i++) {
                auto quat = u.coeff(i);
                fp(m + i * 4) = quat.w();
                fp(m + i * 4 + 1) = quat.x();
                fp(m + i * 4 + 2) = quat.y();
                fp(m + i * 4 + 3) = quat.z();
            }

            return fp;
        }

        //给定位置和方向的平面向量串联，将它们分开
        static void separatePosQuat(const Vector1* concat, Vector1& pos, Orientations& quat, Index separation) {
            Index size = concat->size();
            assert((size - separation) % 4 == 0);
            Index num_quat = (size - separation) / 4;

            pos.setZero(separation);

            quat.resize(num_quat);

            pos = concat->head(separation);

            for (Index i = 0; i < num_quat; i++) {
                Index index = separation + i * 4;
                quat.coeffRef(i) = Quaternion(concat->coeff(index),
                    concat->coeff(index + 1),
                    concat->coeff(index + 2),
                    concat->coeff(index + 3));
            }
        }

        //将给定矢量叠加N次
        static SparseMatrix stackDiagonal(const Vector1 x, Index times) {
            Index size = x.size();
            SparseMatrix concat;
            concat.resize(times * size, times * size);
            concat.setZero();
            for (Index i = 0; i < times; i++) {
                for (Index j = 0; j < size; j++) {
                    concat.coeffRef(i * size + j, i * size + j) = x.coeff(j);
                }
            }
            return concat;
        }

        static float clamp(float x, float min = 0.0f, float max = 1.0f) {
            x = x < min ? min : x;
            x = x > max ? max : x;
            return x;
        }

        //计算给定位置列表的重心
        Vector3 computeCenter(const Positions& positions) {
            Vector3 center = Vector3(0, 0, 0);

            for (Index i = 0; i < positions.rows(); i++) {
                center += positions.row(i);
            }

            return center / positions.rows();
        }

    };

}
