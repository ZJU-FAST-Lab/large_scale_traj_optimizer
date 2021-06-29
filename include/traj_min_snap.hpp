/*
    MIT License

    Copyright (c) 2020 Zhepei Wang (wangzhepei@live.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef TRAJ_MIN_SNAP_HPP
#define TRAJ_MIN_SNAP_HPP

#include "root_finder.hpp"

#include <Eigen/Eigen>

#include <iostream>
#include <cmath>
#include <vector>

namespace min_snap
{

    // Polynomial order and trajectory dimension are fixed here
    constexpr int TrajOrder = 7;
    constexpr int TrajDim = 3;

    // Type for piece boundary condition and coefficient matrix
    typedef Eigen::Matrix<double, TrajDim, TrajOrder + 1> BoundaryCond;
    typedef Eigen::Matrix<double, TrajDim, TrajOrder + 1> CoefficientMat;
    typedef Eigen::Matrix<double, TrajDim, TrajOrder> VelCoefficientMat;
    typedef Eigen::Matrix<double, TrajDim, TrajOrder - 1> AccCoefficientMat;

    // A single piece of a trajectory, which is indeed a polynomial
    class Piece
    {
    private:
        // A piece is totally determined by boundary condition and duration
        // boundCond = [p(0),v(0),a(0),j(0),p(T),v(T),a(T),j(T)]
        BoundaryCond boundCond;
        double duration;

        // The normalized coefficient is generated from boundCond and duration
        // These members CANNOT be accessed unless through normalizedCoeffMat()
        // p(t) = c7*t^7+c6*t^6 + ... + c1*t + c0
        // nCoeffMat = [c7*T^7,c6*T^6,c5*T^5,c4*T^4,c3*T^3,c2*T^2,c1*T,c0*1]
        bool synced;
        CoefficientMat nCoeffMat;
        inline const CoefficientMat &normalizedCoeffMat(void)
        {
            if (!synced)
            {
                double t1 = duration;
                double t2 = t1 * t1;
                double t3 = t2 * t1;

                nCoeffMat.col(0) = (boundCond.col(7) / 6.0 + boundCond.col(3) / 6.0) * t3 +
                                   (-2.0 * boundCond.col(6) + 2.0 * boundCond.col(2)) * t2 +
                                   (10.0 * boundCond.col(5) + 10.0 * boundCond.col(1)) * t1 +
                                   (-20.0 * boundCond.col(4) + 20.0 * boundCond.col(0));
                nCoeffMat.col(1) = (-0.5 * boundCond.col(7) - boundCond.col(3) / 1.5) * t3 +
                                   (6.5 * boundCond.col(6) - 7.5 * boundCond.col(2)) * t2 +
                                   (-34.0 * boundCond.col(5) - 36.0 * boundCond.col(1)) * t1 +
                                   (70.0 * boundCond.col(4) - 70.0 * boundCond.col(0));
                nCoeffMat.col(2) = (0.5 * boundCond.col(7) + boundCond.col(3)) * t3 +
                                   (-7.0 * boundCond.col(6) + 10.0 * boundCond.col(2)) * t2 +
                                   (39.0 * boundCond.col(5) + 45.0 * boundCond.col(1)) * t1 +
                                   (-84.0 * boundCond.col(4) + 84.0 * boundCond.col(0));
                nCoeffMat.col(3) = (-boundCond.col(7) / 6.0 - boundCond.col(3) / 1.5) * t3 +
                                   (2.5 * boundCond.col(6) - 5.0 * boundCond.col(2)) * t2 +
                                   (-15.0 * boundCond.col(5) - 20.0 * boundCond.col(1)) * t1 +
                                   (35.0 * boundCond.col(4) - 35.0 * boundCond.col(0));
                nCoeffMat.col(4) = boundCond.col(3) * t3 / 6.0;
                nCoeffMat.col(5) = boundCond.col(2) * t2 / 2.0;
                nCoeffMat.col(6) = boundCond.col(1) * t1;
                nCoeffMat.col(7) = boundCond.col(0);

                synced = true;
            }

            return nCoeffMat;
        }

    public:
        Piece() = default;

        // Constructor from boundary condition and duration
        Piece(BoundaryCond bdCond, double dur) : boundCond(bdCond), duration(dur), synced(false) {}

        inline int getDim() const
        {
            return TrajDim;
        }

        inline int getOrder() const
        {
            return TrajOrder;
        }

        inline double getDuration() const
        {
            return duration;
        }

        // Get the position at time t in this piece
        inline Eigen::Vector3d getPos(double t)
        {
            // Normalize the time
            t /= duration;
            Eigen::Vector3d pos(0.0, 0.0, 0.0);
            double tn = 1.0;
            for (int i = TrajOrder; i >= 0; i--)
            {
                pos += tn * normalizedCoeffMat().col(i);
                tn *= t;
            }
            // The pos is not affected by normalization
            return pos;
        }

        // Get the velocity at time t in this piece
        inline Eigen::Vector3d getVel(double t)
        {
            // Normalize the time
            t /= duration;
            Eigen::Vector3d vel(0.0, 0.0, 0.0);
            double tn = 1.0;
            int n = 1;
            for (int i = TrajOrder - 1; i >= 0; i--)
            {
                vel += n * tn * normalizedCoeffMat().col(i);
                tn *= t;
                n++;
            }
            // Recover the actual vel
            vel /= duration;
            return vel;
        }

        // Get the acceleration at time t in this piece
        inline Eigen::Vector3d getAcc(double t)
        {
            // Normalize the time
            t /= duration;
            Eigen::Vector3d acc(0.0, 0.0, 0.0);
            double tn = 1.0;
            int m = 1;
            int n = 2;
            for (int i = TrajOrder - 2; i >= 0; i--)
            {
                acc += m * n * tn * normalizedCoeffMat().col(i);
                tn *= t;
                m++;
                n++;
            }
            // Recover the actual acc
            acc /= duration * duration;
            return acc;
        }

        // Get the boundary condition of this piece
        inline const BoundaryCond &getBoundCond() const
        {
            return boundCond;
        }

        // Get the coefficient matrix of the piece
        // Default arg chooses the natural coefficients
        // If normalized version is needed, set the arg true
        inline CoefficientMat getCoeffMat(bool normalized = false)
        {
            CoefficientMat posCoeffsMat;
            double t = 1;
            for (int i = TrajOrder; i >= 0; i--)
            {
                posCoeffsMat.col(i) = normalizedCoeffMat().col(i) / t;
                t *= normalized ? 1.0 : duration;
            }
            return posCoeffsMat;
        }

        // Get the polynomial coefficients of velocity of this piece
        // Default arg chooses the natural coefficients
        // If normalized version is needed, set the arg true
        inline VelCoefficientMat getVelCoeffMat(bool normalized = false)
        {
            VelCoefficientMat velCoeffMat;
            int n = 1;
            double t = 1.0;
            t *= normalized ? 1.0 : duration;
            for (int i = TrajOrder - 1; i >= 0; i--)
            {
                velCoeffMat.col(i) = n * normalizedCoeffMat().col(i) / t;
                n++;
                t *= normalized ? 1.0 : duration;
            }
            return velCoeffMat;
        }

        // Get the polynomial coefficients of acceleration of this piece
        // Default arg chooses the natural coefficients
        // If normalized version is needed, set the arg true
        inline AccCoefficientMat getAccCoeffMat(bool normalized = false)
        {
            AccCoefficientMat accCoeffMat;
            int n = 2;
            int m = 1;
            double t = 1.0;
            t *= normalized ? 1.0 : duration * duration;
            for (int i = TrajOrder - 2; i >= 0; i--)
            {
                accCoeffMat.col(i) = n * m * normalizedCoeffMat().col(i) / t;
                n++;
                m++;
                t *= normalized ? 1.0 : duration;
            }
            return accCoeffMat;
        }

        // Get the max velocity rate of the piece
        inline double getMaxVelRate()
        {
            // Compute normalized squared vel norm polynomial coefficient matrix
            Eigen::MatrixXd nVelCoeffMat = getVelCoeffMat(true);
            Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                                    RootFinder::polySqr(nVelCoeffMat.row(1)) +
                                    RootFinder::polySqr(nVelCoeffMat.row(2));
            int N = coeff.size();
            int n = N - 1;
            for (int i = 0; i < N; i++)
            {
                coeff(i) *= n;
                n--;
            }
            if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON)
            {
                return getVel(0.0).norm();
            }
            else
            {
                // Search an open interval whose boundaries are not zeros
                double l = -0.0625;
                double r = 1.0625;
                while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON)
                {
                    l = 0.5 * l;
                }
                while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON)
                {
                    r = 0.5 * (r + 1.0);
                }
                // Find all stationaries
                std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                          FLT_EPSILON / duration);
                // Check boundary points and stationaries within duration
                candidates.insert(0.0);
                candidates.insert(1.0);
                double maxVelRateSqr = -INFINITY;
                double tempNormSqr;
                for (std::set<double>::const_iterator it = candidates.begin();
                     it != candidates.end();
                     it++)
                {
                    if (0.0 <= *it && 1.0 >= *it)
                    {
                        // Recover the actual time then get the vel squared norm
                        tempNormSqr = getVel((*it) * duration).squaredNorm();
                        maxVelRateSqr = maxVelRateSqr < tempNormSqr ? tempNormSqr : maxVelRateSqr;
                    }
                }
                return sqrt(maxVelRateSqr);
            }
        }

        // Get the max acceleration rate of the piece
        inline double getMaxAccRate()
        {
            // Compute normalized squared acc norm polynomial coefficient matrix
            Eigen::MatrixXd nAccCoeffMat = getAccCoeffMat(true);
            Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                                    RootFinder::polySqr(nAccCoeffMat.row(1)) +
                                    RootFinder::polySqr(nAccCoeffMat.row(2));
            int N = coeff.size();
            int n = N - 1;
            for (int i = 0; i < N; i++)
            {
                coeff(i) *= n;
                n--;
            }
            if (coeff.head(N - 1).squaredNorm() < DBL_EPSILON)
            {
                return getAcc(0.0).norm();
            }
            else
            {
                // Search an open interval whose boundaries are not zeros
                double l = -0.0625;
                double r = 1.0625;
                while (fabs(RootFinder::polyVal(coeff.head(N - 1), l)) < DBL_EPSILON)
                {
                    l = 0.5 * l;
                }
                while (fabs(RootFinder::polyVal(coeff.head(N - 1), r)) < DBL_EPSILON)
                {
                    r = 0.5 * (r + 1.0);
                }
                // Find all stationaries
                std::set<double> candidates = RootFinder::solvePolynomial(coeff.head(N - 1), l, r,
                                                                          FLT_EPSILON / duration);
                // Check boundary points and stationaries within duration
                candidates.insert(0.0);
                candidates.insert(1.0);
                double maxAccRateSqr = -INFINITY;
                double tempNormSqr;
                for (std::set<double>::const_iterator it = candidates.begin();
                     it != candidates.end();
                     it++)
                {
                    if (0.0 <= *it && 1.0 >= *it)
                    {
                        // Recover the actual time then get the acc squared norm
                        tempNormSqr = getAcc((*it) * duration).squaredNorm();
                        maxAccRateSqr = maxAccRateSqr < tempNormSqr ? tempNormSqr : maxAccRateSqr;
                    }
                }
                return sqrt(maxAccRateSqr);
            }
        }

        // Check whether velocity rate of the piece is always less than maxVelRate
        inline bool checkMaxVelRate(double maxVelRate)
        {
            double sqrMaxVelRate = maxVelRate * maxVelRate;
            if (getVel(0.0).squaredNorm() >= sqrMaxVelRate ||
                getVel(duration).squaredNorm() >= sqrMaxVelRate)
            {
                return false;
            }
            else
            {
                Eigen::MatrixXd nVelCoeffMat = getVelCoeffMat(true);
                Eigen::VectorXd coeff = RootFinder::polySqr(nVelCoeffMat.row(0)) +
                                        RootFinder::polySqr(nVelCoeffMat.row(1)) +
                                        RootFinder::polySqr(nVelCoeffMat.row(2));
                // Convert the actual squared maxVelRate to a normalized one
                double t2 = duration * duration;
                coeff.tail<1>()(0) -= sqrMaxVelRate * t2;
                // Directly check the root existence in the normalized interval
                return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
            }
        }

        // Check whether accleration rate of the piece is always less than maxAccRate
        inline bool checkMaxAccRate(double maxAccRate)
        {
            double sqrMaxAccRate = maxAccRate * maxAccRate;
            if (getAcc(0.0).squaredNorm() >= sqrMaxAccRate ||
                getAcc(duration).squaredNorm() >= sqrMaxAccRate)
            {
                return false;
            }
            else
            {
                Eigen::MatrixXd nAccCoeffMat = getAccCoeffMat(true);
                Eigen::VectorXd coeff = RootFinder::polySqr(nAccCoeffMat.row(0)) +
                                        RootFinder::polySqr(nAccCoeffMat.row(1)) +
                                        RootFinder::polySqr(nAccCoeffMat.row(2));
                // Convert the actual squared maxAccRate to a normalized one
                double t2 = duration * duration;
                double t4 = t2 * t2;
                coeff.tail<1>()(0) -= sqrMaxAccRate * t4;
                // Directly check the root existence in the normalized interval
                return RootFinder::countRoots(coeff, 0.0, 1.0) == 0;
            }
        }
    };

    // A whole trajectory which contains multiple pieces
    class Trajectory
    {
    private:
        typedef std::vector<Piece> Pieces;
        Pieces pieces;

    public:
        Trajectory() = default;

        // Constructor from boundary conditions and durations
        Trajectory(const std::vector<BoundaryCond> &bdConds,
                   const std::vector<double> &durs)
        {
            int N = std::min(durs.size(), bdConds.size());
            pieces.reserve(N);
            for (int i = 0; i < N; i++)
            {
                pieces.emplace_back(bdConds[i], durs[i]);
            }
        }

        inline int getPieceNum() const
        {
            return pieces.size();
        }

        // Get durations vector of all pieces
        inline Eigen::VectorXd getDurations() const
        {
            int N = getPieceNum();
            Eigen::VectorXd durations(N);
            for (int i = 0; i < N; i++)
            {
                durations(i) = pieces[i].getDuration();
            }
            return durations;
        }

        // Get total duration of the trajectory
        inline double getTotalDuration() const
        {
            int N = getPieceNum();
            double totalDuration = 0.0;
            for (int i = 0; i < N; i++)
            {
                totalDuration += pieces[i].getDuration();
            }
            return totalDuration;
        }

        inline Eigen::MatrixXd getPositions() const
        {
            int N = getPieceNum();
            Eigen::MatrixXd positions(3, N + 1);
            for (int i = 0; i < N; i++)
            {
                positions.col(i) = pieces[i].getBoundCond().col(0);
            }
            positions.col(N) = pieces[N - 1].getBoundCond().col((TrajOrder + 1) / 2);
            return positions;
        }

        // Reload the operator[] to access the i-th piece
        inline const Piece &operator[](int i) const
        {
            return pieces[i];
        }

        inline Piece &operator[](int i)
        {
            return pieces[i];
        }

        inline void clear(void)
        {
            pieces.clear();
            return;
        }

        inline Pieces::const_iterator begin() const
        {
            return pieces.begin();
        }

        inline Pieces::const_iterator end() const
        {
            return pieces.end();
        }

        inline Pieces::iterator begin()
        {
            return pieces.begin();
        }

        inline Pieces::iterator end()
        {
            return pieces.end();
        }

        inline void reserve(const int &n)
        {
            pieces.reserve(n);
            return;
        }

        // Put another piece at the tail of this trajectory
        inline void emplace_back(const Piece &piece)
        {
            pieces.emplace_back(piece);
            return;
        }

        inline void emplace_back(const BoundaryCond &bdCond, const double &dur)
        {
            pieces.emplace_back(bdCond, dur);
            return;
        }

        // Append another Trajectory at the tail of this trajectory
        inline void append(const Trajectory &traj)
        {
            pieces.insert(pieces.end(), traj.begin(), traj.end());
            return;
        }

        // Find the piece at which the time t is located
        // The index is returned and the offset in t is removed
        inline int locatePieceIdx(double &t) const
        {
            int N = getPieceNum();
            int idx;
            double dur;
            for (idx = 0;
                 idx < N &&
                 t > (dur = pieces[idx].getDuration());
                 idx++)
            {
                t -= dur;
            }
            if (idx == N)
            {
                idx--;
                t += pieces[idx].getDuration();
            }
            return idx;
        }

        // Get the position at time t of the trajectory
        inline Eigen::Vector3d getPos(double t)
        {
            int pieceIdx = locatePieceIdx(t);
            return pieces[pieceIdx].getPos(t);
        }

        // Get the velocity at time t of the trajectory
        inline Eigen::Vector3d getVel(double t)
        {
            int pieceIdx = locatePieceIdx(t);
            return pieces[pieceIdx].getVel(t);
        }

        // Get the acceleration at time t of the trajectory
        inline Eigen::Vector3d getAcc(double t)
        {
            int pieceIdx = locatePieceIdx(t);
            return pieces[pieceIdx].getAcc(t);
        }

        // Get the position at the juncIdx-th waypoint
        inline Eigen::Vector3d getJuncPos(int juncIdx) const
        {
            if (juncIdx != getPieceNum())
            {
                return pieces[juncIdx].getBoundCond().col(0);
            }
            else
            {
                return pieces[juncIdx - 1].getBoundCond().col((TrajOrder + 1) / 2);
            }
        }

        // Get the velocity at the juncIdx-th waypoint
        inline Eigen::Vector3d getJuncVel(int juncIdx) const
        {
            if (juncIdx != getPieceNum())
            {
                return pieces[juncIdx].getBoundCond().col(1);
            }
            else
            {
                return pieces[juncIdx - 1].getBoundCond().col((TrajOrder + 1) / 2 + 1);
            }
        }

        // Get the acceleration at the juncIdx-th waypoint
        inline Eigen::Vector3d getJuncAcc(int juncIdx) const
        {
            if (juncIdx != getPieceNum())
            {
                return pieces[juncIdx].getBoundCond().col(2);
            }
            else
            {
                return pieces[juncIdx - 1].getBoundCond().col((TrajOrder + 1) / 2 + 2);
            }
        }

        // Get the max velocity rate of the trajectory
        inline double getMaxVelRate()
        {
            int N = getPieceNum();
            double maxVelRate = -INFINITY;
            double tempNorm;
            for (int i = 0; i < N; i++)
            {
                tempNorm = pieces[i].getMaxVelRate();
                maxVelRate = maxVelRate < tempNorm ? tempNorm : maxVelRate;
            }
            return maxVelRate;
        }

        // Get the max acceleration rate of the trajectory
        inline double getMaxAccRate()
        {
            int N = getPieceNum();
            double maxAccRate = -INFINITY;
            double tempNorm;
            for (int i = 0; i < N; i++)
            {
                tempNorm = pieces[i].getMaxAccRate();
                maxAccRate = maxAccRate < tempNorm ? tempNorm : maxAccRate;
            }
            return maxAccRate;
        }

        // Check whether the velocity rate of this trajectory exceeds the threshold
        inline bool checkMaxVelRate(double maxVelRate)
        {
            int N = getPieceNum();
            bool feasible = true;
            for (int i = 0; i < N && feasible; i++)
            {
                feasible = feasible && pieces[i].checkMaxVelRate(maxVelRate);
            }
            return feasible;
        }

        // Check whether the acceleration rate of this trajectory exceeds the threshold
        inline bool checkMaxAccRate(double maxAccRate)
        {
            int N = getPieceNum();
            bool feasible = true;
            for (int i = 0; i < N && feasible; i++)
            {
                feasible = feasible && pieces[i].checkMaxAccRate(maxAccRate);
            }
            return feasible;
        }
    };

    // The banded system class is used for solving
    // banded linear system Ax=b efficiently.
    // A is an N*N band matrix with lower band width lowerBw
    // and upper band width upperBw.
    // Banded LU factorization has O(N) time complexity.
    class BandedSystem
    {
    public:
        // The size of A, as well as the lower/upper
        // banded width p/q are needed
        inline void create(const int &n, const int &p, const int &q)
        {
            // In case of re-creating before destroying
            destroy();
            N = n;
            lowerBw = p;
            upperBw = q;
            int rows = lowerBw + upperBw + 1;
            int actualSize = N * rows;
            ptrData = new double[actualSize];
            std::fill_n(ptrData, actualSize, 0.0);
            offset = new double *[rows];
            double *ptrRow = ptrData;
            for (int i = 0; i < rows; i++)
            {
                offset[i] = ptrRow;
                ptrRow += N;
            }
            return;
        }

        inline void destroy()
        {
            if (ptrData != nullptr)
            {
                delete[] ptrData;
                ptrData = nullptr;
            }
            if (offset != nullptr)
            {
                delete[] offset;
                offset = nullptr;
            }
            return;
        }

    private:
        int N;
        int lowerBw;
        int upperBw;
        // Compulsory nullptr initialization here
        double *ptrData = nullptr;
        double **offset = nullptr;

    public:
        // Reset the matrix to zero
        inline void reset(void)
        {
            std::fill_n(ptrData, N * (lowerBw + upperBw + 1), 0.0);
            return;
        }

        // The band matrix is stored as suggested in "Matrix Computation"
        inline const double &operator()(const int &i, const int &j) const
        {
            return offset[i - j + upperBw][j];
        }

        inline double &operator()(const int &i, const int &j)
        {
            return offset[i - j + upperBw][j];
        }

        // This function conducts banded LU factorization in place
        // Note that the matrix "A" MUST NOT HAVE ZERO PIVOTS !!!
        inline void factorizeLU()
        {
            int iM, jM;
            for (int k = 0; k <= N - 2; k++)
            {
                iM = std::min(k + lowerBw, N - 1);
                for (int i = k + 1; i <= iM; i++)
                {
                    operator()(i, k) /= operator()(k, k);
                }
                jM = std::min(k + upperBw, N - 1);
                for (int j = k + 1; j <= jM; j++)
                {
                    for (int i = k + 1; i <= iM; i++)
                    {
                        operator()(i, j) -= operator()(i, k) * operator()(k, j);
                    }
                }
            }
            return;
        }

        // This function solves Ax=b, then stores x in b
        // The input b is required to be N*m, i.e.,
        // m vectors to be solved.
        template <typename EIGENMAT>
        inline void solve(EIGENMAT &b) const
        {
            int iM;
            for (int j = 0; j <= N - 1; j++)
            {
                iM = std::min(j + lowerBw, N - 1);
                for (int i = j + 1; i <= iM; i++)
                {
                    b.row(i) -= operator()(i, j) * b.row(j);
                }
            }
            for (int j = N - 1; j >= 0; j--)
            {
                b.row(j) /= operator()(j, j);
                iM = std::max(0, j - upperBw);
                for (int i = iM; i <= j - 1; i++)
                {
                    b.row(i) -= operator()(i, j) * b.row(j);
                }
            }
            return;
        }
    };

    class SnapOpt
    {
    public:
        SnapOpt() = default;
        ~SnapOpt() { A.destroy(); }

    private:
        int N;
        Eigen::Matrix3Xd Ps;
        Eigen::Matrix3Xd VAJs;
        Eigen::VectorXd T;
        BandedSystem A;
        Eigen::MatrixX3d b;

        // Temp variables
        Eigen::VectorXd t2;
        Eigen::VectorXd t3;
        Eigen::VectorXd t4;
        Eigen::VectorXd t5;
        Eigen::VectorXd t6;
        Eigen::VectorXd t7;
        Eigen::VectorXd cv00, cv01, cv02;
        Eigen::VectorXd cv10, cv11, cv12;
        Eigen::VectorXd cv20, cv21, cv22;
        Eigen::VectorXd cv30, cv31, cv32;
        Eigen::VectorXd ca00, ca01, ca02;
        Eigen::VectorXd ca10, ca11, ca12;
        Eigen::VectorXd ca20, ca21, ca22;
        Eigen::VectorXd ca30, ca31, ca32;
        Eigen::VectorXd cj00, cj01, cj02;
        Eigen::VectorXd cj10, cj11, cj12;
        Eigen::VectorXd cj20, cj21, cj22;
        Eigen::VectorXd cj30, cj31, cj32;

    private:
        inline double evalPieceObjective(const Eigen::Array3d &iP,
                                         const Eigen::Array3d &iV,
                                         const Eigen::Array3d &iA,
                                         const Eigen::Array3d &iJ,
                                         const Eigen::Array3d &fP,
                                         const Eigen::Array3d &fV,
                                         const Eigen::Array3d &fA,
                                         const Eigen::Array3d &fJ,
                                         const double &duration) const
        {
            Eigen::VectorXd coeffsSnpObjective(7);
            coeffsSnpObjective(0) = 8.0 * (2.0 * iJ.square() + iJ * fJ + 2.0 * fJ.square()).sum();
            coeffsSnpObjective(1) = 120.0 * (iA * (2.0 * iJ + fJ) - fA * (iJ + 2.0 * fJ)).sum();
            coeffsSnpObjective(2) = 240.0 * (5.0 * iA.square() - 7.0 * iA * fA + 5.0 * fA.square() + 4.0 * iJ * iV + 3.0 * fJ * iV + 3.0 * iJ * fV + 4.0 * fJ * fV).sum();
            coeffsSnpObjective(3) = 240.0 * (45.0 * iA * iV - 39.0 * fA * iV + 39.0 * iA * fV - 45.0 * fA * fV + 7.0 * (iJ + fJ) * (iP - fP)).sum();
            coeffsSnpObjective(4) = 2880.0 * (9.0 * iV.square() + 17.0 * iV * fV + 9.0 * fV.square() + 7.0 * (iA - fA) * (iP - fP)).sum();
            coeffsSnpObjective(5) = 100800.0 * ((iV + fV) * (iP - fP)).sum();
            coeffsSnpObjective(6) = 100800.0 * (iP - fP).square().sum();

            double t2 = duration * duration;
            double t6 = t2 * t2 * t2;
            double t7 = t6 * duration;

            return RootFinder::polyVal(coeffsSnpObjective, duration) / t7;
        }

    public:
        inline void reset(const Eigen::Matrix<double, 3, 4> &headState,
                          const Eigen::Matrix<double, 3, 4> &tailState,
                          const int &pieceNum)
        {
            N = pieceNum;
            Ps.resize(3, N + 1);
            VAJs.resize(3, 3 * N + 3);
            Ps.leftCols<1>() = headState.leftCols<1>();
            Ps.rightCols<1>() = tailState.leftCols<1>();
            VAJs.leftCols<3>() = headState.rightCols<3>();
            VAJs.rightCols<3>() = tailState.rightCols<3>();
            T.resize(N);
            A.create(3 * N - 3, 5, 5);
            b.resize(3 * N - 3, 3);

            cv00.resize(N - 1);
            cv01.resize(N - 1);
            cv02.resize(N - 1);
            cv10.resize(N - 1);
            cv11.resize(N - 1);
            cv12.resize(N - 1);
            cv20.resize(N - 1);
            cv21.resize(N - 1);
            cv22.resize(N - 1);
            cv30.resize(N - 1);
            cv31.resize(N - 1);
            cv32.resize(N - 1);
            ca00.resize(N - 1);
            ca01.resize(N - 1);
            ca02.resize(N - 1);
            ca10.resize(N - 1);
            ca11.resize(N - 1);
            ca12.resize(N - 1);
            ca20.resize(N - 1);
            ca21.resize(N - 1);
            ca22.resize(N - 1);
            ca30.resize(N - 1);
            ca31.resize(N - 1);
            ca32.resize(N - 1);
            cj00.resize(N - 1);
            cj01.resize(N - 1);
            cj02.resize(N - 1);
            cj10.resize(N - 1);
            cj11.resize(N - 1);
            cj12.resize(N - 1);
            cj20.resize(N - 1);
            cj21.resize(N - 1);
            cj22.resize(N - 1);
            cj30.resize(N - 1);
            cj31.resize(N - 1);
            cj32.resize(N - 1);

            return;
        }

        inline void generate(const Eigen::Matrix3Xd &inPs,
                             const Eigen::VectorXd &ts)
        {
            T = ts;
            const Eigen::VectorXd &t1 = T;
            t2 = t1.cwiseProduct(t1);
            t3 = t2.cwiseProduct(t1);
            t4 = t2.cwiseProduct(t2);
            t5 = t4.cwiseProduct(t1);
            t6 = t3.cwiseProduct(t3);
            t7 = t6.cwiseProduct(t1);

            // Computed nonzero entries in A and b for linear system Ax=b to be solved
            for (int i = 0; i < N - 1; i++)
            {
                cv00(i) = 100800.0 / t6(i);
                cv01(i) = 100800.0 * (1.0 / t6(i + 1) - 1.0 / t6(i));
                cv02(i) = -100800.0 / t6(i + 1);
                cv10(i) = 48960.0 / t5(i);
                cv11(i) = 51840.0 * (1.0 / t5(i + 1) + 1.0 / t5(i));
                cv12(i) = 48960.0 / t5(i + 1);
                cv20(i) = 9360.0 / t4(i);
                cv21(i) = 10800.0 * (1.0 / t4(i + 1) - 1.0 / t4(i));
                cv22(i) = -9360.0 / t4(i + 1);
                cv30(i) = 720.0 / t3(i);
                cv31(i) = 960.0 * (1.0 / t3(i + 1) + 1.0 / t3(i));
                cv32(i) = 720.0 / t3(i + 1);

                ca00(i) = -20160.0 / t5(i);
                ca01(i) = 20160.0 * (1.0 / t5(i + 1) + 1.0 / t5(i));
                ca02(i) = -20160.0 / t5(i + 1);
                ca10(i) = -9360.0 / t4(i);
                ca11(i) = 10800.0 * (1.0 / t4(i + 1) - 1.0 / t4(i));
                ca12(i) = 9360.0 / t4(i + 1);
                ca20(i) = -1680.0 / t3(i);
                ca21(i) = 2400.0 * (1.0 / t3(i + 1) + 1.0 / t3(i));
                ca22(i) = -1680.0 / t3(i + 1);
                ca30(i) = -120.0 / t2(i);
                ca31(i) = 240.0 * (1.0 / t2(i + 1) - 1.0 / t2(i));
                ca32(i) = 120.0 / t2(i + 1);

                cj00(i) = 1680.0 / t4(i);
                cj01(i) = 1680.0 * (1.0 / t4(i + 1) - 1.0 / t4(i));
                cj02(i) = -1680.0 / t4(i + 1);
                cj10(i) = 720.0 / t3(i);
                cj11(i) = 960.0 * (1.0 / t3(i + 1) + 1.0 / t3(i));
                cj12(i) = 720.0 / t3(i + 1);
                cj20(i) = 120.0 / t2(i);
                cj21(i) = 240.0 * (1.0 / t2(i + 1) - 1.0 / t2(i));
                cj22(i) = -120.0 / t2(i + 1);
                cj30(i) = 8.0 / t1(i);
                cj31(i) = 32.0 * (1.0 / t1(i + 1) + 1.0 / t1(i));
                cj32(i) = 8.0 / t1(i + 1);
            }

            Ps.block(0, 1, 3, N - 1) = inPs;

            if (N == 2)
            {
                Eigen::Matrix3d A33, invA, bl;

                A33(0, 0) = cv11(0);
                A33(0, 1) = cv21(0);
                A33(0, 2) = cv31(0);
                A33(1, 0) = ca11(0);
                A33(1, 1) = ca21(0);
                A33(1, 2) = ca31(0);
                A33(2, 0) = cj11(0);
                A33(2, 1) = cj21(0);
                A33(2, 2) = cj31(0);

                invA = A33.inverse();
                bl.row(0) = (-cv00(0) * Ps.col(0) - cv01(0) * Ps.col(1) - cv02(0) * Ps.col(2) - cv10(0) * VAJs.col(0) - cv20(0) * VAJs.col(1) - cv30(0) * VAJs.col(2) - cv12(0) * VAJs.col(6) - cv22(0) * VAJs.col(7) - cv32(0) * VAJs.col(8)).transpose();
                bl.row(1) = (-ca00(0) * Ps.col(0) - ca01(0) * Ps.col(1) - ca02(0) * Ps.col(2) - ca10(0) * VAJs.col(0) - ca20(0) * VAJs.col(1) - ca30(0) * VAJs.col(2) - ca12(0) * VAJs.col(6) - ca22(0) * VAJs.col(7) - ca32(0) * VAJs.col(8)).transpose();
                bl.row(2) = (-cj00(0) * Ps.col(0) - cj01(0) * Ps.col(1) - cj02(0) * Ps.col(2) - cj10(0) * VAJs.col(0) - cj20(0) * VAJs.col(1) - cj30(0) * VAJs.col(2) - cj12(0) * VAJs.col(6) - cj22(0) * VAJs.col(7) - cj32(0) * VAJs.col(8)).transpose();

                VAJs.block(0, 3, 3, 3) = (invA * bl).transpose();
            }
            else
            {
                A.reset();
                b.setZero();

                A(0, 0) = cv11(0);
                A(0, 1) = cv21(0);
                A(0, 2) = cv31(0);
                A(0, 3) = cv12(0);
                A(0, 4) = cv22(0);
                A(0, 5) = cv32(0);
                A(1, 0) = ca11(0);
                A(1, 1) = ca21(0);
                A(1, 2) = ca31(0);
                A(1, 3) = ca12(0);
                A(1, 4) = ca22(0);
                A(1, 5) = ca32(0);
                A(2, 0) = cj11(0);
                A(2, 1) = cj21(0);
                A(2, 2) = cj31(0);
                A(2, 3) = cj12(0);
                A(2, 4) = cj22(0);
                A(2, 5) = cj32(0);
                A(3 * N - 6, 3 * N - 9) = cv10(N - 2);
                A(3 * N - 6, 3 * N - 8) = cv20(N - 2);
                A(3 * N - 6, 3 * N - 7) = cv30(N - 2);
                A(3 * N - 6, 3 * N - 6) = cv11(N - 2);
                A(3 * N - 6, 3 * N - 5) = cv21(N - 2);
                A(3 * N - 6, 3 * N - 4) = cv31(N - 2);
                A(3 * N - 5, 3 * N - 9) = ca10(N - 2);
                A(3 * N - 5, 3 * N - 8) = ca20(N - 2);
                A(3 * N - 5, 3 * N - 7) = ca30(N - 2);
                A(3 * N - 5, 3 * N - 6) = ca11(N - 2);
                A(3 * N - 5, 3 * N - 5) = ca21(N - 2);
                A(3 * N - 5, 3 * N - 4) = ca31(N - 2);
                A(3 * N - 4, 3 * N - 9) = cj10(N - 2);
                A(3 * N - 4, 3 * N - 8) = cj20(N - 2);
                A(3 * N - 4, 3 * N - 7) = cj30(N - 2);
                A(3 * N - 4, 3 * N - 6) = cj11(N - 2);
                A(3 * N - 4, 3 * N - 5) = cj21(N - 2);
                A(3 * N - 4, 3 * N - 4) = cj31(N - 2);

                b.row(0) = (-cv00(0) * Ps.col(0) - cv01(0) * Ps.col(1) - cv02(0) * Ps.col(2) - cv10(0) * VAJs.col(0) - cv20(0) * VAJs.col(1) - cv30(0) * VAJs.col(2)).transpose();
                b.row(1) = (-ca00(0) * Ps.col(0) - ca01(0) * Ps.col(1) - ca02(0) * Ps.col(2) - ca10(0) * VAJs.col(0) - ca20(0) * VAJs.col(1) - ca30(0) * VAJs.col(2)).transpose();
                b.row(2) = (-cj00(0) * Ps.col(0) - cj01(0) * Ps.col(1) - cj02(0) * Ps.col(2) - cj10(0) * VAJs.col(0) - cj20(0) * VAJs.col(1) - cj30(0) * VAJs.col(2)).transpose();
                b.row(3 * N - 6) = (-cv00(N - 2) * Ps.col(N - 2) - cv01(N - 2) * Ps.col(N - 1) - cv02(N - 2) * Ps.col(N) - cv12(N - 2) * VAJs.col(3 * N) - cv22(N - 2) * VAJs.col(3 * N + 1) - cv32(N - 2) * VAJs.col(3 * N + 2)).transpose();
                b.row(3 * N - 5) = (-ca00(N - 2) * Ps.col(N - 2) - ca01(N - 2) * Ps.col(N - 1) - ca02(N - 2) * Ps.col(N) - ca12(N - 2) * VAJs.col(3 * N) - ca22(N - 2) * VAJs.col(3 * N + 1) - ca32(N - 2) * VAJs.col(3 * N + 2)).transpose();
                b.row(3 * N - 4) = (-cj00(N - 2) * Ps.col(N - 2) - cj01(N - 2) * Ps.col(N - 1) - cj02(N - 2) * Ps.col(N) - cj12(N - 2) * VAJs.col(3 * N) - cj22(N - 2) * VAJs.col(3 * N + 1) - cj32(N - 2) * VAJs.col(3 * N + 2)).transpose();

                for (int i = 1; i < N - 2; i++)
                {
                    A(i * 3, i * 3 - 3) = cv10(i);
                    A(i * 3, i * 3 - 2) = cv20(i);
                    A(i * 3, i * 3 - 1) = cv30(i);
                    A(i * 3, i * 3) = cv11(i);
                    A(i * 3, i * 3 + 1) = cv21(i);
                    A(i * 3, i * 3 + 2) = cv31(i);
                    A(i * 3, i * 3 + 3) = cv12(i);
                    A(i * 3, i * 3 + 4) = cv22(i);
                    A(i * 3, i * 3 + 5) = cv32(i);
                    A(i * 3 + 1, i * 3 - 3) = ca10(i);
                    A(i * 3 + 1, i * 3 - 2) = ca20(i);
                    A(i * 3 + 1, i * 3 - 1) = ca30(i);
                    A(i * 3 + 1, i * 3) = ca11(i);
                    A(i * 3 + 1, i * 3 + 1) = ca21(i);
                    A(i * 3 + 1, i * 3 + 2) = ca31(i);
                    A(i * 3 + 1, i * 3 + 3) = ca12(i);
                    A(i * 3 + 1, i * 3 + 4) = ca22(i);
                    A(i * 3 + 1, i * 3 + 5) = ca32(i);
                    A(i * 3 + 2, i * 3 - 3) = cj10(i);
                    A(i * 3 + 2, i * 3 - 2) = cj20(i);
                    A(i * 3 + 2, i * 3 - 1) = cj30(i);
                    A(i * 3 + 2, i * 3) = cj11(i);
                    A(i * 3 + 2, i * 3 + 1) = cj21(i);
                    A(i * 3 + 2, i * 3 + 2) = cj31(i);
                    A(i * 3 + 2, i * 3 + 3) = cj12(i);
                    A(i * 3 + 2, i * 3 + 4) = cj22(i);
                    A(i * 3 + 2, i * 3 + 5) = cj32(i);

                    b.row(i * 3) = (-cv00(i) * Ps.col(i) - cv01(i) * Ps.col(i + 1) - cv02(i) * Ps.col(i + 2)).transpose();
                    b.row(i * 3 + 1) = (-ca00(i) * Ps.col(i) - ca01(i) * Ps.col(i + 1) - ca02(i) * Ps.col(i + 2)).transpose();
                    b.row(i * 3 + 2) = (-cj00(i) * Ps.col(i) - cj01(i) * Ps.col(i + 1) - cj02(i) * Ps.col(i + 2)).transpose();
                }

                // Solve Ax=b using banded LU factorization
                A.factorizeLU();
                // The solution is computed in place.
                A.solve(b);

                VAJs.block(0, 3, 3, 3 * N - 3) = b.transpose();
            }

            return;
        }

        inline double getObjective() const
        {
            double objective = 0.0;

            for (int i = 0; i < N; i++)
            {
                objective += evalPieceObjective(Ps.col(i), VAJs.col(3 * i), VAJs.col(3 * i + 1), VAJs.col(3 * i + 2),
                                                Ps.col(i + 1), VAJs.col(3 * i + 3), VAJs.col(3 * i + 4), VAJs.col(3 * i + 5),
                                                T(i));
            }

            return objective;
        }

        inline Eigen::VectorXd getGradT(void) const
        {
            Eigen::VectorXd grad(N);

            double tempT, tempT8;
            Eigen::Array3d iP, iV, iA, iJ, fP, fV, fA, fJ;
            Eigen::VectorXd coeffsGradT(7);
            for (int i = 0; i < N; i++)
            {
                // Get the information of the piece
                tempT = T(i);
                tempT8 = tempT * tempT;
                tempT8 = tempT8 * tempT8;
                tempT8 = tempT8 * tempT8;

                // Calculate the numerator of dJi(T)/dT without time regularization
                iP = Ps.col(i);
                iV = VAJs.col(3 * i);
                iA = VAJs.col(3 * i + 1);
                iJ = VAJs.col(3 * i + 2);
                fP = Ps.col(i + 1);
                fV = VAJs.col(3 * i + 3);
                fA = VAJs.col(3 * i + 4);
                fJ = VAJs.col(3 * i + 5);

                coeffsGradT(0) = -8.0 * (2.0 * iJ.square() + iJ * fJ + 2.0 * fJ.square()).sum();
                coeffsGradT(1) = -240.0 * (iA * (2.0 * iJ + fJ) - fA * (iJ + 2.0 * fJ)).sum();
                coeffsGradT(2) = -720.0 * (5.0 * iA.square() - 7.0 * iA * fA + 5.0 * fA.square() + 4.0 * iJ * iV + 3.0 * fJ * iV + 3.0 * iJ * fV + 4.0 * fJ * fV).sum();
                coeffsGradT(3) = -960.0 * (45.0 * iA * iV - 39.0 * fA * iV + 39.0 * iA * fV - 45.0 * fA * fV + 7.0 * (iJ + fJ) * (iP - fP)).sum();
                coeffsGradT(4) = -14400.0 * (9.0 * iV.square() + 17.0 * iV * fV + 9.0 * fV.square() + 7.0 * (iA - fA) * (iP - fP)).sum();
                coeffsGradT(5) = 604800.0 * ((iV + fV) * (iP - fP)).sum();
                coeffsGradT(6) = -705600.0 * (iP - fP).square().sum();

                // Calculate the gradient
                grad(i) = RootFinder::polyVal(coeffsGradT, tempT) / tempT8;
            }

            return grad;
        }

        inline Eigen::Matrix3Xd getGradInnerP(void) const
        {
            Eigen::Matrix3Xd grad(3, N - 1);

            for (int i = 1; i < N; i++)
            {
                grad.col(i - 1) = 1680.0 * (-VAJs.col(3 * i - 1) / t4(i - 1) - VAJs.col(3 * i + 2) * (1.0 / t4(i - 1) - 1.0 / t4(i)) + VAJs.col(3 * i + 5) / t4(i)) +
                                  20160.0 * (-VAJs.col(3 * i - 2) / t5(i - 1) + VAJs.col(3 * i + 1) * (1.0 / t5(i - 1) + 1.0 / t5(i)) - VAJs.col(3 * i + 4) / t5(i)) +
                                  100800.0 * (-VAJs.col(3 * i - 3) / t6(i - 1) - VAJs.col(3 * i) * (1.0 / t6(i - 1) - 1.0 / t6(i)) + VAJs.col(3 * i + 3) / t6(i)) +
                                  201600.0 * (-Ps.col(i - 1) / t7(i - 1) + Ps.col(i) * (1.0 / t7(i - 1) + 1.0 / t7(i)) - Ps.col(i + 1) / t7(i));
            }

            return grad;
        }

        inline void getTraj(Trajectory &traj) const
        {
            traj.clear();
            traj.reserve(N);
            BoundaryCond boundCond;
            for (int i = 0; i < N; i++)
            {
                boundCond.col(0) = Ps.col(i);
                boundCond.col(1) = VAJs.col(3 * i);
                boundCond.col(2) = VAJs.col(3 * i + 1);
                boundCond.col(3) = VAJs.col(3 * i + 2);
                boundCond.col(4) = Ps.col(i + 1);
                boundCond.col(5) = VAJs.col(3 * i + 3);
                boundCond.col(6) = VAJs.col(3 * i + 4);
                boundCond.col(7) = VAJs.col(3 * i + 5);

                traj.emplace_back(boundCond, T(i));
            }
            return;
        }
    };

} // namespace min_snap

#endif
