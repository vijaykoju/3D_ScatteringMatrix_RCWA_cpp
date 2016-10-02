#ifndef RCWA_UTILS_H
#define RCWA_UTILS_H

#define pi_const 3.14159265359
#define c_const 299792458.0
#define nanometers 1e-9
#define degrees (pi_const/180.0)

#include <iostream>
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <complex>
#include <cmath>

using namespace Eigen;
using namespace std;

typedef struct g{
	int Nx;
	int Ny;
} Grid;

typedef struct pq{
	int P;
	int Q;
}PQ;

typedef struct sm{
	MatrixXcd s11;
	MatrixXcd s12;
	MatrixXcd s21;
	MatrixXcd s22;
} scatMat;

////////////////////////////////////////////////////////////////////////////////////////////////////
scatMat scattering_matrix_ith(const MatrixXcd &ERC, const MatrixXcd &URC, const MatrixXcd &Kx, const MatrixXcd &Ky, double k0, double L, const MatrixXcd &W0, const MatrixXcd &V0)
{
	int r = ERC.rows();
	MatrixXcd temp(r,r), P(2*r,2*r), Q(2*r,2*r), Omega2(2*r,2*r), W(2*r,2*r), LAM(2*r,2*r);
	MatrixXcd V(2*r,2*r), aa(2*r,2*r), bb(2*r,2*r), A(2*r,2*r), B(2*r,2*r), X(2*r,2*r);
	scatMat S;

	temp = ERC.lu().solve(Ky);	
	P.topLeftCorner(r,r)     = Kx*(temp);
	P.bottomLeftCorner(r,r)  = Ky*(temp) - URC;
	temp = ERC.lu().solve(Kx);	
	P.topRightCorner(r,r)    = URC - Kx*(temp);
	P.bottomRightCorner(r,r) = -Ky*(temp);

	temp = URC.lu().solve(Ky);	
	Q.topLeftCorner(r,r)     = Kx*(temp);
	Q.bottomLeftCorner(r,r)  = Ky*(temp) - ERC;
	temp = URC.lu().solve(Kx);	
	Q.topRightCorner(r,r)    = ERC - Kx*(temp);
	Q.bottomRightCorner(r,r) = -Ky*(temp);

	Omega2 = P*Q;	

	ComplexEigenSolver<MatrixXcd> ces(Omega2);

	LAM = ces.eigenvalues().asDiagonal();
	W   = ces.eigenvectors();
	
	LAM.noalias() = LAM.sqrt();

	V.noalias() = (LAM.transpose().lu().solve((Q*W).transpose())).transpose();

	aa = W.lu().solve(W0);
	bb = V.lu().solve(V0);

	A = aa + bb;
	B = aa - bb;
	X = -LAM*k0*L;
	X.noalias() = X.exp();

	S.s22 = A-X*B*(A.lu().solve(X*B)); 

	S.s11 = S.s22.lu().solve(X*B*(A.lu().solve(X*A))-B);
	S.s12 = (S.s22.lu().solve(X))*(A-B*(A.lu().solve(B)));
	S.s22 = S.s11;
	S.s21 = S.s12;

	return S;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
scatMat redheffer_star_product(const scatMat &SA, const scatMat &SB, const MatrixXd &I)
{
	scatMat S;

	S.s11 = SA.s11 + SA.s12*(((1.0+0.0i)*I-SB.s11*SA.s22).lu().solve(SB.s11*SA.s21));
	S.s12 = SA.s12*(((1.0+0.0i)*I-SB.s11*SA.s22).lu().solve(SB.s12));
	S.s21 = SB.s21*(((1.0+0.0i)*I-SA.s22*SB.s11).lu().solve(SA.s21));
	S.s22 = SB.s22 + SB.s21*(((1.0+0.0i)*I-SA.s22*SB.s11).lu().solve(SA.s22*SB.s12));
	
	return S;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
scatMat scattering_matrix_ref(complex<double> er1, complex<double> ur1, const MatrixXcd &Kx, const MatrixXcd &Ky, const MatrixXcd Kzr, const MatrixXd &Z, const MatrixXd &I, const MatrixXd &II, const MatrixXcd &W0, const MatrixXcd &V0)
{
	int r = Kx.rows();
	scatMat Sr;
	MatrixXcd Q(2*r,2*r), LAM(2*r,2*r), Vref(2*r,2*r);
	MatrixXcd aa(2*r,2*r), bb(2*r,2*r), A(2*r,2*r), B(2*r,2*r);

	Q.topLeftCorner(r,r)     = (1.0/ur1)*Kx*Ky;
	Q.topRightCorner(r,r)    = (1.0/ur1)*ur1*er1*I-Kx*Kx;
	Q.bottomLeftCorner(r,r)  = (1.0/ur1)*Ky*Ky-ur1*er1*I;
	Q.bottomRightCorner(r,r) = -(1.0/ur1)*Ky*Kx;
	
	LAM.topLeftCorner(r,r)     = (0.0-1.0i)*Kzr;
	LAM.topRightCorner(r,r)    = (1.0+0.0i)*Z; 
	LAM.bottomLeftCorner(r,r)  = (1.0+0.0i)*Z; 
	LAM.bottomRightCorner(r,r) = (0.0-1.0i)*Kzr; 

	Vref.noalias() = (LAM.transpose().lu().solve(Q.transpose())).transpose();

	aa = W0.lu().solve(W0);
	bb = V0.lu().solve(Vref);

	A = aa + bb;
	B = aa - bb;

	Sr.s11 = -A.lu().solve(B);
	Sr.s12.noalias() = (A.transpose().lu().solve(((2+0.0i)*II).transpose())).transpose();
	Sr.s21 = 0.5*(A-B*(-Sr.s11));
	Sr.s22.noalias() = (A.transpose().lu().solve(B.transpose())).transpose();

	return Sr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
scatMat scattering_matrix_trn(complex<double> er2, complex<double> ur2, const MatrixXcd &Kx, const MatrixXcd &Ky, const MatrixXcd Kzt, const MatrixXd &Z, const MatrixXd &I, const MatrixXd &II, const MatrixXcd &W0, const MatrixXcd &V0)
{
	int r = Kx.rows();
	scatMat St;
	MatrixXcd Q(2*r,2*r), LAM(2*r,2*r), Vtrn(2*r,2*r);
	MatrixXcd aa(2*r,2*r), bb(2*r,2*r), A(2*r,2*r), B(2*r,2*r);

	Q.topLeftCorner(r,r)     = (1.0/ur2)*Kx*Ky;
	Q.topRightCorner(r,r)    = (1.0/ur2)*ur2*er2*I-Kx*Kx;
	Q.bottomLeftCorner(r,r)  = (1.0/ur2)*Ky*Ky-ur2*er2*I;
	Q.bottomRightCorner(r,r) = -(1.0/ur2)*Ky*Kx;
	
	LAM.topLeftCorner(r,r)     = (0.0+1.0i)*Kzt;
	LAM.topRightCorner(r,r)    = (1.0+0.0i)*Z; 
	LAM.bottomLeftCorner(r,r)  = (1.0+0.0i)*Z; 
	LAM.bottomRightCorner(r,r) = (0.0+1.0i)*Kzt; 
	
	Vtrn.noalias() = (LAM.transpose().lu().solve(Q.transpose())).transpose();

	aa = W0.lu().solve(W0);
	bb = V0.lu().solve(Vtrn);

	A = aa + bb;
	B = aa - bb;

	St.s11.noalias() = (A.transpose().lu().solve(B.transpose())).transpose();
	St.s21.noalias() = (A.transpose().lu().solve(((2+0.0i)*II).transpose())).transpose();
	St.s22 = -A.lu().solve(B);
	St.s12 = 0.5*(A-B*(-St.s22));

	return St;
}

Vector4cd vector_cross_product(const Vector3cd &v1, const Vector3cd &v2)
{
	Vector4cd v;

	v(0) = v1(1)*v2(2)-v1(2)*v2(1);
	v(1) = v1(3)*v2(0)-v1(0)*v2(2);
	v(2) = v1(0)*v2(1)-v1(1)*v2(0);
	v(3) = sqrt(v(0)*v(0)+v(1)*v(1)+v(2)*v(2));

	return v;
}

#endif
