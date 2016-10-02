#ifndef LAYER_H
#define LAYER_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <complex>
#include <cmath>
#include "Eigen_fftw3.h"
#include <Eigen/Sparse>

using namespace std;

//! Layer class.
/*!
	Class to create layers for use in rcwa.
*/


// input data struct for grating layer
typedef struct input{
	int n;
	string l_name[20];
	double x_period[20];
	double y_period[20];
	double thickness[20];
	complex<double> epsilon[20];
	complex<double> mu[20];
} GratingParam;

// input data struct for a plane layer
typedef struct input1{
	int n;
	string l_name;
	double x_period;
	double y_period;
	double thickness;
	complex<double> epsilon;
	complex<double> mu;
} LayerParam;

// class to define geometry of device
class Layer
{
	private:
		int n;
		string l_name[20];
		double x_period[20];
		double y_period[20];
		double thickness[20];
		complex<double> epsilon[20];
		complex<double> mu[20];
		MatrixXcd layerEpsilon;
		MatrixXcd layerMu;
		MatrixXcd ep_convmat;
		MatrixXcd mu_convmat;

	public:
		Layer() {};
		// constructor for a grating layer
		Layer(const GratingParam &gratingParam);
		// constructor for a plane layer
		Layer(const LayerParam &layerParam);
		// function to display layer properties
		void showLayerProperties();
		// function to set epsilon (relative permittivity) in a grid
		void setEpsilonAndMu(int Nx, int Ny);
		// function to display the epsilon grid
		void showEpsilon();

		void showMu();

		void showEpConvmat();

		void showMuConvmat();

		string getName(int indx);

		double getThickness(int indx);

		complex<double> getEp(int indx);

		complex<double> getEp();

		complex<double> getPr();

		// function to get the epsilon matrix
		MatrixXcd getEpsilon();
		// function to get the mu matrix
		MatrixXcd getMu();

		MatrixXcd getEpConvmat();

		MatrixXcd getMuConvmat();
		// function to compute convolution matrix for the epsilon matrix
		void convMat(int P, int Q);
};

// function declarations

//////////////////////////////////////////////////////////////////////
// constructor for a grating layer
inline Layer::Layer(const GratingParam &gratingParam)
{
	if (gratingParam.n > 20)
	{
		cout << "*** ERROR :: IndexParam.n out of bound. Maximum value allowed is 20. ***" << endl;
		exit(EXIT_FAILURE);
	}

	n = gratingParam.n;
	for (int i=0; i<n; i++)
	{
		l_name[i]    = gratingParam.l_name[i];
		x_period[i]  = gratingParam.x_period[i];
		y_period[i]  = gratingParam.y_period[i];
		thickness[i] = gratingParam.thickness[i];
		epsilon[i]   = gratingParam.epsilon[i];
		mu[i]        = gratingParam.mu[i];
	}
}

//////////////////////////////////////////////////////////////////////
// constructor for a plane layer
inline Layer::Layer(const LayerParam &layerParam)
{
	n            = layerParam.n;
	l_name[0]    = layerParam.l_name;
	x_period[0]  = layerParam.x_period;
	y_period[0]  = layerParam.y_period;
	thickness[0] = layerParam.thickness;
	epsilon[0]   = layerParam.epsilon;
	mu[0]        = layerParam.mu;
}

//////////////////////////////////////////////////////////////////////
// function to display layer properties
inline void Layer::showLayerProperties()
{
	if (n==1) // for a plane layer
	{
		cout << "n         = " << n << endl;
		cout << "l_name    = " << l_name[0] << endl;
		cout << "x_period  = " << x_period[0] << endl;
		cout << "y_period  = " << y_period[0] << endl;
		cout << "thickness = " << thickness[0] << endl;
		cout << "epsilon   = " << epsilon[0] << endl;
		cout << "mu        = " << mu[0] << endl;
	}
	else      // for a grating layer
	{
		cout << "n         = " << n << endl;
		cout << "l_name    = { ";
		for (int i=0; i<n; i++)
			cout << l_name[i] << " ";
		cout << "}" << endl;
		cout << "x_period  = { ";
		for (int i=0; i<n; i++)
			cout << x_period[i] << " ";
		cout << "}" << endl;
		cout << "y_period  = { ";
		for (int i=0; i<n; i++)
			cout << y_period[i] << " ";
		cout << "}" << endl;
		cout << "thickness = { ";
		for (int i=0; i<n; i++)
			cout << thickness[i] << " ";
		cout << "}" << endl;
		cout << "epsilon   = { ";
		for (int i=0; i<n; i++)
			cout << epsilon[i] << " ";
		cout << "}" << endl;
		cout << "mu        = { ";
		for (int i=0; i<n; i++)
			cout << mu[i] << " ";
		cout << "}" << endl;
	}
}

inline string Layer::getName(int indx)
{
	return l_name[indx];
}


inline double Layer::getThickness(int indx)
{
	return thickness[indx];
}

inline complex<double> Layer::getEp(int indx)
{
	return epsilon[indx];
}


inline complex<double> Layer::getEp()
{
	return epsilon[0];
}

inline complex<double> Layer::getPr()
{
	return mu[0];
}

//////////////////////////////////////////////////////////////////////
// function to set epsilon (relative permittivity) in a grid
inline void Layer::setEpsilonAndMu(int Nx, int Ny)
{
	layerEpsilon.resize(Nx, Ny);
	layerMu.resize(Nx, Ny);
	// set epsilon for a plane layer
	for (int i=0; i<Nx; i++)
		for (int j=0; j<Ny; j++)
		{
			layerEpsilon(i,j) = epsilon[0];
			layerMu(i,j)      = mu[0];
		}
	// set epsilon for a grating layer
	if (n==3)
	{
		double fillout_factor = x_period[1]/(x_period[0]+x_period[1]+x_period[2]);
		int nx  = (int) round(fillout_factor*Nx);
		int nx1 = (int) floor((Nx-nx)/2);
		int nx2 = nx1 + nx;
		// overlayer grating epsilon on a plane layer
		for (int i=nx1; i<nx2; i++)
			for (int j=0; j<Ny; j++)
			{
				layerEpsilon(i,j) = epsilon[1];
				layerMu(i,j)      = mu[1];
			}
	}
}	

//////////////////////////////////////////////////////////////////////
// function to display the epsilon grid
inline void Layer::showEpsilon()
{
	cout << layerEpsilon << endl;
}

//////////////////////////////////////////////////////////////////////
// function to display the epsilon grid
inline void Layer::showMu()
{
	cout << layerMu << endl;
}

//////////////////////////////////////////////////////////////////////
// function to display the epsilon grid
inline void Layer::showEpConvmat()
{
	cout << ep_convmat << endl;
}

//////////////////////////////////////////////////////////////////////
// function to display the epsilon grid
inline void Layer::showMuConvmat()
{
	cout << mu_convmat << endl;
}

//////////////////////////////////////////////////////////////////////
// function to get the epsilon matrix
inline MatrixXcd Layer::getEpsilon()
{
	return layerEpsilon;
}

//////////////////////////////////////////////////////////////////////
// function to get the epsilon matrix
inline MatrixXcd Layer::getEpConvmat()
{
	return ep_convmat;
}

//////////////////////////////////////////////////////////////////////
// function to get the epsilon matrix
inline MatrixXcd Layer::getMuConvmat()
{
	return mu_convmat;
}

//////////////////////////////////////////////////////////////////////
// function to get the epsilon matrix
inline MatrixXcd Layer::getMu()
{
	return layerMu;
}

//////////////////////////////////////////////////////////////////////
// function to compute convolution matrix for the epsilon matrix
inline void Layer::convMat(int P, int Q)
{
	int Nx, Ny, NH, p0, q0, row, col, pfft, qfft;
	Eigen_fftw3 fftw3_obj;
	MatrixXcd A, B;
	RowVectorXi p, q;

	
	Nx = layerEpsilon.rows(); // # of rows in the epsilon matrix
	Ny = layerEpsilon.cols(); // # of cols in the epsilon matrix

	if (fmod(P,2)==0 && fmod(Q,2)==0)
	{
		NH = (P+1) * (Q+1);
		p.resize(P+1);
		q.resize(Q+1);
	}
	else if (fmod(P,2)!=0 && fmod(Q,2)!=0)
	{
		NH = P * Q;
		p.resize(P);
		q.resize(Q);
	}
	//else if (fmod(P,2)==0 && fmod(Q,2)!=0)
	//{
	//	NH = (P+1) * Q;
	//	p.resize(P+1);
	//	q.resize(Q);
	//}
	//else if (fmod(P,2)!=0 && fmod(Q,2)==0)
	//{
	//	NH = P * (Q+1);
	//	p.resize(P);
	//	q.resize(Q+1);
	//}

	ep_convmat.resize(NH,NH); // Matrix for convolution
	mu_convmat.resize(NH,NH); // Matrix for convolution

	// compute indices of spatial harmonics
	p(0) = -(int)floor(P/2.0);
	for (int i=1; i<p.size(); i++)
		p(i) = p(i-1)+1;

	q(0) = -(int)floor(Q/2.0);
	for (int i=1; i<q.size(); i++)
		q(i) = q(i-1)+1;


	// compute Fourier coefficients of the epsilon matrix
	A = fftw3_obj.fftShift2d(fftw3_obj.fft2d(layerEpsilon))/(Nx*Ny);
	B = fftw3_obj.fftShift2d(fftw3_obj.fft2d(layerMu))/(Nx*Ny);

	//cout << A.sparseView() << endl;
	
	// compute array indices of center harmonics
	p0 = 1 + (int)floor(Nx/2);
	q0 = 1 + (int)floor(Ny/2);

	//cout << p0 << endl;
	//cout << q0 << endl;

	// compute convolution of the epsilon matrix
	for (int qrow=0; qrow<Q; qrow++)
		for (int prow=0; prow<P; prow++)
		{
			row = qrow*P + prow;
			for (int qcol=0; qcol<Q; qcol++)
				for (int pcol=0; pcol<P; pcol++)
				{
					col = qcol*P + pcol;
					pfft = p(prow) - p(pcol);
					qfft = q(qrow) - q(qcol);
					ep_convmat(row,col) = A(p0+pfft-1,q0+qfft-1); 
					mu_convmat(row,col) = B(p0+pfft-1,q0+qfft-1); 
				}
		}
}

#endif
