#ifndef FFTW3_UTILS_H
#define FFTW3_UTILS_H

#include "../include/fftw3.h"
#include <Eigen/Dense>

using namespace Eigen;

// class to compute 2d fft using fftw3 and Eigen.
// The output matrices are in the form of Eigen matrices.
class Eigen_fftw3
{
	public:
		//Eigen_fftw3(MatrixXcd mat);
		// function to compute 2d-fft of a matrix
		MatrixXcd fft2d(const MatrixXcd &A);
		// function to compute 2d-ifft of a matrix
		MatrixXcd ifft2d(const MatrixXcd &A);
		// function to compute 2d-fftshift (like MATLAB) of a matrix
		MatrixXcd fftShift2d(const MatrixXcd &mat);
		// function to compute 2d-ifftshift (like MATLAB) of a matrix
		MatrixXcd ifftShift2d(const MatrixXcd &mat);

};

// function declarations

//inline Eigen_fftw3::Eigen_fftw3(MatrixXcd mat){ A = mat; }

//////////////////////////////////////////////////////////////////////
// function to compute 2d-fft of a matrix
inline MatrixXcd Eigen_fftw3::fft2d(const MatrixXcd &A)
{
	int nx, ny;
	fftw_complex *in;
	fftw_complex *out;
	fftw_plan plan_forward;

	nx = A.rows();
	ny = A.cols();

	// matrix to store 2d fft data
	MatrixXcd fft2d_mat(nx,ny);

	// allocate 1d array to store matrix data
	in = (fftw_complex*)fftw_malloc (sizeof (fftw_complex) * nx * ny );

	// convert a 2d-matrix into a 1d array
	for (int i=0; i<nx; i++)
		for (int j=0; j<ny; j++)
		{
			in[i*ny+j][0] = A(i,j).real();
			in[i*ny+j][1] = A(i,j).imag();
			//cout << i << " " << j << " " << in[i*ny+j][0] << " " << in[i*ny+j][1] << endl; 
		}

	// allocate 1d array to store fft data
	out = (fftw_complex*)fftw_malloc (sizeof (fftw_complex) * nx * ny );
	
	// prepare to do forward fft
	plan_forward = fftw_plan_dft_2d (nx, ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	
	// compute forward fft
	fftw_execute(plan_forward); 

	// covert 1d fft data to a 2d-matrix
	for (int i=0; i<nx; i++)
		for (int j=0; j<ny; j++)
		{
			//cout << i << " " << j << " " << out[i*ny+j][0] << " " << out[i*ny+j][1] << endl; 
			fft2d_mat(i,j).real() = out[i*ny+j][0];	
			fft2d_mat(i,j).imag() = out[i*ny+j][1];	
		}

	// destroy fft objects
	fftw_destroy_plan(plan_forward);
	fftw_free(in);
	fftw_free(out);

	// return the 2d fft matrix
	return fft2d_mat;
}

//////////////////////////////////////////////////////////////////////
// function to compute 2d-ifft of a matrix
inline MatrixXcd Eigen_fftw3::ifft2d(const MatrixXcd &A)
{
	int nx, ny;
	fftw_complex *in;
	fftw_complex *out;
	fftw_plan plan_backward;

	nx = A.rows();
	ny = A.cols();

	// matrix to store 2d fft data
	MatrixXcd ifft2d_mat(nx,ny);

	// allocate 1d array to store matrix data
	in = (fftw_complex*)fftw_malloc (sizeof (fftw_complex) * nx * ny );

	// convert a 2d-matrix into a 1d array
	for (int i=0; i<nx; i++)
		for (int j=0; j<ny; j++)
		{
			in[i*ny+j][0] = A(i,j).real();
			in[i*ny+j][1] = A(i,j).imag();
			//cout << i << " " << j << " " << in[i*ny+j][0] << " " << in[i*ny+j][1] << endl; 
		}

	// allocate 1d array to store ifft data
	out = (fftw_complex*)fftw_malloc (sizeof (fftw_complex) * nx * ny );
	
	// prepare to do backward fft (ifft)
	plan_backward = fftw_plan_dft_2d (nx, ny, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

	// compute backward fft (ifft)	
	fftw_execute(plan_backward); 

	// covert 1d ifft data to a 2d-matrix
	for (int i=0; i<nx; i++)
		for (int j=0; j<ny; j++)
		{
			//cout << i << " " << j << " " << out[i*ny+j][0] << " " << out[i*ny+j][1] << endl; 
			ifft2d_mat(i,j).real() = out[i*ny+j][0]/(double)(nx*ny);	
			ifft2d_mat(i,j).imag() = out[i*ny+j][1]/(double)(nx*ny);	
		}

	// destroy ifft objects
	fftw_destroy_plan(plan_backward);
	fftw_free(in);
	fftw_free(out);

	// return the 2d ifft matrix
	return ifft2d_mat;
}


//////////////////////////////////////////////////////////////////////
// function to compute 2d-fftshift (like MATLAB) of a matrix
inline MatrixXcd Eigen_fftw3::fftShift2d(const MatrixXcd &mat)
{
	int m, n, p, q;
	m = mat.rows();
	n = mat.cols();
	
	// matrix to store fftshift data
	MatrixXcd mat_fftshift(m,n);

	// odd # of rows and cols
	if ((int)fmod(m,2)==1)
	{
		p = (int) floor(m/2.0);
		q = (int) floor(n/2.0);
	}
	else // even # of rows and cols
	{
		p = (int) ceil(m/2.0);
		q = (int) ceil(n/2.0);
	}

	// vectors to store swap indices
	RowVectorXi indx(m), indy(n);

	// compute swap indices
	if ((int)fmod(m,2)==1) // # of rows odd
	{
		for (int i=0; i<m-p-1; i++)
			indx(i) = (m-p)+i;
		for (int i=m-p-1; i<m; i++)
			indx(i) = i-(m-p-1);
	}
	else // # of rows even
	{
		for (int i=0; i<m-p; i++)
			indx(i) = p+i;
		for (int i=m-p; i<m; i++)
			indx(i) = i-(m-p);
	}

	if ((int)fmod(n,2)==1) // # of cols odd
	{
		for (int i=0; i<n-q-1; i++)
			indy(i) = (n-q)+i;
		for (int i=n-q-1; i<n; i++)
			indy(i) = i-(n-q-1);
	}
	else // # of cols even
	{
		for (int i=0; i<n-q; i++)
			indy(i) = q+i;
		for (int i=n-q; i<n; i++)
			indy(i) = i-(n-q);
	}

	// rearrange the matrix elements by swapping the elements
	// according to the indices computed above.
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			mat_fftshift(i,j) = mat(indx(i),indy(j));
	
	// return fftshift matrix
	return mat_fftshift;
}

//////////////////////////////////////////////////////////////////////
// function to compute 2d-ifftshift (like MATLAB) of a matrix
inline MatrixXcd Eigen_fftw3::ifftShift2d(const MatrixXcd &mat)
{
	int m, n, p, q;
	m = mat.rows();
	n = mat.cols();
	
	// matrix to store ifftshift data
	MatrixXcd mat_ifftshift(m,n);

	p = (int) ceil(m/2.0);
	q = (int) ceil(n/2.0);

	// vectors to store swap indices
	RowVectorXi indx(m), indy(n);
	
	// compute swap indices
	if ((int)fmod(m,2)==1) // # of rows odd
	{
		for (int i=0; i<=m-p; i++)
			indx(i) = (m-p)+i;
		for (int i=m-p+1; i<m; i++)
			indx(i) = i-(m-p+1);
	}
	else // # of rows even
	{
		for (int i=0; i<m-p; i++)
			indx(i) = p+i;
		for (int i=m-p; i<m; i++)
			indx(i) = i-(m-p);
	}

	if ((int)fmod(n,2)==1) // # of cols odd
	{
		for (int i=0; i<=n-q; i++)
			indy(i) = (n-q)+i;
		for (int i=n-q+1; i<n; i++)
			indy(i) = i-(n-q+1);
	}
	else // # of cols odd
	{
		for (int i=0; i<n-q; i++)
			indy(i) = q+i;
		for (int i=n-q; i<n; i++)
			indy(i) = i-(n-q);
	}

	// rearrange the matrix elements by swapping the elements
	// according to the indices computed above.
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			mat_ifftshift(i,j) = mat(indx(i),indy(j));
	
	// return ifftshift matrix	
	return mat_ifftshift;
}

#endif
