#ifndef SOURCE_H
#define SOURCE_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <Eigen/Dense>

using namespace Eigen;

using namespace std;

// class to define Input source
class Source
{
	private:
		string s_type;
		string s_pol;
		int    n_wavelength;
		RowVectorXd wavelength;
		int    n_theta;
		RowVectorXd theta;
		int    n_phi;
		RowVectorXd phi;

	public:
		Source(string type = "PlaneWave", string pol = "TE");
		void setSourceType(string type);
		void setSourcePolarization(string pol);
		void setWavelength(double wavelen);
		void setWavelengthRange(double min_wavelength, double max_wavelength, int num_wavelength);
		void setTheta(double th);
		void setThetaRange(double min_theta, double max_theta, int num_theta);
		void setPhi(double ph);
		void setPhiRange(double min_phi, double max_phi, int num_phi);
		string getSourceType();
		string getSourcePolarization();
		double getWavelength();
		double getWavelength(int index);
		double getTheta();
		double getTheta(int index);
		double getPhi();
		double getPhi(int index);
		void   showWavelength();
		void   showTheta();
		void   showPhi();
};

// function declarations

//////////////////////////////////////////////////////////////////////
// constructor
inline Source::Source(string type, string pol) 
{
	s_type = type;
	s_pol = pol;
};

//////////////////////////////////////////////////////////////////////
// set source type
inline void Source::setSourceType(string type) { s_type = type; }

//////////////////////////////////////////////////////////////////////
// set source polarization
inline void Source::setSourcePolarization(string pol) { s_pol = pol; }

//////////////////////////////////////////////////////////////////////
// set a single wavlength
inline void Source::setWavelength(double wavelen) 
{
	wavelength.resize(1);
	wavelength(0) = wavelen;
}

//////////////////////////////////////////////////////////////////////
// set source wavelength
inline void Source::setWavelengthRange(double min_wavelength, double max_wavelength, int num_wavelength)
{
	n_wavelength = num_wavelength;
	wavelength.resize(n_wavelength);
	for (int i=0; i<n_wavelength; i++)
		wavelength(i) = min_wavelength + i*(max_wavelength-min_wavelength)/(n_wavelength-1);
}

//////////////////////////////////////////////////////////////////////
// set a single wavlength
inline void Source::setTheta(double th) 
{
	theta.resize(1);
	theta(0) = th;
}

//////////////////////////////////////////////////////////////////////
// set source theta
inline void Source::setThetaRange(double min_theta, double max_theta, int num_theta)
{
	n_theta = num_theta;
	theta.resize(n_theta);
	for (int i=0; i<n_theta; i++)
		theta(i) = min_theta + i*(max_theta-min_theta)/(n_theta-1);
}

//////////////////////////////////////////////////////////////////////
// set a single wavlength
inline void Source::setPhi(double ph) 
{
	phi.resize(1);
	phi(0) = ph;
}

//////////////////////////////////////////////////////////////////////
// set source phi
inline void Source::setPhiRange(double min_phi, double max_phi, int num_phi)
{
	n_phi = num_phi;
	phi.resize(n_phi);
	for (int i=0; i<n_phi; i++)
		phi(i) = min_phi + i*(max_phi-min_phi)/(n_phi-1);
}

//////////////////////////////////////////////////////////////////////
// get source type
inline string Source::getSourceType() { return s_type; }

//////////////////////////////////////////////////////////////////////
// get source polarization
inline string Source::getSourcePolarization() { return s_pol; }

//////////////////////////////////////////////////////////////////////
// get source wavelength
inline double Source::getWavelength()
{
	return wavelength(0);
}

//////////////////////////////////////////////////////////////////////
// get source wavelength range
inline double Source::getWavelength(int index)
{
	if (index > n_wavelength-1)
	{
		cout << "*** ERROR :: Wavelength range :: Out of bound. ***" << endl;
		exit(EXIT_FAILURE);
	}

	return wavelength(index);
}

//////////////////////////////////////////////////////////////////////
// get source theta
inline double Source::getTheta()
{
	return theta(0);
}

//////////////////////////////////////////////////////////////////////
// get source theta range
inline double Source::getTheta(int index)
{
	if (index > n_theta-1)
	{
		cout << "*** ERROR :: Theta range :: Out of bound. ***" << endl;
		exit(EXIT_FAILURE);
	}

	return theta(index);
}

//////////////////////////////////////////////////////////////////////
// get source phi
inline double Source::getPhi()
{
	return phi(0);
}

//////////////////////////////////////////////////////////////////////
// get source phi range
inline double Source::getPhi(int index)
{
	if (index > n_phi-1)
	{
		cout << "*** ERROR :: Phi range :: Out of bound. ***" << endl;
		exit(EXIT_FAILURE);
	}

	return phi(index);
}

//////////////////////////////////////////////////////////////////////
// print wavelenth values
inline void Source::showWavelength()
{
	if (wavelength.size() == 1)
		cout << "Wavelength : " << endl;
	else
		cout << "Wavelength range : " << endl;
	cout << wavelength << endl;
}

//////////////////////////////////////////////////////////////////////
// print theta values
inline void Source::showTheta()
{
	if (theta.size() == 1)
		cout << "Angle : " << endl;
	else
		cout << "Angle range : " << endl;
	cout << theta << endl;
}

//////////////////////////////////////////////////////////////////////
// print phi values
inline void Source::showPhi()
{
	if (phi.size() == 1)
		cout << "Phi : " << endl;
	else
		cout << "Phi range : " << endl;
	cout << phi << endl;
}

#endif
