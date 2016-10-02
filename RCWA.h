#ifndef RCWA_H
#define RCWA_H

#include "Source.h"
#include "Layer.h"
#include "RCWA_utils.h"

//! RCWA Class
/*!
	Create a RCWA object that intializes the variables, vectors, and matrices used
	during the computation.
	It defines methods to set up device layers, set material properties, set source,
	and run rigorous coupled wave analysis.
*/

class RCWA
{
	private:
		string name;
		double x_period;
		double y_period;
		Source lightSource;
		int    numLayers;
		Layer  *layers;
		PQ     pq;
		Grid   grid;
		double reflection;
		double transmission;

		MatrixXd I, II;
		MatrixXd Z, ZZ;
		complex<double> n1;
		double k0;

		Vector3cd kinc;
		RowVectorXi p, q;
		RowVectorXcd kx, ky;

		MatrixXcd Kx, Ky, Kzr, Kzt, Q, W0, V0;
		SparseMatrix<complex<double> > KKx, KKy;

		int r1, r2, r3, r4, c1, c2, c3, c4;
		scatMat S, Sg, Sr, St, Sgg;

		VectorXcd delta, esrc, eref, etrn, csrc, cref, ctrn;
		VectorXcd rx, ry, rz, tx, ty, tz;
		VectorXd ref, trn, RR, TT;
		Vector3cd n_hat, a_hat_te, a_hat_tm, EP;
		Vector4cd vec_cross;

		int pte, ptm;


	public:
		//RCWA(string nm = "rcwa", int nl = 2, PQ p = {11,11}, Grid g = {512, 512}, double Lx = 500e-9, double Ly = 500e-9);
		RCWA(string nm, int nl, const PQ &p, const Grid &g, double Lx, double Ly);
		void setLayer(const LayerParam &lparam);
		void addGratingLayer(int n, const GratingParam &gparam);
		void addPlaneLayer(int n, const LayerParam &lparam);
		void displayErConvmat(int n);
		void displayUrConvmat(int n);
		void displayLayerProperties(int n);
		void displayEr(int n);
		void displayUr(int n);
		double getXPeriod();
		double getYPeriod();
		double getReflection();
		double getTransmission();
		void applyLightSource(string type, double wavelen, double theta, double phi);
		MatrixXcd getEr(int n);
		MatrixXcd getUr(int n);
		MatrixXcd getErConvmat(int n);
		MatrixXcd getUrConvmat(int n);
		void allocate();
		void run();

};

// function declarations

///////////////////////////////////////////////////////////////////////////////////////////////
inline RCWA::RCWA(string nm, int nl, const PQ &p, const Grid &g, double Lx, double Ly)
{
	name      = nm;
	numLayers = nl;
	pq        = p;
	grid      = g;
	x_period  = Lx;
  y_period  = Ly;
	r1        = pq.P;
	r2        = pq.P+1;
	r3        = pq.P*pq.P;
	r4        = (pq.P+1)*(pq.P+1);
	c1        = pq.Q;
	c2        = pq.Q+1;
	c3        = pq.Q*pq.Q;
	c4        = (pq.Q+1)*(pq.Q+1);
	layers    = new Layer [numLayers+2];
}


///////////////////////////////////////////////////////////////////////////////////////////////
inline void RCWA::setLayer(const LayerParam &lparam)
{
	if (lparam.l_name.compare("in")==0)
	{
		layers[0]	= lparam;
		layers[0].setEpsilonAndMu(grid.Nx, grid.Ny);
		layers[0].convMat(r1, c1);
	}
	else if (lparam.l_name.compare("out")==0)
	{
		layers[numLayers+1]	= lparam;
		layers[numLayers+1].setEpsilonAndMu(grid.Nx, grid.Ny);
		layers[numLayers+1].convMat(r1, c1);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline void RCWA::addGratingLayer(int n, const GratingParam &gparam)
{
	layers[n] = gparam;
	layers[n].setEpsilonAndMu(grid.Nx, grid.Ny);
	layers[n].convMat(r1, c1);
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline void RCWA::addPlaneLayer(int n, const LayerParam &lparam)
{
	layers[n] = lparam;
	layers[n].setEpsilonAndMu(grid.Nx, grid.Ny);
	layers[n].convMat(r1, c1);	
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline void RCWA::displayErConvmat(int n)
{
	layers[n].showEpConvmat();
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline void RCWA::displayUrConvmat(int n)
{
	layers[n].showMuConvmat();
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline void RCWA::displayLayerProperties(int n)
{
	layers[n].showLayerProperties();
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline void RCWA::displayEr(int n)
{
	layers[n].showEpsilon();
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline void RCWA::displayUr(int n)
{
	layers[n].showMu();
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline double RCWA::getXPeriod()
{
	return x_period;
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline double RCWA::getYPeriod()
{
	return y_period;
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline double RCWA::getReflection()
{
	return reflection;
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline double RCWA::getTransmission()
{
	return transmission;
}
///////////////////////////////////////////////////////////////////////////////////////////////
inline void RCWA::applyLightSource(string type, double wavelen, double theta, double phi)
{
	lightSource.setSourcePolarization(type);
	lightSource.setWavelength(wavelen);
	lightSource.setTheta(theta);
	lightSource.setPhi(phi);
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline MatrixXcd RCWA::getEr(int n)
{
	return layers[n].getEpsilon();
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline MatrixXcd RCWA::getUr(int n)
{
	return layers[n].getMu();
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline MatrixXcd RCWA::getErConvmat(int n)
{
	return layers[n].getEpConvmat();
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline MatrixXcd RCWA::getUrConvmat(int n)
{
	return layers[n].getMuConvmat();
}

///////////////////////////////////////////////////////////////////////////////////////////////
inline void RCWA::allocate()
{
	if (fmod(r1,2)==0 && fmod(c1,2)==0)
	{
		p.resize(r2);
		q.resize(c2);
		kx.resize(r2);
		ky.resize(c2);
		I.resize(r4, c4);
		Z.resize(r4, c4);
		Kx.resize(r4, c4);
		Ky.resize(r4, c4);
		Kzr.resize(r4, c4);
		Kzt.resize(r4, c4);
		Q.resize(2*r4, 2*c4);
		W0.resize(2*r4, 2*c4);
		V0.resize(2*r4, 2*c4);
		II.resize(2*r4, 2*c4);
		ZZ.resize(2*r4, 2*c4);
		S.s11.resize(2*r4, 2*c4);
		S.s12.resize(2*r4, 2*c4);
		S.s21.resize(2*r4, 2*c4);
		S.s22.resize(2*r4, 2*c4);
		Sg.s11.resize(2*r4, 2*c4);
		Sg.s12.resize(2*r4, 2*c4);
		Sg.s21.resize(2*r4, 2*c4);
		Sg.s22.resize(2*r4, 2*c4);
		Sr.s11.resize(2*r4, 2*c4);
		Sr.s12.resize(2*r4, 2*c4);
		Sr.s21.resize(2*r4, 2*c4);
		Sr.s22.resize(2*r4, 2*c4);
		St.s11.resize(2*r4, 2*c4);
		St.s12.resize(2*r4, 2*c4);
		St.s21.resize(2*r4, 2*c4);
		St.s22.resize(2*r4, 2*c4);
		Sgg.s11.resize(2*r4, 2*c4);
		Sgg.s12.resize(2*r4, 2*c4);
		Sgg.s21.resize(2*r4, 2*c4);
		Sgg.s22.resize(2*r4, 2*c4);
		delta.resize(r4);
		esrc.resize(2*r4);
		eref.resize(2*r4);
		etrn.resize(2*r4);
		csrc.resize(2*r4);
		cref.resize(2*r4);
		ctrn.resize(2*r4);
		rx.resize(r4);
		ry.resize(r4);
		rz.resize(r4);
		tx.resize(r4);
		ty.resize(r4);
		tz.resize(r4);
		ref.resize(r4);
		RR.resize(r4);
		trn.resize(r4);
		TT.resize(r4);
	}
	else if (fmod(r1,2)!=0 && fmod(c1,2)!=0)
	{
		p.resize(r1);
		q.resize(c1);
		kx.resize(r1);
		ky.resize(c1);
		I.resize(r3, c3);
		Z.resize(r3, c3);
		Kx.resize(r3, c3);
		Ky.resize(r3, c3);
		Kzr.resize(r3, c3);
		Kzt.resize(r3, c3);
		Q.resize(2*r3, 2*c3);
		W0.resize(2*r3, 2*c3);
		V0.resize(2*r3, 2*c3);
		II.resize(2*r3, 2*c3);
		ZZ.resize(2*r3, 2*c3);
		S.s11.resize(2*r3, 2*c3);
		S.s12.resize(2*r3, 2*c3);
		S.s21.resize(2*r3, 2*c3);
		S.s22.resize(2*r3, 2*c3);
		Sg.s11.resize(2*r3, 2*c3);
		Sg.s12.resize(2*r3, 2*c3);
		Sg.s21.resize(2*r3, 2*c3);
		Sg.s22.resize(2*r3, 2*c3);
		Sr.s11.resize(2*r3, 2*c3);
		Sr.s12.resize(2*r3, 2*c3);
		Sr.s21.resize(2*r3, 2*c3);
		Sr.s22.resize(2*r3, 2*c3);
		St.s11.resize(2*r3, 2*c3);
		St.s12.resize(2*r3, 2*c3);
		St.s21.resize(2*r3, 2*c3);
		St.s22.resize(2*r3, 2*c3);
		Sgg.s11.resize(2*r3, 2*c3);
		Sgg.s12.resize(2*r3, 2*c3);
		Sgg.s21.resize(2*r3, 2*c3);
		Sgg.s22.resize(2*r3, 2*c3);
		delta.resize(r3);
		esrc.resize(2*r3);
		eref.resize(2*r3);
		etrn.resize(2*r3);
		csrc.resize(2*r3);
		cref.resize(2*r3);
		ctrn.resize(2*r3);
		rx.resize(r3);
		ry.resize(r3);
		rz.resize(r3);
		tx.resize(r3);
		ty.resize(r3);
		tz.resize(r3);
		ref.resize(r3);
		RR.resize(r3);
		trn.resize(r3);
		TT.resize(r3);
	}

	I.setIdentity();
	Z.setZero();
	II.setIdentity();
	ZZ.setZero();

	p(0) = -(int)floor(r1/2.0);
	for (int i=1; i<p.size(); i++)
		p(i) = p(i-1)+1;

	q(0) = -(int)floor(c1/2.0);
	for (int i=1; i<q.size(); i++)
		q(i) = q(i-1)+1;


	delta.setZero();
	delta(ceil(delta.size()/2)) = 1.0;

	n_hat << 0, 0, -1;

	esrc.setZero();
	//cout << delta << endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////
inline void RCWA::run()
{
	reflection   = 0;
	transmission = 0; 	
	n1 = sqrt(layers[0].getEp());
	k0 = 2*pi_const/lightSource.getWavelength();

	kinc << n1*sin(lightSource.getTheta())*cos(lightSource.getPhi()),
					n1*sin(lightSource.getTheta())*sin(lightSource.getPhi()),
					n1*cos(lightSource.getTheta());

	for (int i=0; i<p.size(); i++)	
		kx(i) = kinc(0) - 2*pi_const*p(i)/(k0*x_period);
	for (int i=0; i<q.size(); i++)	
		ky(i) = kinc(1) - 2*pi_const*q(i)/(k0*y_period);

	for (int i=0; i<p.size(); i++)
		for (int j=0; j<q.size(); j++)
			Kx(q.size()*i+j,q.size()*i+j) = kx(j);
	for (int i=0; i<p.size(); i++)
		for (int j=0; j<q.size(); j++)
			Ky(q.size()*i+j,q.size()*i+j) = ky(i);
 
	Kzr = -(layers[0].getPr()*layers[0].getEp()*I-Kx*Kx-Ky*Ky).array().sqrt().conjugate();
	Kzt = (layers[numLayers+1].getPr()*layers[numLayers+1].getEp()*I-Kx*Kx-Ky*Ky).array().sqrt().conjugate();

	Q.topLeftCorner(p.size()*p.size(),q.size()*q.size())      = Kx*Ky;	
	Q.topRightCorner(p.size()*p.size(),q.size()*q.size())     = (1.0+0.0i)*I+Ky*Ky;	
	Q.bottomLeftCorner(p.size()*p.size(),q.size()*q.size())   = -((1.0+0.0i)*I+Kx*Kx);	
	Q.bottomRightCorner(p.size()*p.size(),q.size()*q.size())  = -Kx*Ky;	

	W0.topLeftCorner(p.size()*p.size(),q.size()*q.size())     = (1.0+0.0i)*I;	
	W0.topRightCorner(p.size()*p.size(),q.size()*q.size())    = (1.0+0.0i)*Z;	
	W0.bottomLeftCorner(p.size()*p.size(),q.size()*q.size())  = (1.0+0.0i)*Z;	
	W0.bottomRightCorner(p.size()*p.size(),q.size()*q.size()) = (1.0+0.0i)*I;	

	V0 = -1.0i*Q;

	Sg.s11 = (1.0+0.0i)*ZZ;
	Sg.s12 = (1.0+0.0i)*II;
	Sg.s21 = (1.0+0.0i)*II;
	Sg.s22 = (1.0+0.0i)*ZZ;

	for (int i=1; i<=numLayers; i++)
	{
		S = scattering_matrix_ith(layers[i].getEpConvmat(), layers[i].getMuConvmat(), Kx, Ky, k0, layers[i].getThickness(0), W0, V0);
		Sg = redheffer_star_product(Sg, S, II);
	}

	Sr = scattering_matrix_ref(layers[0].getEp(), layers[0].getPr(), Kx, Ky, Kzr, Z, I, II, W0, V0);
	St = scattering_matrix_trn(layers[numLayers+1].getEp(), layers[numLayers+1].getPr(), Kx, Ky, Kzt, Z, I, II, W0, V0);

	Sgg  = redheffer_star_product(Sr, Sg, II);
	Sgg  = redheffer_star_product(Sgg, St, II);

	vec_cross = vector_cross_product(kinc/k0, n_hat);

	if (lightSource.getTheta() == 0)
		a_hat_te << 0, 1, 0;
	else
		a_hat_te << vec_cross(0)/vec_cross(3), vec_cross(1)/vec_cross(3), vec_cross(2)/vec_cross(3);

	vec_cross = vector_cross_product(a_hat_te, kinc/k0);

	a_hat_tm << vec_cross(0)/vec_cross(3), vec_cross(1)/vec_cross(3), vec_cross(2)/vec_cross(3);

	if (lightSource.getSourcePolarization().compare("TE")==0)
	{
		pte = 1;
		ptm = 0;
	}
	else if (lightSource.getSourcePolarization().compare("TM")==0)
	{
		pte = 0;
		ptm = 1;
	}

	EP = pte*a_hat_te + ptm*a_hat_tm;
	
	esrc.segment(0,esrc.size()/2) = EP(0)*delta;	
	esrc.segment(esrc.size()/2, esrc.size()/2) = EP(1)*delta;	

	csrc = W0.lu().solve(esrc);

	cref = Sgg.s11*csrc;
	ctrn = Sgg.s21*csrc;

	eref = W0*cref;
	etrn = W0*ctrn;

	rx = eref.segment(0,eref.size()/2);
	ry = eref.segment(eref.size()/2 ,eref.size()/2);
	rz = -Kzr.lu().solve(Kx*rx + Ky*ry);

	tx = etrn.segment(0,etrn.size()/2);
	ty = etrn.segment(etrn.size()/2 ,etrn.size()/2);
	tz = -Kzt.lu().solve(Kx*tx + Ky*ty);

	ref = rx.array().abs()*rx.array().abs() + ry.array().abs()*ry.array().abs() + rz.array().abs()*rz.array().abs();
	RR = (-Kzr/kinc(2)).real()*ref;
	reflection = RR.array().sum();

	trn = tx.array().abs()*tx.array().abs() + ty.array().abs()*ty.array().abs() + tz.array().abs()*tz.array().abs();
	TT = ((layers[0].getPr()/layers[numLayers+1].getPr())*Kzt/kinc(2)).real()*trn;
	transmission = TT.array().sum();
}

#endif
