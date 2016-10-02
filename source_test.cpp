#include "Source.h"
#include "Grid.h"
#include "Layer.h"

int main()
{
	// UNITS
	double pi = 3.14159;
	double meters = 1.0;
	double nanometers = 1e-9 * meters;
	double degrees = pi/180; 

	// SOURCE PARAMETERS
	double lam0 = 632.8 * nanometers;
	double th_left = 0 * degrees;
	double th_right = 90 * degrees;
	int    th_num = 1000;
	double phi0 = 0 * degrees;
	Source s1;
	s1.setWavelength(lam0);
	s1.showWavelength();
	s1.setThetaRange(th_left,th_right,th_num);
	s1.setPhi(phi0);

	// Resolution in x and y direction
	Grid grid = {512,512};
	// # of harmonics in x and y direction
	PQ pq = {5,5};

	// DEVICE PARAMETERS
	double ur1 = 1.0;         // permeability in reflection region
	double er1 = pow(1.33,2); // permittivity in reflection region
	double ur2 = 1.0;         // permeability in transmission region
	double er2 = pow(1.5,2);  // permittivity in transmission region

	int n_layers_unitcell = 2;
	int b_bilayers = 8;

	RowVectorXd urd(n_layers_unitcell);
	RowVectorXcd erd(n_layers_unitcell);
	RowVectorXd L(n_layers_unitcell), Ld(n_layers_unitcell);

	urd << 1.0, 1.0;
	erd << 2.1316 - 0.0001i, 4.84 - 0.0007i;

	double d1 = 205.41 * nanometers;
	double d2 = 126.13 * nanometers;
	double df = 280.03 * nanometers;

	L  << d1, d2;
	Ld << df, d2;

	double ff = 0.5;
	double rf = (1-ff)/2;
	double gh = 70 *nanometers;

	double Lx = 336 * nanometers;	
	double Ly = 336 * nanometers;	

	// grating layer
	Layer grating((GratingParam){3,{"Substrate","SiO2","Substrate"},{rf*Lx,ff*Lx,rf*Lx},{Ly,Ly,Ly},{gh,gh,gh},{er1,erd(0),er1},{ur1,urd(0),ur2}});
	grating.setEpsilonAndMu(grid.Nx,grid.Ny);
	//grating.showLayerProperties();
	//grating.showEpsilon();

	// SiO2 layer
	Layer SiO2((LayerParam){1,"SiO2",Lx,Ly,L(0),erd(0),urd(0)});
	SiO2.setEpsilonAndMu(grid.Nx,grid.Ny);
	//SiO2.showLayerProperties();
	//SiO2.showEpsilon();

	// TiO2 layer
	Layer TiO2((LayerParam){1,"TiO2",Lx,Ly,L(1),erd(1),urd(0)});
	TiO2.setEpsilonAndMu(grid.Nx,grid.Ny);
	//TiO2.showLayerProperties();
	//TiO2.showEpsilon();

	// permeability matrix
	// I am gonna use this for now.
	// Need to set a function setMu() and add properties in Layer class
	// for permeability of layers.


	grating.convMat(pq.P, pq.Q);
	grating.showEpConvmat();
	//cout << grating_convmat << endl;
	SiO2.convMat(pq.P, pq.Q);
	//cout << sio2_convmat << endl;
	TiO2.convMat(pq.P, pq.Q);
	//cout << tio2_convmat << endl;

	//cout << grating_convmat.rows() << " " << grating_convmat.cols() << endl;



	return 0;
}
