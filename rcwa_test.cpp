#include "RCWA.h"

int main()
{
	//****************************** SOURCE PARAMETERS *****************************//
	// wavelength range 
	int    wl_num = 1;
	double wl_left = 632.8 * nanometers;
	double wl_right = 632.8 * nanometers;

	// theta range
	int    th_num = 5;
	double th_left = 0.0 * degrees;
	double th_right = 90 * degrees;

	// phi range
	int    ph_num = 1;
	double ph_left = 0 * degrees;
	double ph_right = 0 * degrees;

	VectorXd theta, wavelength, phi; 
	theta = VectorXd::LinSpaced(Sequential,th_num,th_left,th_right).transpose();
	wavelength = VectorXd::LinSpaced(Sequential,wl_num,wl_left,wl_right).transpose();
	phi = VectorXd::LinSpaced(Sequential,ph_num,ph_left,ph_right).transpose();
	//******************************************************************************//

	//****************************** DEVICE PARAMETERS ****************************//
	complex<double> ur1 = 1.0;  // permeability in reflection region
	complex<double> er1 = 1.77; // permittivity in reflection region
	complex<double> ur2 = 1.0;  // permeability in transmission region
	complex<double> er2 = 2.5;  // permittivity in transmission region

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
	//******************************************************************************//

	//****************************** RCWA OBJECT ***********************************//
	RCWA rcwa("rw", 2*b_bilayers, (PQ){5,5}, (Grid){512,512}, Lx, Ly);
	rcwa.allocate();

	//****************************** CREATE DEVICE *********************************//
	rcwa.setLayer((LayerParam){1,"in",Lx,Ly,0,er1,ur1});
	rcwa.setLayer((LayerParam){1,"out",Lx,Ly,0,er2,ur2});
	rcwa.addGratingLayer(1,(GratingParam){3,{"Substrate","SiO2","Substrate"},
											{rf*Lx,ff*Lx,rf*Lx},{Ly,Ly,Ly},{gh,gh,gh},{er1,erd(0),er1},
											{ur1,urd(0),ur2}});
	rcwa.addPlaneLayer(2,(LayerParam){1,"SiO2",Lx,Ly,Ld(0),erd(0),urd(0)});
	rcwa.addPlaneLayer(3,(LayerParam){1,"TiO2",Lx,Ly,Ld(1),erd(1),urd(1)});
	rcwa.addPlaneLayer(4,(LayerParam){1,"SiO2",Lx,Ly,L(0),erd(0),urd(0)});
	rcwa.addPlaneLayer(5,(LayerParam){1,"TiO2",Lx,Ly,L(1),erd(1),urd(1)});
	rcwa.addPlaneLayer(6,(LayerParam){1,"SiO2",Lx,Ly,L(0),erd(0),urd(0)});
	rcwa.addPlaneLayer(7,(LayerParam){1,"TiO2",Lx,Ly,L(1),erd(1),urd(1)});
	rcwa.addPlaneLayer(8,(LayerParam){1,"SiO2",Lx,Ly,L(0),erd(0),urd(0)});
	rcwa.addPlaneLayer(9,(LayerParam){1,"TiO2",Lx,Ly,L(1),erd(1),urd(1)});
	rcwa.addPlaneLayer(10,(LayerParam){1,"SiO2",Lx,Ly,L(0),erd(0),urd(0)});
	rcwa.addPlaneLayer(11,(LayerParam){1,"TiO2",Lx,Ly,L(1),erd(1),urd(1)});
	rcwa.addPlaneLayer(12,(LayerParam){1,"SiO2",Lx,Ly,L(0),erd(0),urd(0)});
	rcwa.addPlaneLayer(13,(LayerParam){1,"TiO2",Lx,Ly,L(1),erd(1),urd(1)});
	rcwa.addPlaneLayer(14,(LayerParam){1,"SiO2",Lx,Ly,L(0),erd(0),urd(0)});
	rcwa.addPlaneLayer(15,(LayerParam){1,"TiO2",Lx,Ly,L(1),erd(1),urd(1)});
	rcwa.addPlaneLayer(16,(LayerParam){1,"SiO2",Lx,Ly,L(0),erd(0),urd(0)});
	//******************************************************************************//


	for (int i=0; i<theta.size(); i++)
	{
		rcwa.applyLightSource("TE",wavelength(0),theta(i),phi(0));
		rcwa.run();
		//cout << theta(i)/degrees << " " << rcwa.getReflection() << ";" << endl;
	}

	return 0;
}
