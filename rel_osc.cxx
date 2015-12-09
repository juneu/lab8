// Lorenz.cxx
// Runge Kutta 4th order
// GL, 4.12.2015
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------

//Vorlage fürs Programm von der Lösung von hw6
void f(double* const y0, const double x);
void RKstep(double* const yn, const double* const y0, const double x, const double dx,
            double* k1, double* k2, double* k3, double* k4);
//--------------------
using namespace std;
//--------------------
int main(void)
{
        
	ofstream out("solution");
        const int dim = 2; // Dimension 2
	double dx = 0.1,x=0;
	const double L = 20;
        double y0[dim]; //Anfangsbedingung
  
        double k1[dim], k2[dim], k3[dim], k4[dim]; //damit ks auch in der Main funktion existieren
	double yn[dim];

	
	for( double p0=0.1; p0<5; p0+=0.1){
	x=0;
	y0[0]=p0;
	y0[1]=0;
        //out << x << "\t" << y0[0] << "\t" << y0[1] << "\t" << endl;
	while(x<=L)
	{
		x += dx;
		RKstep(yn, y0, x, dx, k1,k2,k3,k4);
		
		if ((y0[1]>0) && (yn[1]<0)) break; // es müssen beide Bedingungen zutreffen ;   
		  
                for(int i=0; i<dim; i++) y0[i] = yn[i];  // beide gerade berechnete Werte werden in y0[i] gespeichert
	//	out << x << "\t" << y0[0] << "\t" << y0[1] << "\t" << endl;
	}
	
	double Genauigkeit= 0.00001; 
	double thetaM= 0.5;
	double ftheta=y0[1];
	double thetaL=0.0;
	double thetaR=1.0;
	

	while(abs(ftheta)>Genauigkeit) // sobald das kleiner ist als 0.00001
	{  
	
	  double b1=thetaM- (3*thetaM*thetaM)/2+(2*thetaM*thetaM*thetaM)/3;   //Gegebene b-Werte aus Aufgabenstellung
	  double b2=thetaM*thetaM-(2*thetaM*thetaM*thetaM)/3;
	  double b4=-(thetaM*thetaM)/2+(2*thetaM*thetaM*thetaM)/3;
	  
	ftheta=y0[1]+dx*(b1*k1[1]+b2*k2[1]+b2*k3[1]+b4*k4[1]);  //Interpolation zwischen y0 und yn mit verschiedenen theta-Werten
								      
	if(ftheta>0)      // hier tastet er das theta ab: wenn theta>0, setze thetaL=thetaM --> Ändern des Intervalls   --> annähren thetaM an Nullstelle              
	  thetaL=thetaM;  
	else thetaR=thetaM;  
	
	thetaM=(thetaL+thetaR)/2;
	}
	
	out << p0 << "\t" <<x+thetaM*dx<<"\t"<< 1/(x+thetaM*dx) << endl; // x+thetaM*dx: gibt uns die Periodenlänge an   ,  1/(x+thetaM*dx) gibt die Frequenz an
	}
	out.close();
	return(0);
}
//-------------------
void RKstep(double* const yn, const double* const y0,
            const double x, const double dx,
	    double* k1, double* k2, double* k3, double* k4
 	  )
{
	const int dim = 2;
	

  for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + 0.5 * dx * k1[i];
  f(k2, x+0.5*dx);

	for(int i=0;i<dim; i++) k3[i] = y0[i] + 0.5 * dx * k2[i];
	f(k3, x+0.5*dx);

  for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * k3[i];
	f(k4,  x+dx);

	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + 1./6.*dx*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
	
	
	
	    
	    
	  
	  
}
//-------------------
// 
void f(double* const y0, const double x)
{
	
	double y[2] = { y0[0], y0[1]};

  y0[0] = y[1];
  y0[1] = -y[0]/(sqrt(y[0]*y[0]+1));
	
}
