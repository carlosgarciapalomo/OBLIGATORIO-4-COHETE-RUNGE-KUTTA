//PROGRAMA OBLIGATORIO NUMERO 4: COHETE. Carlos García Palomo
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>


static const double G = 6.67e-11;//Constante de gravitacion
static const double Mt = 5.9736e24;//Masa de la tierra
static const double Ml = 0.07349e24;//Masa de la luna
static const double dtl = 3.844e8;//Distancia tierra luna
static const double w = 2.6617e-6;//Velocidad angular de la tierra
static const double Rt = 6.378160e6;//Radio de la tierra
static const double Rl = 1.7374e6;//Radio de la luna
static const double delta = G*Mt/(dtl*dtl*dtl);
static const double mu = Ml/Mt;

void actualizak (double f[],double h, double t ,double Y[], double k1[],double k2[],double k3[],double k4[]);
void actualizay (double Y[],double k1[],double k2[],double k3[],double k4[]);
void actualizaf (double f[],double r,double phi, double pr,double pphi,double t);
double rprima (double r,double phi,double t);
double rpunto (double t,double pr);
double phipunto (double t,double pphi, double r);
double prprima (double t, double r, double pphi,double phi,double rp );
double pphiprima (double t, double r,double phi, double rp);
double polaresx (double r, double ang);
double polaresy (double r, double ang);




using namespace std;

//Para que funcione se han dado las condiciones iniciales: phi=0, theta=0.128, v=1.1*Vescape

int main()
{
    int h=60;
    int i=0,j;//Contador de pasos
    int pasos;//Numero de pasos
    double Y[4];//Vector que contiene los 4 valores de las ecuaciones diferenciales ordenados como: r,phi,pr y pphy
    double m=20000; //Masa del satelite
    double f[4];//Venctor que contiene los resultados de las funciones
    double k1[4],k2[4],k3[4],k4[4];//Vectores para las ks
    double v,theta;//Este valor nos facilitara el calculo de las condiciones iniciales
    double t=0.0;
    double hamiltoniano,conservacion;


    ofstream fich, fich2;

    fich.open("Hamiltoniano.txt");
    fich2.open("Tierra-luna.txt");


    //Condiciones iniciales
    Y[0]=Rt/dtl;
    Y[1]=0;
    theta=0.128;
    v=1.1*sqrt(2*G*Mt/Rt);
    pasos=10000;


    Y[2]=(v/dtl)*cos(theta-Y[1]);//Calculo de los momentos conjugados
    Y[3]=Y[0]*(v/dtl)*sin(theta-Y[1]);


    //Implemento el algoritmo
    while (i<pasos)
    {
        actualizak(f,h,t,Y,k1,k2,k3,k4);
        actualizay(Y,k1,k2,k3,k4);
         
         hamiltoniano=pow(m,2)*dtl*dtl*(0.5*pow(Y[2],2)*pow(Y[3]/Y[0],2)-delta*(pow(Y[0],-1)+mu/rprima(Y[0],Y[1],t)));
         conservacion=hamiltoniano-w*Y[3];
         fich<< t <<" "<<conservacion<<endl;
   

        if (i%10==0)
        {
            fich2<<polaresx(Y[0],Y[1])<<" "<<polaresy(Y[0],Y[1])<<" ";
            fich2<<cos(w*t)<<" "<<sin(w*t)<<endl;
            
        }

        t=t+h;

        i++;
    }

    fich.close();
    fich2.close();

    return 0;
}

void actualizak (double f[],double h, double t ,double Y[], double k1[],double k2[],double k3[],double k4[])
{
    int i;

    actualizaf(f,Y[0],Y[1],Y[2],Y[3],t);
    for (i=0;i<4;i++)
    {
        k1[i]=h*f[i];
    }

    actualizaf(f,Y[0]+k1[0]*1.0/2.0,Y[1]+k1[1]*1.0/2.0,Y[2]+k1[2]*1.0/2,Y[3]+k1[3]*1.0/2,t+h/2);
    for (i=0;i<4;i++)
    {
        k2[i]=h*f[i];
    }

    actualizaf(f,Y[0]+k2[0]*1.0/2.0,Y[1]+k2[1]*1.0/2.0,Y[2]+k2[2]*1.0/2.0,Y[3]+k2[3]*1.0/2.0,t+h/2);
    for (i=0;i<4;i++)
    {
        k3[i]=h*f[i];
    }

    actualizaf(f,Y[0]+k3[0],Y[1]+k3[1],Y[2]+k3[2],Y[3]+k3[3],t+h);
    for (i=0;i<4;i++)
    {
        k4[i]=h*f[i];
    }

    return;
}

void actualizaf (double f[],double r,double phi, double pr,double pphi,double t)
{
    double rp;
    rp=rprima(r,phi,t);

    f[0]=rpunto(t,pr);
    f[1]=phipunto(t,pphi,r);
    f[2]=prprima(t,r,pphi,phi,rp);
    f[3]=pphiprima(t,r,phi,rp);
}

void actualizay (double Y[],double k1[],double k2[],double k3[],double k4[])
{
    int i;

    for(i=0;i<4;i++)
        Y[i]=Y[i]+(1.0/6.0)*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
}

double rprima (double r,double phi,double t)
{
    return sqrt(1.0+r*r-2.0*r*cos(phi-w*t));
}

double rpunto (double t,double pr)
{
    return pr;
}

double phipunto (double t,double pphi, double r)
{
    return pphi/(r*r);
}

double prprima (double t, double r, double pphi,double phi,double rp )
{
   return (pphi*pphi/(r*r*r))-delta*(1.0/(r*r)+(mu/(rp*rp*rp))*(r-cos(phi-w*t)));
}

double pphiprima (double t, double r,double phi, double rp)
{
    return -1.0*(delta*mu*r)/(rp*rp*rp)*sin(phi-w*t);
}

double polaresx (double r, double ang)
{
    return r*cos(ang);
}

double polaresy (double r, double ang)
{
    return r*sin(ang);
}

