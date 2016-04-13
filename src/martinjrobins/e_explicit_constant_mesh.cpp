
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include "sinusoidal_voltammetry.hpp"
#include <iostream>
#include <exception>



void e_explicit_constant_mesh(map& params, vector& Itot, vector& t) {
    const double k0 = get(params,std::string("k0"),35.0);
    const double alpha = get(params,std::string("alpha"),0.5);
    const double Cdl = get(params,std::string("Cdl"),0.0037);
    const double Ru = get(params,std::string("Ru"),2.74);
    const double E0 = get(params,std::string("E0"),0.0);
    const double dE = get(params,std::string("dE"),0.1);
    const int Nx = get(params,std::string("Nx"),300.0);
    const double Estart = get(params,std::string("Estart"),-10.0);
    const double Ereverse = get(params,std::string("Ereverse"),10.0);
    const double pi = boost::math::constants::pi<double>();
    const double omega = get(params,std::string("omega"),2*pi);
    const double phase = get(params,std::string("phase"),0.0);

    std::cout << "Running ec_impexp_thomas with parameters:"<<std::endl;
    std::cout << "\tk0 = "<<k0<<std::endl;
    std::cout << "\talpha = "<<alpha<<std::endl;
    std::cout << "\tCdl = "<<Cdl<<std::endl;
    std::cout << "\tRu = "<<Ru<<std::endl;
    std::cout << "\tE0 = "<<E0<<std::endl;
    std::cout << "\tdE = "<<dE<<std::endl;
    std::cout << "\tNx = "<<Nx<<std::endl;
    std::cout << "\tEstart = "<<Estart<<std::endl;
    std::cout << "\tEreverse = "<<Ereverse<<std::endl;
    std::cout << "\tomega = "<<omega<<std::endl;
    std::cout << "\tphase = "<<phase<<std::endl;
    std::cout << "\tdE= "<<dE<<std::endl;


    //set up spatial mesh
    const double Xmax = 20;
    const double dx = Xmax/Nx;

    //set up temporal mesh
    const int Nt = 20001;
    const double Tmax = std::abs(Ereverse-Estart)*2;
    const double dt = Tmax/Nt;
    Itot.resize(Nt,0);
    t.resize(Nt);
    for (int i=0; i<Nt; i++) {
        t[i] = i*dt;
    }

    //t=linspace(0,20,Nt);
    //dt=t(2);
    const double mu = dt/pow(dx,2);

    Efun E(Estart,Ereverse,dE,omega,phase,dt);

    const double a = mu;
    const double b = 1+2*mu;
    const double c = mu;
    vector d(Nx);
    vector e(Nx,0);
    vector f(Nx,0);
    f[Nx-1] = 1;

    for (int i=Nx-2; i >= 0; i--) {
        e[i]=a/(b-c*e[i+1]);
    }
    
    vector U(Nx,1);
    //vector V(Nx,0);
    //Itot(1)=Cdl+Cdl*dE*omega/(1+(omega*Ru*Cdl)^2);
    
    for (int n=0; n<Nt-1; n++) {
        for (int i=0; i<Nx; i++) {
            d[i] = U[i];
        }
        for (int i=Nx-2; i>=0; i--) {
            f[i]=(d[i+1]+f[i+1]*c)/(b-c*e[i+1]);
        }
        const double b0=1+dx*k0*(std::exp((1-alpha)*(E[n]-E0-Itot[n]*Ru))+std::exp(-alpha*(E[n]-E0-Itot[n]*Ru)));
        const double d0=dx*k0*exp(-alpha*(E[n]-E0-Itot[n]*Ru));
        U[0]=(d0+f[0])/(b0-e[0]);

        for (int i=0; i<Nx-1; i++) {
            U[i+1]=f[i]+e[i]*U[i];
        }
        Itot[n+1]=(dt*(U[2]-U[1])/dx+Cdl*(E[n+1]-E[n])+Ru*Cdl*Itot[n]) /
                        (dt+Ru*Cdl);
    }
}


