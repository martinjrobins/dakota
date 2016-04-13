
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>
#include <math.h> 
#include "sinusoidal_voltammetry.hpp"
#include <iostream>
#include <exception>
#include <boost/numeric/odeint.hpp>

//typedef std::vector<double> state_type;
typedef boost::numeric::ublas::vector< double > state_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

struct seq_elec_fun {
    Efun E;
    double Cdl,CdlE,CdlE2,CdlE3;
    double E01,E02;
    double Ru;
    double k01,k02;
    double alpha1,alpha2;

    seq_elec_fun ( const Efun E,
                    const double Cdl,
                    const double CdlE,
                    const double CdlE2,
                    const double CdlE3,
                    const double E01,
                    const double E02,
                    const double Ru,
                    const double k01,
                    const double k02,
                    const double alpha1,
                    const double alpha2) : 
        E(E),Cdl(Cdl),CdlE(CdlE),CdlE2(CdlE2),CdlE3(CdlE3),E01(E01),E02(E02),Ru(Ru),k01(k01),k02(k02),alpha1(alpha1),alpha2(alpha2) { }

    void operator() ( const state_type &x , state_type &dxdt , const double t ) {

        const double Et = E(t);
        const double dEdt = E.ddt(t);
    
        const double E2 = Et*Et;
        const double E3 = E2*Et;
        const double Cdlp = Cdl*(1.0 + CdlE*Et + CdlE2*E2 + CdlE3*E3);
        const double expval1 = Et - E01 - Ru*x[2];
        const double expval2 = Et - E02 - Ru*x[2];

        dxdt[0] = k01*((1.0-x[0]-x[1])*std::exp((1.0-alpha1)*expval1) - x[0]*std::exp(-alpha1*expval1));
        dxdt[1] = k02*((1.0-x[0]-x[1])*std::exp(-alpha2*expval2) - x[1]*std::exp((1.0-alpha2)*expval2));
        dxdt[2] = (1.0/(Cdlp*Ru))*(Cdlp*dEdt - x[2] + dxdt[0] - dxdt[1]);
    }

    void init(state_type &x) {
        x.resize(3);
        const double Et = E(0);
        const double dEdt = E.ddt(0);
    
        const double E2 = Et*Et;
        const double E3 = E2*Et;
        const double Cdlp = Cdl*(1.0 + CdlE*Et + CdlE2*E2 + CdlE3*E3);

        x[0] = 1.0;
        x[1] = 0.0;
        x[2] = Cdlp*dEdt;
    }
};

struct seq_elec_fun_jacobi {
    Efun E;
    double Cdl,CdlE,CdlE2,CdlE3;
    double E01,E02;
    double Ru;
    double k01,k02;
    double alpha1,alpha2;

    seq_elec_fun_jacobi ( const Efun E,
                    const double Cdl,
                    const double CdlE,
                    const double CdlE2,
                    const double CdlE3,
                    const double E01,
                    const double E02,
                    const double Ru,
                    const double k01,
                    const double k02,
                    const double alpha1,
                    const double alpha2) : 
        E(E),Cdl(Cdl),CdlE(CdlE),CdlE2(CdlE2),CdlE3(CdlE3),E01(E01),E02(E02),Ru(Ru),k01(k01),k02(k02),alpha1(alpha1),alpha2(alpha2) { }

void operator()( const state_type &x , matrix_type &J , const double t, state_type &dfdt ) {

        const double Et = E(t);
        const double dEdt = E.ddt(t);
        const double dEdt2 = E.ddt2(t);
    
        const double E2 = Et*Et;
        const double E3 = E2*Et;
        const double Cdlp = Cdl*(1.0 + CdlE*Et + CdlE2*E2 + CdlE3*E3);
        const double expval1 = Et - E01 - Ru*x[2];
        const double expval2 = Et - E02 - Ru*x[2];
        const double exp11 = std::exp((1.0-alpha1)*expval1);
        const double exp12 = std::exp(-alpha1*expval1);
        const double exp21 = std::exp(-alpha2*expval2);
        const double exp22 = std::exp((1.0-alpha2)*expval2);

        J(0,0) = k01*(-exp11-exp12);
        J(0,1) = k01*(-exp11);
        J(0,2) = k01*((1-x[0]-x[1])*(Ru*(alpha1-1))*exp11 - x[0]*alpha1*Ru*exp12);

        J(1,0) = k02*(-exp21);
        J(1,1) = k02*(-exp21-exp22);
        J(1,2) = k02*((1-x[0]-x[1])*(Ru*alpha2)*exp21 + x[1]*(1-alpha2)*Ru*exp22);

        J(2,0) = (1.0/(Cdlp*Ru))*(J(0,0)+J(1,0));
        J(2,1) = -(1.0/(Cdlp*Ru))*(J(0,1)+J(1,1));
        J(2,2) = (1.0/(Cdlp*Ru))*(-1.0 + J(0,2) - J(1,2));

        dfdt[0] = k01*dEdt*((1-x[0]-x[1])*(1-alpha1)*exp11 + x[0]*alpha1*exp12);
        dfdt[1] = k02*dEdt*(-(1-x[0]-x[1])*alpha2*exp21 - x[1]*(1-alpha2)*exp22);
        dfdt[2] = (1.0/(Cdlp*Ru))*(Cdlp*dEdt2 + dfdt[0] - dfdt[1]);
    }
};

struct write_output 
{
    vector::iterator i;

    write_output( vector::iterator &i)
    : i(i) { }

    void operator()( const state_type &x , double t ) {
        *i = x[2];
        i++;
    }
};

void seq_electron_transfer2(map& params, vector& Itot, vector& t) {
    const double k01 = get(params,std::string("k01"),35.0);
    const double k02 = get(params,std::string("k02"),65.0);
    const double alpha1 = get(params,std::string("alpha1"),0.5);
    const double alpha2 = get(params,std::string("alpha2"),0.5);
    const double E01 = get(params,std::string("E01"),0.25);
    const double E02 = get(params,std::string("E02"),-0.25);
    const double Ru = get(params,std::string("Ru"),2.74);
    const double Cdl = get(params,std::string("Cdl"),0.0037);
    const double CdlE = get(params,std::string("CdlE"),0.0);
    const double CdlE2 = get(params,std::string("CdlE2"),0.0);
    const double CdlE3 = get(params,std::string("CdlE3"),0.0);
    const double Estart = get(params,std::string("Estart"),-10.0);
    const double Ereverse = get(params,std::string("Ereverse"),10.0);
    const int Nt = get(params,std::string("Nt"),200.0);

    const double pi = boost::math::constants::pi<double>();
    const double omega = get(params,std::string("omega"),2*pi);
    const double phase = get(params,std::string("phase"),0.0);
    const double dE = get(params,std::string("dE"),0.1);
    const double reverse = 0;
    

#ifndef NDEBUG
    std::cout << "Running seq_electron_transfer with parameters:"<<std::endl;
    std::cout << "\tk01 = "<<k01<<std::endl;
    std::cout << "\tk02 = "<<k02<<std::endl;
    std::cout << "\talpha1 = "<<alpha1<<std::endl;
    std::cout << "\talpha2 = "<<alpha2<<std::endl;
    std::cout << "\tE01 = "<<E01<<std::endl;
    std::cout << "\tE02 = "<<E02<<std::endl;
    std::cout << "\tRu = "<<Ru<<std::endl;
    std::cout << "\tCdl = "<<Cdl<<std::endl;
    std::cout << "\tCdlE = "<<CdlE<<std::endl;
    std::cout << "\tCdlE2 = "<<CdlE2<<std::endl;
    std::cout << "\tCdlE3 = "<<CdlE3<<std::endl;
    std::cout << "\tEstart = "<<Estart<<std::endl;
    std::cout << "\tEreverse = "<<Ereverse<<std::endl;
    std::cout << "\tomega = "<<omega<<std::endl;
    std::cout << "\tdE= "<<dE<<std::endl;
    std::cout << "\tNt= "<<Nt<<std::endl;
#endif
    if (Ereverse < Estart) throw std::runtime_error("Ereverse must be greater than Estart");

    //set up temporal mesh
    const double dt = (1.0/Nt)*2*pi/omega;
    if (t.size()==0) {
        const double Tmax = std::abs(Ereverse-Estart)*2;
        const int Nt = Tmax/dt;
        std::cout << "\tNt= "<<Nt<<std::endl;
        Itot.resize(Nt,0);
        t.resize(Nt);
        for (int i=0; i<Nt; i++) {
            t[i] = i*dt;
        }
    } else {
#ifndef NDEBUG
        std::cout << "\thave "<<t.size()<<" samples from "<<t[0]<<" to "<<t[t.size()-1]<<std::endl;
#endif
        Itot.resize(t.size(),0);
    }
    

    Efun E(Estart,Ereverse,dE,omega,phase,dt);
    seq_elec_fun dxdt(E,Cdl,CdlE,CdlE2,CdlE3,E01,E02,Ru,k01,k02,alpha1,alpha2);
    seq_elec_fun_jacobi dxdt_jacobi(E,Cdl,CdlE,CdlE2,CdlE3,E01,E02,Ru,k01,k02,alpha1,alpha2);

    using namespace boost::numeric::odeint;

    state_type x;
    dxdt.init(x);

    vector::iterator i = Itot.begin();
    size_t num_of_steps = integrate_times( 
            make_dense_output< rosenbrock4< double > >( 1.0e-4 , 1.0e-4 ) ,
            std::make_pair( dxdt , dxdt_jacobi ) ,
            x , t.begin() , t.end() , dt ,
            write_output(i)
            );
}
