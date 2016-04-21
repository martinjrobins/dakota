
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>
#include <math.h> 
#include "sinusoidal_voltammetry.hpp"
#include <iostream>
#include <exception>


struct seq_elec_fun {
    double E,dE;
    double Cdl;
    double E01,E02;
    double Ru;
    double k01,k02;
    double alpha1,alpha2;
    double In0,u1n0,u2n0;
    double dt;

    double exp11,exp12,exp21,exp22;
    double dexp11,dexp12,dexp21,dexp22;
    double u1n1,u2n1;
    double du1n1,du2n1;

    seq_elec_fun ( 
                    const double E,
                    const double dE,
                    const double Cdl,
                    const double E01,
                    const double E02,
                    const double Ru,
                    const double k01,
                    const double k02,
                    const double alpha1,
                    const double alpha2,
                    const double In0,
                    const double u1n0,
                    const double u2n0,
                    const double dt
                    
                    ) : 
        E(E),dE(dE),Cdl(Cdl),E01(E01),E02(E02),Ru(Ru),k01(k01),k02(k02),alpha1(alpha1),alpha2(alpha2),In0(In0),u1n0(u1n0),u2n0(u2n0),dt(dt) { }

    boost::math::tuple<double,double> operator()(const double In1) {
        update_temporaries(In1);
        return boost::math::make_tuple(residual(In1),residual_gradient(In1));
    }

    double residual(const double In1) const {
        return Cdl*(dt*dE-Ru*(In1-In0)) - dt*In1 + (u1n1-u1n0) - (u2n1-u2n0);
    }
    double residual_gradient(const double In1) const {
        return -Cdl*Ru - dt + du1n1 - du2n1;
    }

    void update_temporaries(const double In1) {
        const double expval1 = E - E01 - Ru*In1;
        const double expval2 = E - E02 - Ru*In1;
        exp11 = std::exp((1.0-alpha1)*expval1);
        exp12 = std::exp(-alpha1*expval1);
        exp21 = std::exp(-alpha2*expval2);
        exp22 = std::exp((1.0-alpha2)*expval2);

        dexp11 = -Ru*(1-alpha1)*exp11;
        dexp12 = Ru*alpha1*exp11;
        dexp21 = Ru*alpha2*exp21;
        dexp22 = -Ru*(1-alpha2)*exp22;

        const double u1n1_top = dt*k01*(dt*k02*exp21 + u2n0)*exp11 - (dt*k01*exp11 + u1n0)*(dt*k02*exp21 + dt*k02*exp22 + 1);
        const double u2n1_top = dt*k02*(dt*k01*exp11 + u1n0)*exp21 - (dt*k02*exp21 + u2n0)*(dt*k01*exp11 + dt*k01*exp12 + 1);

        const double du1n1_top = dt*(-dt*k01*k02*exp11*dexp22 - dt*k01*k02*exp22*dexp11 + k01*u2n0*dexp11 - k01*dexp11 - k02*u1n0*dexp21 - k02*u1n0*dexp22);
        const double du2n1_top = dt*(-dt*k01*k02*exp12*dexp21 - dt*k01*k02*exp21*dexp12 - k01*u2n0*dexp11 - k01*u2n0*dexp12 + k02*u1n0*dexp21 - k02*dexp21);

        const double denom = pow(dt,2)*k01*k02*exp11*exp21 - (dt*k01*exp11 + dt*k01*exp12 + 1)*(dt*k02*exp21 + dt*k02*exp22 + 1);
        const double ddenom = dt*(dt*k01*k02*exp11*dexp21 + dt*k01*k02*exp21*dexp11 - k01*(dexp11 + dexp12)*(dt*k02*exp21 + dt*k02*exp22 + 1) - k02*(dexp21+dexp22)*(dt*k01*exp11 + dt*k01*exp12 + 1));

        const double tmp = 1.0/denom;
        const double tmp2 = pow(tmp,2);
        u1n1 = u1n1_top*tmp;
        u2n1 = u2n1_top*tmp;
        du1n1 = -(u1n1_top*ddenom + du1n1_top*denom)*tmp2;
        du2n1 = -(u2n1_top*ddenom + du2n1_top*denom)*tmp2;
    }
};

void seq_electron_transfer(map& params, vector& Itot, vector& t) {
    const double k01 = get(params,std::string("k01"),35.0);
    const double k02 = get(params,std::string("k02"),65.0);
    const double alpha1 = get(params,std::string("alpha1"),0.5);
    const double alpha2 = get(params,std::string("alpha2"),0.5);
    const double E01 = get(params,std::string("E01"),0.25);
    const double E02 = get(params,std::string("E02"),-0.25);
    const double Ru = get(params,std::string("Ru"),0.001);
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
    
    const int digits_accuracy = std::numeric_limits<double>::digits*0.5;
    const double max_iterations = 100;

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
    

    Efun Eeq(Estart,Ereverse,dE,omega,phase,dt);

    double Itot0,Itot1;
    double u1n0,u2n0;
    double t1 = 0;
    Itot0 = 0;
    Itot1 = 0;
    u1n0 = 1.0;
    u2n0 = 0.0;
    const double Itot_bound = std::max(10*Cdl*dE*omega/Nt,1.0);
    for (int n_out = 0; n_out < t.size(); n_out++) {
        while (t1 < t[n_out]) {
            Itot0 = Itot1;
            const double E = Eeq(t1+dt);
            const double dE = Eeq.ddt(t1+0.5*dt);
            const double Edc = Eeq.dc(t1+dt);
            const double Cdlp = Cdl*(1.0 + CdlE*E + CdlE2*pow(E,2) + CdlE3*pow(E,3));
            seq_elec_fun bc(E,dE,Cdlp,E01,E02,Ru,k01,k02,alpha1,alpha2,Itot0,u1n0,u2n0,dt);

            boost::uintmax_t max_it = max_iterations;
            Itot1 = boost::math::tools::newton_raphson_iterate(bc, Itot0,Itot0-Itot_bound,Itot0+Itot_bound, digits_accuracy, max_it);
            //Itot1 = boost::math::tools::bisect(bc,Itot0-1.1,Itot0+1.1,tol,max_it).first;
            if (max_it == max_iterations) throw std::runtime_error("non-linear solve for Itot[n+1] failed, max number of iterations reached");

            //std::cout << "residual is "<<bc.residual(Itot1)<<std::endl;
            //std::cout << "max_it "<<max_it<<std::endl;
            //std::cout << "residual gradient is "<<bc.residual_gradient(Itot1)<<std::endl;
            //std::cout << "If is "<<bc.If(Itot1)<<std::endl;
            //std::cout << "If2 is "<<bc.If2(Itot1)<<std::endl;
            bc.update_temporaries(Itot1);
            u1n0 = bc.u1n1;
            u2n0 = bc.u2n1;
            t1 += dt;
        }
        Itot[n_out] = (Itot1-Itot0)*(t[n_out]-t1+dt)/dt + Itot0;
    }
}
