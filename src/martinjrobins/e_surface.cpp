
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>
#include <math.h> 
#include "sinusoidal_voltammetry.hpp"
#include <iostream>
#include <exception>


struct e_surface_fun {
    double E,dE;
    double Cdl,CdlE,CdlE2,CdlE3;
    double E0;
    double Ru;
    double k0;
    double alpha;
    double In0,u1n0;
    double dt;

    double exp11,exp12;
    double dexp11,dexp12;
    double u1n1;
    double du1n1;
    double Cdlp;

    e_surface_fun ( 
                    const double E,
                    const double dE,
                    const double Cdl,
                    const double CdlE,
                    const double CdlE2,
                    const double CdlE3,
                    const double E0,
                    const double Ru,
                    const double k0,
                    const double alpha,
                    const double In0,
                    const double u1n0,
                    const double dt
                    
                    ) : 
        E(E),dE(dE),Cdl(Cdl),CdlE(CdlE),CdlE2(CdlE2),CdlE3(CdlE3),E0(E0),Ru(Ru),k0(k0),alpha(alpha),In0(In0),u1n0(u1n0),dt(dt) { }

    boost::math::tuple<double,double> operator()(const double In1) {
        update_temporaries(In1);
        return boost::math::make_tuple(residual(In1),residual_gradient(In1));
    }

    double residual(const double In1) const {
        return Cdlp*(dt*dE-Ru*(In1-In0)) - dt*In1 + (u1n1-u1n0);
    }
    double residual_gradient(const double In1) const {
        return -Cdlp*Ru - dt + du1n1;
    }

    void update_temporaries(const double In1) {
        const double Ereduced = E - Ru*In1;
        const double Ereduced2 = pow(Ereduced,2);
        const double Ereduced3 = Ereduced*Ereduced2;
        const double expval1 = Ereduced - E0;
        exp11 = std::exp((1.0-alpha)*expval1);
        exp12 = std::exp(-alpha*expval1);

        dexp11 = -Ru*(1-alpha)*exp11;
        dexp12 = Ru*alpha*exp11;

        const double u1n1_top = dt*k0*exp11 - u1n0;
        const double du1n1_top = dt*k0*dexp11;
        const double denom = (dt*k0*exp11 + dt*k0*exp12 + 1);
        const double ddenom = dt*k0*(dexp11 + dexp12);
        const double tmp = 1.0/denom;
        const double tmp2 = pow(tmp,2);
        u1n1 = u1n1_top*tmp;
        du1n1 = -(u1n1_top*ddenom + du1n1_top*denom)*tmp2;

        Cdlp = Cdl*(1.0 + CdlE*Ereduced + CdlE2*Ereduced2 + CdlE3*Ereduced3);
    }
};

void e_surface(map& params, vector& Itot, vector& t) {
    const double k0 = get(params,std::string("k0"),35.0);
    const double alpha = get(params,std::string("alpha"),0.5);
    const double E0 = get(params,std::string("E0"),0.25);
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
    std::cout << "Running e_surface with parameters:"<<std::endl;
    std::cout << "\tk0 = "<<k0<<std::endl;
    std::cout << "\talpha = "<<alpha<<std::endl;
    std::cout << "\tE0 = "<<E0<<std::endl;
    std::cout << "\tRu = "<<Ru<<std::endl;
    std::cout << "\tCdl = "<<Cdl<<std::endl;
    std::cout << "\tCdlE = "<<CdlE<<std::endl;
    std::cout << "\tCdlE2 = "<<CdlE2<<std::endl;
    std::cout << "\tCdlE3 = "<<CdlE3<<std::endl;
    std::cout << "\tEstart = "<<Estart<<std::endl;
    std::cout << "\tEreverse = "<<Ereverse<<std::endl;
    std::cout << "\tomega = "<<omega<<std::endl;
    std::cout << "\tphase = "<<phase<<std::endl;
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
    double u1n0;
    double t1 = 0;
    Itot0 = 0;
    Itot1 = 0;
    u1n0 = 1.0;
    const double Itot_bound = std::max(10*Cdl*dE*omega/Nt,1.0);
    for (int n_out = 0; n_out < t.size(); n_out++) {
        while (t1 < t[n_out]) {
            Itot0 = Itot1;
            const double E = Eeq(t1+dt);
            const double dE = Eeq.ddt(t1+0.5*dt);
            e_surface_fun bc(E,dE,Cdl,CdlE,CdlE2,CdlE3,E0,Ru,k0,alpha,Itot0,u1n0,dt);

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
            t1 += dt;
        }
        Itot[n_out] = (Itot1-Itot0)*(t[n_out]-t1+dt)/dt + Itot0;
    }
}
