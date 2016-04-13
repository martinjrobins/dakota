#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <cmath>
#include <map>

template <typename K, typename V>
V get(const  std::map <K,V> & m, const K & key, const V & defval ) {
    typename std::map<K,V>::const_iterator it = m.find( key );
    if ( it == m.end() ) {
        return defval;
    } else {
        return it->second;
    }
}

struct Efun {
    Efun() {};
    Efun(const double Estart, const double Ereverse, const double dE, const double omega, const double phase, const double dt):
        Estart(Estart),Ereverse(Ereverse),dE(dE),omega(omega),phase(phase),dt(dt),
        treverse(std::abs(Estart-Ereverse)),upfirst(Ereverse>Estart) {
        };
    double operator[](const int n) const {
        const double t = n*dt;
        if (t<treverse) {
            return Estart + t + dE*std::sin(omega*t+phase);
        } else {
            return Ereverse - (t-treverse) + dE*std::sin(omega*t+phase);
        }
    }
    double operator()(const double t) const {
         if (t<treverse) {
            return Estart + t + dE*std::sin(omega*t+phase);
        } else {
            return Ereverse - (t-treverse) + dE*std::sin(omega*t+phase);
        }
    }
    double dc(const double t) const {
         if (t<treverse) {
            return Estart + t;
        } else {
            return Ereverse - (t-treverse);
        }
    }
    double ddt(const double t) const {
        if (t>treverse) {
            return -1 + omega*dE*std::cos(omega*t+phase);
        } else {
            return 1 + omega*dE*std::cos(omega*t+phase);
        }
    }
    double ddt2(const double t) const {
        return -std::pow(omega,2)*dE*std::sin(omega*t+phase);
    }
    void set_phase(const double new_phase) { phase = new_phase; }
    double Estart,Ereverse,dE,omega,dt,treverse;
    double phase;
    bool upfirst;
};


#endif
