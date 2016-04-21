#ifndef SINUSOIDAL_VOLTAMMETRY_HPP
#define SINUSOIDAL_VOLTAMMETRY_HPP

#include <vector>
#include <map>
#include <string>

#include "utilities.hpp"

typedef std::vector<double> vector; 
typedef std::map<std::string,double> map; 

void e_solution(map& params, vector& Itot, vector& t);
void seq_electron_transfer(map& params, vector& Itot, vector& t);
void e_surface(map& params, vector& Itot, vector& t);

#endif //SINUSOIDAL_VOLTAMMETRY_HPP
