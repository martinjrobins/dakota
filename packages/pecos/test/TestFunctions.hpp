#include "pecos_data_types.hpp"

using namespace Pecos;

void get_genz_coefficients( int num_dims, Real factor, int c_type, RealVector &c, RealVector &w )
{
  c.resize( num_dims );
  w.resize( num_dims );
  switch ( c_type )
    {
    case 0: // Linear decay
      {
	Real csum = 0.0;
	for ( int d = 0; d < num_dims; d++ )
	  {
	    w[d] = 0.0;
	    c[d] = ( (Real)d + 0.5 ) / (Real)num_dims;
	    csum += c[d]; 
	  }
	for ( int d = 0; d < num_dims; d++ )
	  {
	    c[d] *= ( factor / csum );
	  }
	break;
      }
    case 1: // Quadratic decay
      {
	Real csum = 0.0;
	for ( int d = 0; d < num_dims; d++ )
	  {
	    w[d] = 0.0;
	    c[d] = 1.0 / (Real)( ( d + 1 ) * ( d + 1 ) );
	    csum += c[d]; 
	  }
	for ( int d = 0; d < num_dims; d++ )
	  {
	    c[d] *= ( factor / csum );
	  }
	break;
      }
    case 2: // Exponential decay
      {
	Real csum = 0.0;
	for ( int d = 0; d < num_dims; d++ )
	  {
	    w[d] = 0.;
	    c[d] = std::exp( (Real)(d+1)*std::log( 1.e-8 ) / (Real)num_dims );
	    csum += c[d];
	  }
	for ( int d = 0; d < num_dims; d++ )
	  {
	    c[d] *= ( factor / csum );
	  }
	break;
      }
    default:
      throw(std::runtime_error("GetCoefficients() ensure type in [0,1,2]"));
    };
}

double genz(String an_comp, RealVector &xC)
{
  int coeff_type, fn_type;
  Real factor;
  if (an_comp == "os1")
    { coeff_type = 0; fn_type = 0; factor = 4.5; }
  else if (an_comp == "os2")
    { coeff_type = 1; fn_type = 0; factor = 4.5; }
  else if (an_comp == "os3")
    { coeff_type = 2; fn_type = 0; factor = 4.5; }
  else if (an_comp == "cp1")
    { coeff_type = 0; fn_type = 1; factor = .25; }
  else if (an_comp == "cp2")
    { coeff_type = 1; fn_type = 1; factor = .25; }
  else if (an_comp == "cp3")
    { coeff_type = 2; fn_type = 1; factor = .25; }
  else {
    std::cerr << "Error: option for genz() is not yet implemented : "<<an_comp
	 << std::endl;
    std::terminate();
  } 

  RealVector c, w;
  int numVars = xC.length();
  get_genz_coefficients( numVars, factor, coeff_type, c, w );
  Real fnVals, pi = 4.0 * std::atan( 1.0 );

  // **** f:
  switch (fn_type) {
    case 0: {
      fnVals = 2.0 * pi * w[0];
      for ( int d = 0; d < numVars; d++ ){
	fnVals += c[d] * xC[d];
      }
      fnVals = std::cos( fnVals );
      break;
    }
    case 1:{
      fnVals = 1.0;
      for ( int d = 0; d < numVars; d++ )
	{
	  fnVals += c[d]* xC[d];
	}
      fnVals = 1.0 / std::pow( fnVals, (Real)(numVars+1) );
      break;
    }
  }

  return fnVals; // no failure
}
