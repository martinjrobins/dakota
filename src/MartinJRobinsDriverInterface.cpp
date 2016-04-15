/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:        MartinJRobinsDriverInterface
//- Description:  Class implementation
//- Owner:        Mike Eldred, Brian Adams

#include "MartinJRobinsDriverInterface.hpp"
#include "ParallelLibrary.hpp"
#include "DataMethod.hpp"  // for output levels
//#include <unistd.h> // for sleep(int)
#ifdef DAKOTA_MODELCENTER
#include "PHXCppApi.h"
#endif
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/assign.hpp>
#include <vector>
#include "Teuchos_SerialDenseHelpers.hpp"
#include "NonDLHSSampling.hpp"
#include "spectral_diffusion.hpp"

namespace Dakota {


MartinJRobinsDriverInterface::MartinJRobinsDriverInterface(const ProblemDescDB& problem_db)
  : DirectApplicInterface(problem_db)
{
  // register this class' analysis driver types with the string to enum map
  // at the base class
  driverTypeMap["e_surface"]              = E_SURFACE;
  driverTypeMap["e_solution"]             = E_SOLUTION;
  driverTypeMap["seq_electron_transfer"]  = SEQ_ELECTRON_TRANSFER;

  // convert strings to enums for analysisDriverTypes, iFilterType, oFilterType
  analysisDriverTypes.resize(numAnalysisDrivers);
  std::map<String, driver_t>::iterator sd_iter;
  for (size_t i=0; i<numAnalysisDrivers; ++i) {
    sd_iter = driverTypeMap.find(analysisDrivers[i]);//toLower(Drivers[i]));
    if (sd_iter == driverTypeMap.end()) {
      if (outputLevel > NORMAL_OUTPUT)
	Cerr << "Warning: analysis_driver \"" << analysisDrivers[i] << "\" not "
	     << "available at construct time in MartinJRobinsDriverInterface.\n       "
	     << "  Subsequent interface plug-in may resolve." << std::endl;
      analysisDriverTypes[i] = NO_DRIVER;
    }
    else
      analysisDriverTypes[i] = sd_iter->second;
  }

  sd_iter = driverTypeMap.find(iFilterName);  //toLower(iFilterName));
  if (sd_iter == driverTypeMap.end()) {
    if (outputLevel > NORMAL_OUTPUT)
      Cerr << "Warning: input filter \"" << iFilterName << "\" not available at"
	   << " construct time in MartinJRobinsDriverInterface.\n         Subsequent "
	   << "interface plug-in may resolve." << std::endl;
    iFilterType = NO_DRIVER;
  }
  else
    iFilterType = sd_iter->second;

  sd_iter = driverTypeMap.find(oFilterName);  //toLower(oFilterName));
  if (sd_iter == driverTypeMap.end()) {
    if (outputLevel > NORMAL_OUTPUT)
      Cerr << "Warning: output filter \"" << oFilterName << "\" not available "
	   << "at construct time in MartinJRobinsDriverInterface.\n         Subsequent"
	   << " interface plug-in may resolve." << std::endl;
    oFilterType = NO_DRIVER;
  }
  else
    oFilterType = sd_iter->second;

  // define localDataView from analysisDriverTypes,
  // overriding any base class constructor setting
  localDataView = 0;
  for (size_t i=0; i<numAnalysisDrivers; ++i)
    switch (analysisDriverTypes[i]) {
    case E_SURFACE: case E_SOLUTION:
    case SEQ_ELECTRON_TRANSFER: 
      localDataView |= VARIABLES_VECTOR; break;
    }

  // define varTypeMap for analysis drivers based on xCM/XDM maps
  if (localDataView & VARIABLES_MAP) {
      varTypeMap["k0"] = VAR_k0; varTypeMap["alpha"] = VAR_alpha;
      varTypeMap["E0"] = VAR_E0; varTypeMap["Cdl"]  = VAR_Cdl;
  }
}


MartinJRobinsDriverInterface::~MartinJRobinsDriverInterface()
{
  // No tear-down required for now
}


/** Derived map to evaluate a particular built-in test analysis function */
int MartinJRobinsDriverInterface::derived_map_ac(const String& ac_name)
{
  // NOTE: a Factory pattern might be appropriate in the future to manage the
  // conditional presence of linked subroutines in MartinJRobinsDriverInterface.

#ifdef MPI_DEBUG
    Cout << "analysis server " << analysisServerId << " invoking " << ac_name
         << " within MartinJRobinsDriverInterface." << std::endl;
#endif // MPI_DEBUG
  int fail_code = 0;
  std::map<String, driver_t>::iterator sd_iter = driverTypeMap.find(ac_name);
  driver_t ac_type
    = (sd_iter!=driverTypeMap.end()) ? sd_iter->second : NO_DRIVER;
  switch (ac_type) {
  case E_SURFACE:
    fail_code = e_surface(); break;
  case E_SOLUTION:
    fail_code = e_solution(); break;
  case SEQ_ELECTRON_TRANSFER:
    fail_code = seq_electron_transfer(); break;
  default: {
    Cerr << "Error: analysis_driver '" << ac_name << "' is not available in "
	 << "the direct interface." << std::endl;
    abort_handler(INTERFACE_ERROR);
  }
  }

  // Failure capturing
  if (fail_code) {
    std::string err_msg("Error evaluating direct analysis_driver ");
    err_msg += ac_name;
    throw FunctionEvalFailure(err_msg);
  }

  return 0;
}


// -----------------------------------------
// Begin direct interfaces to test functions
// -----------------------------------------

int MartinJRobinsDriverInterface::e_surface(){

  if (multiProcAnalysisFlag) {
    Cerr << "Error: damped oscillator direct fn does not support "
	 << "multiprocessor analyses." << std::endl;
    abort_handler(-1);
  }
  if (numVars < 1 || numVars > 6 || numADIV || numADRV) {
    Cerr << "Error: Bad variable types in e_surface direct fn."
	 << std::endl;
    abort_handler(INTERFACE_ERROR);
  }
  if (numFns < 1) {
    Cerr << "Error: Bad number of functions in e_surface direct fn."
	 << std::endl;
    abort_handler(INTERFACE_ERROR);
  }
  if (hessFlag || gradFlag) {
    Cerr << "Error: Gradients and Hessians not supported in e_surface"
	 << "direct fn." << std::endl;
    abort_handler(INTERFACE_ERROR);
  }

  Real initial_time = 0., final_time = 20., //delta_t = 0.3;
    delta_t = (final_time - initial_time) / numFns;
  //int num_time_steps = numFns;

  Real pi = 4.0 * std::atan( 1.0 );

  Real k0 = xC[0], alpha = 0.5, E0 = 0.0, Cdl = 0.0037, Ru = 2.74;
  if ( numVars >= 2 ) alpha  = xC[1];
  if ( numVars >= 3 ) E0  = xC[2];
  if ( numVars >= 4 ) Cdl  = xC[3];
  if ( numVars >= 5 ) Ru = xC[4];

  std::map<std::string,double> params;
  params["k0"] = k0;
  params["alpha"] = alpha;
  params["Cdl"] = Cdl;
  params["Ru"] = Ru;
  params["E0"] = E0;
  // params for e28q 9 b current_cv_current
  params["dE"] = 5.838264;
  params["Estart"] = -33.083494;
  params["Ereverse"] = -3.892176;
  params["omega"] = 51.764579;
  std::vector<double> Itot,t;

  e_surface(params,Itot,t);

  // response at initial time isn't included as it isn't a fn of some params.
  Real time = initial_time, y_stead, y_trans;
  for (size_t i=0; i<numFns; ++i) {
    time += delta_t;
    if (directFnASV[i] & 1) {
      // Steady state solution (y_stead) for rhs = 0
      y_stead = F * std::sin( w*time + phi ) / zeta;
      // Compute transient (y_trans) component of solution
      y_trans = std::exp( -g*time )*( A1 * std::cos( wd*time ) +
				      A2 * std::sin( wd*time ) );
      fnVals[i] = y_stead + y_trans;
    }
  }

  return 0; // no failure
}
}  // namespace Dakota
