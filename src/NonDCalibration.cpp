/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	 NonDCalibration
//- Description: Base class for nondeterministic calibration
//- Owner:       Laura Swiler
//- Checked by:
//- Version:

#include <boost/math/special_functions/round.hpp>
#include "NonDCalibration.hpp"
#include "DakotaModel.hpp"
#include "ProblemDescDB.hpp"

static const char rcsId[]="@(#) $Id$";


namespace Dakota {

/** This constructor is called for a standard letter-envelope iterator 
    instantiation.  In this case, set_db_list_nodes has been called and 
    probDescDB can be queried for settings from the method specification. */
NonDCalibration::NonDCalibration(ProblemDescDB& problem_db, Model& model):
  NonD(problem_db, model),
  calibrationData(probDescDB.get_bool("responses.calibration_data") ||
    !probDescDB.get_string("responses.scalar_data_filename").empty()),
  numExperiments(probDescDB.get_sizet("responses.num_experiments")),
  numExpConfigVars(probDescDB.get_sizet("responses.num_config_vars")),
  varianceTypesRead(probDescDB.get_sa("responses.variance_type")),
  expData(problem_db, iteratedModel.current_response().shared_data(), 
	  outputLevel),
  numTotalCalibTerms(model.num_primary_fns()),
  continuousConfigVars(0), discreteIntConfigVars(0), discreteRealConfigVars(0),
  continuousConfigStart(0), discreteIntConfigStart(0),
  discreteRealConfigStart(0)
{ 
  bool found_error = false;
  if (numExpConfigVars > 0) {

    // would need to trap error if all_variables were later allowed

    // available config var counts
    continuousConfigVars =
      probDescDB.get_sizet("variables.continuous_state"); 
    discreteIntConfigVars =
      probDescDB.get_sizet("variables.discrete_state_range") + 
      probDescDB.get_sizet("variables.discrete_state_set_int");
    discreteRealConfigVars =
      probDescDB.get_sizet("variables.discrete_state_set_real"); 

    size_t total_state_config_vars = 
      continuousConfigVars + discreteIntConfigVars + discreteRealConfigVars;

    if (numExpConfigVars > 0 && numExpConfigVars != total_state_config_vars) {
      Cerr << "\nError (NonDCalibration): experimental_config_variables = "
	   << numExpConfigVars << " read from file must equal total state "
	   << "variables = " << total_state_config_vars << std::endl;
      found_error = true;
    }
    else {

      // set indices of the state variables within the all* variable arrays
      bool index_found = false;
      if (continuousConfigVars > 0) {
	index_found = 
	  find_state_index(CONTINUOUS_STATE, 
			   iteratedModel.all_continuous_variable_types(), 
			   "continuous state", continuousConfigStart);
	if (!index_found) found_error = true;
      }

      // this assumes that DISCRETE_STATE_RANGE and DISCRETE_STATE_SET_INT
      // are continguous
      if (discreteIntConfigVars > 0) {
	index_found = 
	  find_state_index(DISCRETE_STATE_RANGE, 
			   iteratedModel.all_discrete_int_variable_types(), 
			   "discrete state range", discreteIntConfigStart);
	if (!index_found) found_error = true;
      }
      
      if (discreteRealConfigVars > 0) {
	index_found = 
	  find_state_index(DISCRETE_STATE_SET_REAL, 
			   iteratedModel.all_discrete_real_variable_types(), 
			   "discrete state set real", discreteRealConfigStart);
	if (!index_found) found_error = true;
      }

    }
    
  } // numExpConfigVars

  if (found_error)
    abort_handler(-1);

  // Read in all of the experimental data, including any x configuration 
  // variables, y observations, and covariance information if available 
  //if (outputLevel > NORMAL_OUTPUT)
  //  Cout << "Read data from file " << calibrationData << '\n';
  if (calibrationData) {
    expData.load_data("NonDCalibration");
    numTotalCalibTerms = expData.num_total_exppoints();
  } else {
    Cout << "No experiment data from files.\nCalibration is assuming the "
	 << "simulation is returning the residuals" << std::endl;
  }
}


NonDCalibration::~NonDCalibration()
{ }


//void NonDCalibration::print_results(std::ostream& s)
//{ Cout << "Posterior sample results " << '\n'; }

void NonDCalibration::
set_configuration_vars(Model& model, const RealVector& config_vars) {

  // config_vars consists of [continuous, discreteInt, discreteReal] in order

  // current index into configuration variables
  size_t config_ind = 0;
  // where to insert in the all variables array
  size_t state_ins = continuousConfigStart;
  for ( ; config_ind < continuousConfigVars; ++config_ind, ++state_ins)
    model.all_continuous_variable(config_vars[config_ind], state_ins);

  // don't reset the config index
  state_ins = discreteIntConfigStart;
  // TODO: read integer config vars as integers instead of rounding
  //       (need file reader supporting partitioned matrices)
  for ( ; config_ind < discreteIntConfigVars; ++config_ind, ++state_ins) {
    int int_config_var = boost::math::iround(config_vars[config_ind]);
    model.all_discrete_int_variable(int_config_var, state_ins);
  }

  // don't reset the config index
  state_ins = discreteRealConfigStart;
  for ( ; config_ind < discreteRealConfigVars; ++config_ind, ++state_ins)
    model.all_discrete_real_variable(config_vars[config_ind], state_ins);

}


bool NonDCalibration::
find_state_index(unsigned short state_type,
		 UShortMultiArrayConstView variable_types,
		 std::string context_message,
		 size_t& start_index)
{
  UShortMultiArray::const_iterator cit
    = std::find(variable_types.begin(), variable_types.end(), state_type);
  if (cit == variable_types.end()) { 
    Cerr << "\nError looking up " << context_message << " state variable index "
	 << "in (NonDCalibration)" << std::endl;
    start_index = 0;
    return(false);
  }
  start_index = std::distance(variable_types.begin(), cit);
  return(true);
}

} // namespace Dakota
