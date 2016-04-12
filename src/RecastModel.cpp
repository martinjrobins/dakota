/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:       RecastModel
//- Description: Implementation code for the RecastModel class
//- Owner:       Mike Eldred
//- Checked by:

#include "dakota_system_defs.hpp"
#include "RecastModel.hpp"

static const char rcsId[]="@(#) $Id: RecastModel.cpp 7029 2010-10-22 00:17:02Z mseldre $";


using namespace std;

namespace Dakota {

// define special values for componentParallelMode
#define SUB_MODEL 2
//#define DEBUG


/** Default recast model constructor.  Requires full definition of the
    transformation.  Parameter vars_comps_totals indicates the number
    of each type of variable {4 types} x {3 domains} in the recast
    variable space.  Note: recast_secondary_offset is the start index
    for equality constraints, typically num nonlinear ineq constraints. */
RecastModel::
RecastModel(const Model& sub_model, const Sizet2DArray& vars_map_indices,
	    const SizetArray& vars_comps_totals, const BitArray& all_relax_di,
	    const BitArray& all_relax_dr, bool nonlinear_vars_mapping,
	    void (*variables_map)      (const Variables& recast_vars,
					Variables& sub_model_vars),
	    void (*set_map)            (const Variables& recast_vars,
					const ActiveSet& recast_set,
					ActiveSet& sub_model_set),
	    const Sizet2DArray& primary_resp_map_indices,
	    const Sizet2DArray& secondary_resp_map_indices,
	    size_t recast_secondary_offset, short recast_resp_order,
	    const BoolDequeArray& nonlinear_resp_mapping,
	    void (*primary_resp_map)   (const Variables& sub_model_vars,
					const Variables& recast_vars,
					const Response& sub_model_response,
					Response& recast_response),
	    void (*secondary_resp_map) (const Variables& sub_model_vars,
					const Variables& recast_vars,
					const Response& sub_model_response,
					Response& recast_response)):
  Model(LightWtBaseConstructor(), sub_model.problem_description_db(),
	sub_model.parallel_library()),
  subModel(sub_model), varsMapIndices(vars_map_indices),
  nonlinearVarsMapping(nonlinear_vars_mapping), variablesMapping(variables_map),
  setMapping(set_map), primaryRespMapIndices(primary_resp_map_indices),
  secondaryRespMapIndices(secondary_resp_map_indices),
  nonlinearRespMapping(nonlinear_resp_mapping),
  primaryRespMapping(primary_resp_map),
  secondaryRespMapping(secondary_resp_map),// inverseMapFlag(false),
  invVarsMapping(NULL), invSetMapping(NULL), invPriRespMapping(NULL),
  invSecRespMapping(NULL)
{
  modelType = "recast"; supportsEstimDerivs = false;

  // synchronize output level and grad/Hess settings with subModel
  initialize_data_from_submodel();

  // recasting of variables
  const Variables& sub_model_vars = subModel.current_variables();
  bool reshape_vars; // only reshape if change in variable type counts
  SharedVariablesData recast_svd;

  // BMA TODO: it's possible vars_comp_totals is resized, but no mapping...

  // variables are not mapped: deep copy of vars to allow independence, but 
  // shallow copy of svd since types/labels/ids can be kept consistent
  if (variablesMapping == NULL) {
    currentVariables = sub_model_vars.copy(); // shared svd
    reshape_vars = false;
  }
  // variables are mapped but not resized: deep copy of vars and svd, since
  // types may change in transformed space
  else {
    const SharedVariablesData& svd = sub_model_vars.shared_data();
    if ( ( vars_comps_totals.empty() ||
	   svd.components_totals()         == vars_comps_totals ) &&
	 ( all_relax_di.empty() ||
	   svd.all_relaxed_discrete_int()  == all_relax_di )      &&
	 ( all_relax_dr.empty() || 
	   svd.all_relaxed_discrete_real() == all_relax_dr ) ) {
      currentVariables = sub_model_vars.copy(true); // independent svd
      reshape_vars = false;
    }
    else { // variables are resized
      recast_svd = SharedVariablesData(sub_model_vars.view(), vars_comps_totals,
				       all_relax_di, all_relax_dr);
      currentVariables = Variables(recast_svd);
      reshape_vars = true;
    }
  }
  // propagate number of active continuous vars to deriv vars
  numDerivVars = currentVariables.cv();

  respMapping = (primaryRespMapping || secondaryRespMapping);
  size_t num_recast_primary_fns = primaryRespMapIndices.size(),
    num_recast_secondary_fns = secondaryRespMapIndices.size(),
    num_recast_fns = num_recast_primary_fns + num_recast_secondary_fns;
  if (nonlinearRespMapping.size() != num_recast_fns) {
    Cerr << "Error: size mismatch in response mapping configuration." << endl;
    abort_handler(-1);
  }

  // recasting of response
  const Response& sub_model_resp = subModel.current_response();
  currentResponse = sub_model_resp.copy();
  if (respMapping) {
    numFns = num_recast_fns;
    bool grad_flag = (recast_resp_order & 2),
         hess_flag = (recast_resp_order & 4),
         sm_grad_flag = !sub_model_resp.function_gradients().empty(),
         sm_hess_flag = !sub_model_resp.function_hessians().empty();
    if ( sub_model_vars.cv()            != numDerivVars ||
	 sub_model_resp.num_functions() != numFns       ||
	 grad_flag != sm_grad_flag || hess_flag != sm_hess_flag )
      currentResponse.reshape(numFns, numDerivVars, grad_flag, hess_flag);
  }
  else
    numFns = currentResponse.num_functions();

  // recasting of constraints
  const Constraints& sub_model_cons = subModel.user_defined_constraints();
  userDefinedConstraints = (reshape_vars) ?
    Constraints(recast_svd) : sub_model_cons.copy();
  if (secondaryRespMapping) {
    // the recast_secondary_offset cannot in general be inferred from the
    // contributing fns in secondaryRespMapIndices (recast constraints may be
    // defined, e.g., with no contributing fns), and must therefore be passed.
    size_t num_recast_nln_ineq = recast_secondary_offset,
      num_recast_nln_eq = num_recast_secondary_fns - num_recast_nln_ineq;
    if ( num_recast_nln_ineq != sub_model_cons.num_nonlinear_ineq_constraints()
      || num_recast_nln_eq   != sub_model_cons.num_nonlinear_eq_constraints() )
      userDefinedConstraints.reshape(num_recast_nln_ineq, num_recast_nln_eq,
        sub_model_cons.num_linear_ineq_constraints(),
        sub_model_cons.num_linear_eq_constraints());
  }
}


/** This alternate constructor defers initialization of the function
    pointers until a separate call to initialize(), and accepts the
    minimum information needed to construct currentVariables,
    currentResponse, and userDefinedConstraints.  The resulting model
    is sufficiently complete for passing to an Iterator.  Parameter
    vars_comps_totals indicates the number of each type of variable {4
    types} x {3 domains} in the recast variable space. Note:
    recast_secondary_offset is the start index for equality
    constraints, typically num nonlinear ineq constraints. */
RecastModel::
RecastModel(const Model& sub_model, //size_t num_deriv_vars,
	    const SizetArray& vars_comps_totals, const BitArray& all_relax_di,
	    const BitArray& all_relax_dr,    size_t num_recast_primary_fns,
	    size_t num_recast_secondary_fns, size_t recast_secondary_offset,
	    short recast_resp_order):
  Model(LightWtBaseConstructor(), sub_model.problem_description_db(),
	sub_model.parallel_library()),
  subModel(sub_model), nonlinearVarsMapping(false), respMapping(false),
  variablesMapping(NULL), setMapping(NULL), primaryRespMapping(NULL),
  secondaryRespMapping(NULL),// inverseMapFlag(false),
  invVarsMapping(NULL), invSetMapping(NULL), invPriRespMapping(NULL),
  invSecRespMapping(NULL)
{
  modelType = "recast"; supportsEstimDerivs = false;

  // synchronize output level and grad/Hess settings with subModel
  initialize_data_from_submodel();

  // recasting of variables
  const Variables& sub_model_vars = subModel.current_variables();
  bool reshape_vars = false; // only reshape if change in variable type counts
  SharedVariablesData recast_svd;
  // variables may be mapped, but the mapping hasn't been provided
  // yet; may need to resize based on vars_comps_totals

  // variables are mapped but not resized: deep copy of vars and svd, since
  // types may change in transformed space
  const SharedVariablesData& svd = sub_model_vars.shared_data();
  if ( ( vars_comps_totals.empty() || 
         svd.components_totals()         == vars_comps_totals ) &&
       ( all_relax_di.empty() || 
         svd.all_relaxed_discrete_int()  == all_relax_di )      &&
       ( all_relax_dr.empty() || 
         svd.all_relaxed_discrete_real() == all_relax_dr ) ) {
    currentVariables = sub_model_vars.copy(true); // independent svd
    reshape_vars = false;
  }
  else { // variables are resized
    recast_svd = SharedVariablesData(sub_model_vars.view(), vars_comps_totals,
                                     all_relax_di, all_relax_dr);
    currentVariables = Variables(recast_svd);
    reshape_vars = true;
  }
  // propagate number of active continuous vars to deriv vars
  numDerivVars = currentVariables.cv();

  // recasting of response
  const Response& sub_model_resp = subModel.current_response();
  currentResponse = sub_model_resp.copy();
  numFns = num_recast_primary_fns + num_recast_secondary_fns;
  bool grad_flag = (recast_resp_order & 2), hess_flag = (recast_resp_order & 4),
       sm_grad_flag = !sub_model_resp.function_gradients().empty(),
       sm_hess_flag = !sub_model_resp.function_hessians().empty();
  if ( sub_model_vars.cv()            != numDerivVars ||
       sub_model_resp.num_functions() != numFns       ||
       sm_grad_flag != grad_flag || sm_hess_flag != hess_flag )
    currentResponse.reshape(numFns, numDerivVars, grad_flag, hess_flag);

  // recasting of constraints
  const Constraints& sub_model_cons = subModel.user_defined_constraints();
  userDefinedConstraints = (reshape_vars) ?
    Constraints(recast_svd) : sub_model_cons.copy();
  // the recast_secondary_offset cannot in general be inferred from the
  // contributing fns in secondaryRespMapIndices (recast constraints may be
  // defined, e.g., with no contributing fns), and must therefore be passed.
  const size_t& num_recast_nln_ineq = recast_secondary_offset;
  size_t num_recast_nln_eq = num_recast_secondary_fns - num_recast_nln_ineq;
  if ( num_recast_nln_ineq != sub_model_cons.num_nonlinear_ineq_constraints()
    || num_recast_nln_eq   != sub_model_cons.num_nonlinear_eq_constraints() )
    userDefinedConstraints.reshape(num_recast_nln_ineq, num_recast_nln_eq,
      sub_model_cons.num_linear_ineq_constraints(),
      sub_model_cons.num_linear_eq_constraints());
}


/** This function is used for late initialization of the recasting
    functions.  It is used in concert with the alternate constructor. */
void RecastModel::
initialize(const Sizet2DArray& vars_map_indices,
	   bool nonlinear_vars_mapping,
	   void (*variables_map)      (const Variables& recast_vars,
				       Variables& sub_model_vars),
	   void (*set_map)            (const Variables& recast_vars,
				       const ActiveSet& recast_set,
				       ActiveSet& sub_model_set),
	   const Sizet2DArray& primary_resp_map_indices,
	   const Sizet2DArray& secondary_resp_map_indices,
	   const BoolDequeArray& nonlinear_resp_mapping,
	   void (*primary_resp_map)   (const Variables& sub_model_vars,
				       const Variables& recast_vars,
				       const Response& sub_model_response,
				       Response& recast_response),
	   void (*secondary_resp_map) (const Variables& sub_model_vars,
				       const Variables& recast_vars,
				       const Response& sub_model_response,
				       Response& recast_response))
{
  varsMapIndices          = vars_map_indices;
  nonlinearVarsMapping    = nonlinear_vars_mapping;
  variablesMapping        = variables_map;
  setMapping              = set_map;
  primaryRespMapIndices   = primary_resp_map_indices;
  secondaryRespMapIndices = secondary_resp_map_indices;
  nonlinearRespMapping    = nonlinear_resp_mapping;
  primaryRespMapping      = primary_resp_map;
  secondaryRespMapping    = secondary_resp_map;

  respMapping = (primaryRespMapping || secondaryRespMapping);

  if (nonlinearRespMapping.size() != primaryRespMapIndices.size() +
      secondaryRespMapIndices.size()) {
    Cerr << "Error: size mismatch in response mapping configuration." << endl;
    abort_handler(-1);
  }
}


void RecastModel::inverse_mappings(
    void (*inv_vars_map)     (const Variables& recast_vars,
			      Variables& sub_model_vars),
    void (*inv_set_map)      (const Variables& recast_vars,
			      const ActiveSet& recast_set,
			      ActiveSet& sub_model_set),
    void (*inv_pri_resp_map) (const Variables& sub_model_vars,
			      const Variables& recast_vars,
			      const Response& sub_model_resp,
			      Response& recast_resp),
    void (*inv_sec_resp_map) (const Variables& sub_model_vars,
			      const Variables& recast_vars,
			      const Response& sub_model_resp,
			      Response& recast_resp))
{
  //inverseMapFlag  = true;
  invVarsMapping    = inv_vars_map;     invSetMapping     = inv_set_map;
  invPriRespMapping = inv_pri_resp_map; invSecRespMapping = inv_sec_resp_map;
}


/** The RecastModel is evaluated by an Iterator for a recast problem
    formulation.  Therefore, the currentVariables, incoming active set,
    and output currentResponse all correspond to the recast inputs/outputs. */
void RecastModel::derived_compute_response(const ActiveSet& set)
{
  // transform from recast (Iterator) to sub-model (user) variables
  transform_variables(currentVariables, subModel.current_variables());

  // the incoming set is for the recast problem, which must be converted
  // back to the underlying response set for evaluation by the subModel.
  ActiveSet sub_model_set;
  transform_set(currentVariables, set, sub_model_set);

  // evaluate the subModel in the original fn set definition.  Doing this here 
  // eliminates the need for eval tracking logic within the separate eval fns.
  subModel.compute_response(sub_model_set);

  // recast the subModel response ("user space") into the currentResponse
  // ("iterator space")
  currentResponse.active_set(set);
  if (respMapping)
    transform_response(currentVariables, subModel.current_variables(),
		       subModel.current_response(), currentResponse);
  else
    currentResponse.update(subModel.current_response());

#ifdef DEBUG
  Cout << "Recast variables:\n"   << currentVariables
       << "subModel variables:\n" << subModel.current_variables()
       << "subModel response:\n"  << subModel.current_response()
       << "Recast response:\n"    << currentResponse;
#endif
}


void RecastModel::derived_asynch_compute_response(const ActiveSet& set)
{
  // transform from recast (Iterator) to sub-model (user) variables
  transform_variables(currentVariables, subModel.current_variables());

  // the incoming set is for the recast problem, which must be converted
  // back to the underlying response set for evaluation by the subModel.
  ActiveSet sub_model_set;
  transform_set(currentVariables, set, sub_model_set);

  // evaluate the subModel in the original fn set definition.  Doing this here 
  // eliminates the need for eval tracking logic within the separate eval fns.
  subModel.asynch_compute_response(sub_model_set);

  // bookkeep variables for use in primaryRespMapping/secondaryRespMapping
  if (respMapping) {
    int eval_id = subModel.evaluation_id();
    recastSetMap[eval_id]  = set;
    recastVarsMap[eval_id] = currentVariables.copy();
    if (variablesMapping)
      subModelVarsMap[eval_id] = subModel.current_variables().copy();
  }
}


const IntResponseMap& RecastModel::derived_synchronize()
{
  const IntResponseMap& orig_resp_map = subModel.synchronize();
  if (respMapping) {
    recastResponseMap.clear();
    IntASMIter     rsm_it = recastSetMap.begin();
    IntVarsMIter   rvm_it = recastVarsMap.begin();
    IntVarsMIter  smvm_it = (variablesMapping) ? subModelVarsMap.begin()
                                               : recastVarsMap.begin();
    IntRespMCIter map_cit = orig_resp_map.begin();
    for (; map_cit != orig_resp_map.end();
	 ++map_cit, ++rsm_it, ++rvm_it, ++smvm_it) {
      currentResponse.active_set(rsm_it->second);
      transform_response(rvm_it->second, smvm_it->second, map_cit->second,
			 currentResponse);
      recastResponseMap[map_cit->first] = currentResponse.copy();
    }
    recastSetMap.clear();
    recastVarsMap.clear();
    if (variablesMapping)
      subModelVarsMap.clear();
    return recastResponseMap;
  }
  else
    return orig_resp_map;
}


const IntResponseMap& RecastModel::derived_synchronize_nowait()
{
  const IntResponseMap& orig_resp_map = subModel.synchronize_nowait();
  if (respMapping) {
    recastResponseMap.clear();
    for (IntRespMCIter map_cit = orig_resp_map.begin();
	 map_cit != orig_resp_map.end(); ++map_cit) {
      int eval_id = map_cit->first;
      // IntResponseMap from subModel.synchronize_nowait() must be
      // consistent with subModel.evaluation_id() used above
      const Variables& sub_model_vars = (variablesMapping) ?
	subModelVarsMap[eval_id] : recastVarsMap[eval_id];
      currentResponse.active_set(recastSetMap[eval_id]);
      transform_response(recastVarsMap[eval_id], sub_model_vars,
			 map_cit->second, currentResponse);
      recastResponseMap[eval_id] = currentResponse.copy();
      recastSetMap.erase(eval_id);
      recastVarsMap.erase(eval_id);
      if (variablesMapping)
	subModelVarsMap.erase(eval_id);
    }
    return recastResponseMap;
  }
  else
    return orig_resp_map;
}


void RecastModel::
transform_variables(const Variables& recast_vars, Variables& sub_model_vars)
{
  // typical flow: mapping from recast variables ("iterator space")
  // into the sub-model variables ("user space")
  if (variablesMapping) variablesMapping(recast_vars, sub_model_vars);
  else                  sub_model_vars.active_variables(recast_vars);
}


void RecastModel::
inverse_transform_variables(const Variables& sub_model_vars,
			    Variables& recast_vars)
{
  // atypical flow: mapping from sub-model variables ("user space")
  // into the recast variables ("iterator space")
  if (invVarsMapping) invVarsMapping(sub_model_vars, recast_vars);
  else                recast_vars.active_variables(sub_model_vars);
}


void RecastModel::
transform_set(const Variables& recast_vars, const ActiveSet& recast_set,
	      ActiveSet& sub_model_set)
{
  // typical flow: mapping from recast set ("iterator space") into the
  // sub-model set ("user space")

  size_t i, j, num_recast_primary_fns = primaryRespMapIndices.size(),
    num_recast_secondary_fns = secondaryRespMapIndices.size(),
    num_recast_fns = num_recast_primary_fns + num_recast_secondary_fns;
  const ShortArray& recast_asv = recast_set.request_vector();
  if (recast_asv.size() != num_recast_fns) {
    Cerr << "Error: inconsistent asv sizing in RecastModel::transform_set().\n"
	 << "       recast asv size = " << recast_asv.size() << '\n'
	 << "       recast functions = " << num_recast_fns << endl;
    abort_handler(-1);
  }

  // Define default request vector and derivative vector mappings:
  // For the ASV, project each recast_asv request onto the contributing
  // set of functions within the sub_model_asv.  In the case of nonlinear
  // input/output mappings, the recast_asv request is augmented with
  // additional data requirements derived from chain rule differentiation.
  // The default sub-model DVV is just a copy of the recast DVV.
  ShortArray sub_model_asv(subModel.num_functions(), 0);
  for (i=0; i<num_recast_fns; i++) {
    short asv_val = recast_asv[i];
    // For nonlinear variable mappings, gradient required to transform Hessian.
    // A single nonlinear variable mapping affects all function derivatives.
    if (nonlinearVarsMapping && (asv_val & 4))
      asv_val |= 2;
    // assign the asv_val to each contributing sub-model function
    const SizetArray& recast_fn_contributors = (i<num_recast_primary_fns) ?
      primaryRespMapIndices[i] :
      secondaryRespMapIndices[i-num_recast_primary_fns];
    size_t num_contributors = recast_fn_contributors.size();
    for (j=0; j<num_contributors; j++) {
      short sub_model_asv_val = asv_val;
      // Bit deletions: for NLS recasting for full Newton without LeastSq term
      // Hessians, could remove 4 bit based on {gradient,hessian}Type, but this
      // is better accomplished from an Iterator's configuration using the
      // setMapping plug-in below (e.g., see Optimizer::gnewton_set_recast()). 

      // Bit additions: for nonlinear resp mappings, derivatives require all
      // lower order data. The nonlinearity of each fn contribution is employed.
      if (nonlinearRespMapping[i][j]) {
	if (asv_val & 4)
	  sub_model_asv_val |= 3;
	else if (asv_val & 2)
	  sub_model_asv_val |= 1;
      }
      sub_model_asv[recast_fn_contributors[j]] |= sub_model_asv_val;
    }
  }
  sub_model_set.request_vector(sub_model_asv);
  sub_model_set.derivative_vector(recast_set.derivative_vector()); // copy

  // a setMapping (provided in the RecastModel ctor or initialize()) augments
  // the standard mappings.  Current examples include NonD::set_u_to_x_mapping,
  // NonDReliability::PMA2_set_mapping, and Optimizer::gauss_newton_set_recast.
  // This follows the standard mappings so that provided mappings don't get
  // overwritten by the standard logic.  However, this means that any provided
  // additions will not be automatically augmented by nonlinear mapping logic
  // above.  This should not be a significant problem, since the provided
  // additions have case-specific context whereas the logic above is generic.
  // It would be preferable if provided mappings focused on updating the
  // sub_model_set rather than generating it from recast_set.
  if (setMapping)
    setMapping(recast_vars, recast_set, sub_model_set);
}


void RecastModel::
inverse_transform_set(const Variables& sub_model_vars,
		      const ActiveSet& sub_model_set, ActiveSet& recast_set)
{
  // atypical flow: mapping from sub-model set ("user space") into the
  // recast set ("iterator space")

  /* TO DO: modify mapping below from forward to inverse

  size_t i, j, num_recast_primary_fns = primaryRespMapIndices.size(),
    num_recast_secondary_fns = secondaryRespMapIndices.size(),
    num_recast_fns = num_recast_primary_fns + num_recast_secondary_fns;
  const ShortArray& recast_asv = recast_set.request_vector();
  if (recast_asv.size() != num_recast_fns) {
    Cerr << "Error: inconsistent asv sizing in RecastModel::"
         << "inverse_transform_set()." << std::endl;
    abort_handler(-1);
  }

  // Define default request vector and derivative vector mappings:
  // For the ASV, project each recast_asv request onto the contributing
  // set of functions within the sub_model_asv.  In the case of nonlinear
  // input/output mappings, the recast_asv request is augmented with
  // additional data requirements derived from chain rule differentiation.
  // The default sub-model DVV is just a copy of the recast DVV.
  ShortArray sub_model_asv(subModel.num_functions(), 0);
  for (i=0; i<num_recast_fns; i++) {
    short asv_val = recast_asv[i];
    // For nonlinear variable mappings, gradient required to transform Hessian.
    // A single nonlinear variable mapping affects all function derivatives.
    if (nonlinearVarsMapping && (asv_val & 4))
      asv_val |= 2;
    // assign the asv_val to each contributing sub-model function
    const SizetArray& recast_fn_contributors = (i<num_recast_primary_fns) ?
      primaryRespMapIndices[i] :
      secondaryRespMapIndices[i-num_recast_primary_fns];
    size_t num_contributors = recast_fn_contributors.size();
    for (j=0; j<num_contributors; j++) {
      short sub_model_asv_val = asv_val;
      // Bit deletions: for NLS recasting for full Newton without LeastSq term
      // Hessians, could remove 4 bit based on {gradient,hessian}Type, but this
      // is better accomplished from an Iterator's configuration using the
      // setMapping plug-in below (e.g., see Optimizer::gnewton_set_recast()). 

      // Bit additions: for nonlinear resp mappings, derivatives require all
      // lower order data. The nonlinearity of each fn contribution is employed.
      if (nonlinearRespMapping[i][j]) {
	if (asv_val & 4)
	  sub_model_asv_val |= 3;
	else if (asv_val & 2)
	  sub_model_asv_val |= 1;
      }
      sub_model_asv[recast_fn_contributors[j]] |= sub_model_asv_val;
    }
  }
  sub_model_set.request_vector(sub_model_asv);
  sub_model_set.derivative_vector(recast_set.derivative_vector()); // copy
  */

  // an invSetMapping (provided in inverse_mappings()) augments the standard
  // mappings above, such that the provided mappings don't get overwritten by
  // the standard logic.
  if (invSetMapping)
    invSetMapping(sub_model_vars, sub_model_set, recast_set);
}


void RecastModel::
transform_response(const Variables& recast_vars,
		   const Variables& sub_model_vars,
		   const Response& sub_model_resp, Response& recast_resp)
{
  // typical flow: mapping from sub-model response ("user space") into
  // the recast response ("iterator space")

  size_t num_recast_1_fns = primaryRespMapIndices.size();

  if (primaryRespMapping)
    primaryRespMapping(sub_model_vars, recast_vars,
		       sub_model_resp, recast_resp);
  else // number of recast primary = number of sub-model primary
    recast_resp.update_partial(0, num_recast_1_fns, sub_model_resp, 0);

  if (secondaryRespMapping)
    secondaryRespMapping(sub_model_vars, recast_vars,
			 sub_model_resp, recast_resp);
  else {
    // number of recast secondary = number of sub-model secondary,
    // but primary offsets may differ
    size_t num_recast_2_fns = secondaryRespMapIndices.size(),
           num_sm_1_fns     = sub_model_resp.num_functions() - num_recast_2_fns;
    recast_resp.update_partial(num_recast_1_fns, num_recast_2_fns,
			       sub_model_resp, num_sm_1_fns);
  }
}


void RecastModel::
inverse_transform_response(const Variables& sub_model_vars,
			   const Variables& recast_vars,
			   const Response& recast_resp,
			   Response& sub_model_resp)
{
  // atypical flow: mapping from the recast response ("iterator space")
  // into the sub-model response ("user space")

  size_t num_recast_1_fns = primaryRespMapIndices.size();

  if (invPriRespMapping)
    invPriRespMapping(recast_vars, sub_model_vars, recast_resp, sub_model_resp);
  else // number of recast primary = number of sub-model primary
    sub_model_resp.update_partial(0, num_recast_1_fns, recast_resp, 0);

  if (invSecRespMapping)
    invSecRespMapping(recast_vars, sub_model_vars, recast_resp, sub_model_resp);
  else {
    // number of recast secondary = number of sub-model secondary,
    // but primary offsets may differ
    size_t num_recast_2_fns = secondaryRespMapIndices.size(),
           num_sm_1_fns     = sub_model_resp.num_functions() - num_recast_2_fns;
    sub_model_resp.update_partial(num_sm_1_fns, num_recast_2_fns,
				  recast_resp, num_recast_1_fns);
  }
}


void RecastModel::initialize_data_from_submodel()
{
  componentParallelMode = SUB_MODEL;
  outputLevel           = subModel.output_level();

  gradientType          = subModel.gradient_type();
  methodSource          = subModel.method_source();
  ignoreBounds          = subModel.ignore_bounds();
  centralHess	          = subModel.central_hess();
  intervalType          = subModel.interval_type();
  fdGradStepSize        = subModel.fd_gradient_step_size();
  fdGradStepType        = subModel.fd_gradient_step_type();
  gradIdAnalytic        = subModel.gradient_id_analytic();
  gradIdNumerical       = subModel.gradient_id_numerical();

  hessianType           = subModel.hessian_type();
  quasiHessType         = subModel.quasi_hessian_type();
  fdHessByFnStepSize    = subModel.fd_hessian_by_fn_step_size();
  fdHessByGradStepSize  = subModel.fd_hessian_by_grad_step_size();
  fdHessStepType        = subModel.fd_hessian_step_type();
  hessIdAnalytic        = subModel.hessian_id_analytic();
  hessIdNumerical       = subModel.hessian_id_numerical();
  hessIdQuasi           = subModel.hessian_id_quasi();

  scalingOpts           = subModel.scaling_options();
}


/** Update inactive values and labels in currentVariables and inactive
    bound constraints in userDefinedConstraints from variables and
    constraints data within subModel. */
void RecastModel::update_from_sub_model()
{
  currentVariables.inactive_continuous_variables(
    subModel.inactive_continuous_variables());
  currentVariables.inactive_discrete_int_variables(
    subModel.inactive_discrete_int_variables());
  currentVariables.inactive_discrete_real_variables(
    subModel.inactive_discrete_real_variables());

  userDefinedConstraints.inactive_continuous_lower_bounds(
    subModel.inactive_continuous_lower_bounds());
  userDefinedConstraints.inactive_continuous_upper_bounds(
    subModel.inactive_continuous_upper_bounds());
  userDefinedConstraints.inactive_discrete_int_lower_bounds(
    subModel.inactive_discrete_int_lower_bounds());
  userDefinedConstraints.inactive_discrete_int_upper_bounds(
    subModel.inactive_discrete_int_upper_bounds());
  userDefinedConstraints.inactive_discrete_real_lower_bounds(
    subModel.inactive_discrete_real_lower_bounds());
  userDefinedConstraints.inactive_discrete_real_upper_bounds(
    subModel.inactive_discrete_real_upper_bounds());

  currentVariables.inactive_continuous_variable_labels(
    subModel.inactive_continuous_variable_labels());
  currentVariables.inactive_discrete_int_variable_labels(
    subModel.inactive_discrete_int_variable_labels());
  currentVariables.inactive_discrete_real_variable_labels(
    subModel.inactive_discrete_real_variable_labels());

  if (invVarsMapping) {
    invVarsMapping(subModel.current_variables(), currentVariables);
    // BMA TODO: there may be cases where we also want to update the
    // constraints and values, but there's currently no mechanism to
    // do so.  The client of a RecastModel must manage this.
  } else if (variablesMapping) {
    // no reasonable default

    // can't just apply variables mapping to values/bounds, since need inverse
    // of variablesMapping to go from subModel vars to currentVariables

    // any label, uncertain variable distributions, and linear
    // constraint mappings must be performed explicitly

    // for partial mapping of distribution parameters that are unmodified by a
    // variable transformation, see NonDExpansion::initialize_expansion()
  }
  else {
    // variable values
    currentVariables.continuous_variables(subModel.continuous_variables());
    currentVariables.discrete_int_variables(subModel.discrete_int_variables());
    currentVariables.discrete_real_variables(
      subModel.discrete_real_variables());
    // variable bounds
    userDefinedConstraints.continuous_lower_bounds(
      subModel.continuous_lower_bounds());
    userDefinedConstraints.continuous_upper_bounds(
      subModel.continuous_upper_bounds());
    userDefinedConstraints.discrete_int_lower_bounds(
      subModel.discrete_int_lower_bounds());
    userDefinedConstraints.discrete_int_upper_bounds(
      subModel.discrete_int_upper_bounds());
    userDefinedConstraints.discrete_real_lower_bounds(
      subModel.discrete_real_lower_bounds());
    userDefinedConstraints.discrete_real_upper_bounds(
      subModel.discrete_real_upper_bounds());
    // variable labels
    currentVariables.continuous_variable_labels(
      subModel.continuous_variable_labels());
    currentVariables.discrete_int_variable_labels(
      subModel.discrete_int_variable_labels());
    currentVariables.discrete_real_variable_labels(
      subModel.discrete_real_variable_labels());

    if (!subModel.discrete_design_set_int_values().empty())
      discreteDesignSetIntValues = subModel.discrete_design_set_int_values();
    if (!subModel.discrete_design_set_real_values().empty())
      discreteDesignSetRealValues = subModel.discrete_design_set_real_values();

    // uncertain variable distribution data
    aleatDistParams.update(subModel.aleatory_distribution_parameters());
    epistDistParams.update(subModel.epistemic_distribution_parameters());

    if (!subModel.discrete_state_set_int_values().empty())
      discreteStateSetIntValues = subModel.discrete_state_set_int_values();
    if (!subModel.discrete_state_set_real_values().empty())
      discreteStateSetRealValues = subModel.discrete_state_set_real_values();

    // linear constraints
    if (subModel.num_linear_ineq_constraints()) {
      userDefinedConstraints.linear_ineq_constraint_coeffs(
        subModel.linear_ineq_constraint_coeffs());
      userDefinedConstraints.linear_ineq_constraint_lower_bounds(
        subModel.linear_ineq_constraint_lower_bounds());
      userDefinedConstraints.linear_ineq_constraint_upper_bounds(
        subModel.linear_ineq_constraint_upper_bounds());
    }
    if (subModel.num_linear_eq_constraints()) {
      userDefinedConstraints.linear_eq_constraint_coeffs(
        subModel.linear_eq_constraint_coeffs());
      userDefinedConstraints.linear_eq_constraint_targets(
        subModel.linear_eq_constraint_targets());
    }
  }

  if (primaryRespMapping) {
    // response mappings are in opposite direction from variables
    // mappings, so primaryRespMapping could potentially be used to
    // update currentResponse from subModel primary fns
  }
  else {
    // primary response function weights
    primaryRespFnWts = subModel.primary_response_fn_weights();
    // primary response function sense (min or max)
    primaryRespFnSense = subModel.primary_response_fn_sense();

    // primary response function labels
    const StringArray& sm_resp_labels = subModel.response_labels();
    size_t i, num_primary = numFns 
      - userDefinedConstraints.num_nonlinear_eq_constraints()
      - userDefinedConstraints.num_nonlinear_ineq_constraints();
    for (i=0; i<num_primary; i++)
      currentResponse.shared_data().function_label(sm_resp_labels[i], i);
  }

  if (secondaryRespMapping) {
    // response mappings are in opposite direction from variables
    // mappings, so secondaryRespMapping could potentially be used to
    // update currentResponse from subModel secondary fns
  }
  else {
    // secondary response function labels
    const StringArray& sm_resp_labels = subModel.response_labels();
    size_t i,
      num_nln_con = userDefinedConstraints.num_nonlinear_eq_constraints() +
        userDefinedConstraints.num_nonlinear_ineq_constraints(),
      num_primary    = numFns - num_nln_con,
      num_sm_primary = subModel.num_functions() - num_nln_con;
    for (i=0; i<num_nln_con; i++)
      currentResponse.shared_data().function_label(
	sm_resp_labels[num_sm_primary+i], num_primary+i);

    // nonlinear constraint bounds/targets
    if (subModel.num_nonlinear_ineq_constraints()) {
      userDefinedConstraints.nonlinear_ineq_constraint_lower_bounds(
        subModel.nonlinear_ineq_constraint_lower_bounds());
      userDefinedConstraints.nonlinear_ineq_constraint_upper_bounds(
        subModel.nonlinear_ineq_constraint_upper_bounds());
    }
    if (subModel.num_nonlinear_eq_constraints())
      userDefinedConstraints.nonlinear_eq_constraint_targets(
        subModel.nonlinear_eq_constraint_targets());
  }
}


bool RecastModel::
db_lookup(const Variables& search_vars, const ActiveSet& search_set,
	  Response& found_resp)
{
  // transform from recast (Iterator) to sub-model (user) variables;
  // making copy to avoid modifying submodel state during the lookup
  Variables sub_model_vars(subModel.current_variables().copy());
  transform_variables(search_vars, sub_model_vars);

  // the incoming set is for the recast problem, which must be converted
  // back to the underlying response set for evaluation by the subModel.
  ActiveSet sub_model_set;
  transform_set(search_vars, search_set, sub_model_set);

  // invoke default implementation for the lookup; making copy to
  // avoid modifying submodel state during the lookup
  Response sub_model_resp(subModel.current_response().copy());
  bool eval_found = Model::db_lookup(sub_model_vars, search_set, sub_model_resp);
  if (!eval_found)
    return false;

  // recast the subModel response ("user space") into the "iterator space"
  found_resp.active_set(search_set);
  if (respMapping)
    transform_response(search_vars, sub_model_vars, sub_model_resp, found_resp);
  else
    found_resp.update(sub_model_resp);

  return eval_found;
}


} // namespace Dakota
