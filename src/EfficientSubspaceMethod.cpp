/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:	 EfficientSubspaceMethod
//- Description: Implementation code for EfficientSubspaceMethod class
//- Owner:       Brian Adams
//- Version:

//- Last edited by: Brian Adams on 01/16/2013

#include "ProblemDescDB.hpp"
#include "EfficientSubspaceMethod.hpp"
#include "NonDLHSSampling.hpp"
#include "RecastModel.hpp"
#include "dakota_linear_algebra.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "BootstrapSampler.hpp"
#include <vector>
#include <cmath>

// TODO: Ultimately (like NonDAdaptiveSampling) want this method to
// result in a approximation (reduced model) that can be used in other
// contexts.

// TODO: Review use of linear algebra and dimensions throughout

// TODO: Limit reduced rank to size of original space (done?)

// TODO: consider using a single tolerance and using the derivative to
// estimate the relationship between the ytol and svdtol

// TODO: probably for now should set the nullspaceTol based on
// convergenceTol and make the user singular value more generous

// TODO: error if convergenceTol < macheps
// TODO: recon tol needs to be bounded below too

// TODO: Add a convergence manager class

// TODO: subtract UQ samples from the original function budget

namespace Dakota {

EfficientSubspaceMethod* EfficientSubspaceMethod::esmInstance(NULL);

EfficientSubspaceMethod::
EfficientSubspaceMethod(ProblemDescDB& problem_db, Model& model):
  NonD(problem_db, model),
  seedSpec(probDescDB.get_int("method.random_seed")),
  initialSamples(probDescDB.get_int("method.samples")),    // default 0
  batchSize(probDescDB.get_int("method.nond.batch_size")), // default 0
  subspaceSamples(probDescDB.get_int("method.nond.emulator_samples")), // def 0
  currIter(0), totalSamples(0), totalEvals(0), userSVTol(convergenceTol),
  nullspaceTol(convergenceTol/1.0e3), svRatio(0.0), reducedRank(0)
{
  // the Iterator initializes:
  //   maxIterations    (default -1)
  //   convergenceTol   (default 1.0e-4); tolerance before checking recon error
  //   maxFunctionEvals (default 1000)
  validate_inputs();

  // initialize the fullspace Monte Carlo derivative sampler; this
  // will configure it to perform initialSamples
  init_fullspace_sampler();

  // would be useful here for setting concurrency, but don't know the
  // sizes yet
  //init_reduced_sampler(subspaceSamples);

  // The same model is used in multiple contexts, with varying degrees
  // of concurrency: initialSamples, batchSize, verif_samples, and
  // finally subspaceSamples.  For now, we configure max eval
  // concurrency for the initialSamples and allow batches in
  // subsequent phases; could instead take max of all these if
  // desired.
  maxEvalConcurrency *= initialSamples;
}


EfficientSubspaceMethod::~EfficientSubspaceMethod()
{
  /* empty destructor */
}


void EfficientSubspaceMethod::quantify_uncertainty()
{
  esmInstance = this;

  Cout << "ESM: Performing sampling to build reduced space" << std::endl;

  // whether user singular value tolerance met
  bool user_svtol_met = false;
  // whether macheps singular value tolerance met
  bool mach_svtol_met = false;
  // whether reconstruction error tolerance met
  bool recon_tol_met = false;

  // iterate until numerically singular or small reconstruction error
  while (!mach_svtol_met && !recon_tol_met &&
	 currIter < maxIterations && totalEvals < maxFunctionEvals &&
	 reducedRank < numContinuousVars) {

    // Run the inner iteration until user tol met (if possible), then check
    // reconstruction error, then continue outer if not tight enough
    if (currIter == 0) {

      // initially do this loop until user tolerance met
      // once met (userSVTol will always be met),
      while (!user_svtol_met && !mach_svtol_met &&
	     currIter < maxIterations && totalEvals < maxFunctionEvals
	     && reducedRank < numContinuousVars) {

	++currIter;  // any addition of batch of points counts as an iteration
	Cout << "\nESM: Iteration " << currIter << "." << std::endl;
	expand_basis(mach_svtol_met, user_svtol_met);
	print_svd_stats();

      } // until SVD converged

      if (user_svtol_met || mach_svtol_met) {
	Cout << "\nESM: SVD converged to tolerance.\n     Proceeding to "
	     << "reconstruction with reduced rank = " << reducedRank << "."
	     << std::endl;
      }
      else {
	Cout << "\nESM: SVD not converged within budget.";
	if (currIter >= maxIterations)
	  Cout << "\n     Maximum iterations reached.";
	if (totalEvals >= maxFunctionEvals)
	  Cout << " \n    Maximum function evaluations reached.";
	Cout << "\n    Proceeding to reconstruction with reducedRank = "
	     << reducedRank << std::endl;
      }

    }
    else {
	++currIter;  // any addition of batch of points counts as an iteration
	Cout << "\nESM: Iteration " << currIter << "." << std::endl;
	expand_basis(mach_svtol_met, user_svtol_met);
	print_svd_stats();
    }

    // update the reducedBasis
    // the reduced basis is dimension N x r and stored in the first r
    // cols of derivativeMatrix; extract it instead of using BLAS directly
    // TODO: could probably do a View and avoid the Copy
    RealMatrix reduced_basis_U(Teuchos::Copy, derivativeMatrix,
			       numContinuousVars, reducedRank);
    reducedBasis = reduced_basis_U;
    if (outputLevel >= DEBUG_OUTPUT) {
      Cout << "\nESM: Reduced basis is";
      reducedBasis.print(Cout);
    }

    // evaluate the fidelity of the reconstruction via orthogonal subspace
    // evaluations (constrained to stay in bounds)
    assess_reconstruction(recon_tol_met);

    if (!recon_tol_met)
      Cout << "\nESM: Reconstruction tolerance not met." << std::endl;

  } // until reconstruction converged

  Cout << "\n --- ESM Final Build Statistics ---"
       << "\n  total iterations: " << currIter
       << "\n  small/large singular value: " << svRatio
       << "\n  build samples: " << totalSamples
       << "\n  reduced rank: " << reducedRank
       << std::endl;

  Cout << "\n --- ESM Build Convergence Criteria ---"
       << "\n  user tolerance on SVD met?: " << user_svtol_met
       << "\n  macheps tolerance on SVD met?: " << mach_svtol_met
       << "\n  reconstruction tolerance met?: " << recon_tol_met
       << "\n  max_iterations reached: " << (bool) (currIter >= maxIterations)
       << "\n  max_evals reached: " << (bool) (totalEvals >= maxFunctionEvals)
       << "\n  max rank reached: " << (bool) (reducedRank >= numContinuousVars)
       << std::endl;

  // perform the reduced space UQ
  Cout << "\nESM: Performing reduced-space UQ" << std::endl;
  reduced_space_uq();

}


void EfficientSubspaceMethod::validate_inputs()
{
  bool error_flag = false;

  // validate iteration controls

  // set default initialSamples, with lower bound of 2
  // TODO: allow other user control of initial sample rule?
  if (initialSamples <= 0) {
    initialSamples =
      (unsigned int) std::ceil( (double) numContinuousVars / 100.0 );
    initialSamples = std::max(2, initialSamples);
    Cout << "\nInfo: Efficient subspace method setting (initial) samples = "
	 << initialSamples << "." << std::endl;
  }
  else if (initialSamples < 2) {
    initialSamples = 2;
    Cout << "\nWarning: Efficient subspace method resetting samples to minimum "
	 << "allowed = " << initialSamples << "." << std::endl;
  }

  if (initialSamples > maxFunctionEvals) {
    error_flag = true;
    Cerr << "\nError: Efficient subspace method build samples exceeds function "
	 << "budget." << std::endl;
  }

  if (batchSize <= 0) {
    // default is to add one point at a time
    batchSize = 1;
  }
  else if (batchSize > initialSamples) {
    batchSize = initialSamples;
    Cout << "\nWarning: batch_size = " << batchSize << " exceeds (initial) "
	 << "samples = " << initialSamples << ";\n        resetting batch_size "
	 << "= " << batchSize << "." << std::endl;
  }

  // maxIterations controls the number of build iterations
  if (maxIterations < 0) {
    maxIterations = 1;
    Cout << "\nInfo: Efficient subspace method setting max_iterations = "
	 << maxIterations << "." << std::endl;
  }

  // emulator samples don't count toward the sample budget
  if (subspaceSamples <= 0) {
    subspaceSamples = 10*initialSamples;
    Cout << "\nInfo: Efficient subspace method setting emulator_samples = "
	 << subspaceSamples << "\n      (10*samples specified)." << std::endl;
  }


  // validate variables specification
  if (numContinuousVars != numNormalVars
      || numDiscreteIntVars > 0 || numDiscreteRealVars > 0) {
    error_flag = true;
    Cerr << "\nError: Efficient subspace method only supports normal uncertain "
         << "variables;\n       remove other variable specifications."
	 << std::endl;
  }

  // validate response data
  if (iteratedModel.gradient_type() == "none") {
    error_flag = true;
    Cerr << "\nError: Efficient subspace method requires gradients.\n"
         << "       Please select numerical, analytic (recommended), or mixed "
	 << "gradients." << std::endl;
  }

  if (error_flag)
    abort_handler(-1);
}


void EfficientSubspaceMethod::init_fullspace_sampler()
{
  // ----
  // Instantiate the DACE iterator, using default RNG and samples
  // ----
  int mc_seed = (seedSpec) ? seedSpec : generate_system_seed();
  if (seedSpec)
    Cout << "ESM: build seed (user-specified) = " << mc_seed << std::endl;
  else
    Cout << "ESM: build seed (system-generated) = " << mc_seed << std::endl;

  // use Monte Carlo due to iterative growth process
  unsigned short sample_type = SUBMETHOD_RANDOM;
  std::string rng; // use default random number generator

  // configure this sampler initially to work with initialSamples
  Analyzer* ndlhss =
    new NonDLHSSampling(iteratedModel, sample_type, initialSamples, mc_seed,
			rng, ACTIVE_UNIFORM);
  // allow random number sequence to span multiple calls to run()
  ndlhss->vary_pattern(true);
  fullSpaceSampler.assign_rep(ndlhss, false);

  // TODO: review whether this is needed
  //fullSpaceSampler.sub_iterator_flag(true);
}


/**  This specialization is because the model is used in multiple
     contexts in this iterator, depending on build phase.  Note that
     this overrides the default behavior at Iterator which recurses
     into any submodels. */
void EfficientSubspaceMethod::derived_init_communicators(ParLevLIter pl_iter)
{
  // In ESM, the same model will be used in multiple contexts:
  //  - fullSpaceSampler(iteratedModel) with initialSamples
  //  - fullSpaceSampler(iteratedModel) with batchSamples
  //  - direct compute_response() of verif_samples, one at a time
  //  - reduced_space_sampler(RecastModel(iteratedModel)) with subspaceSamples

  // Instead of using the helper iterators convenience function, we
  // directly set up communicators for each of the contexts that will
  // be encountered at run-time
  iteratedModel.init_communicators(pl_iter, maxEvalConcurrency);//initialSamples
  iteratedModel.init_communicators(pl_iter, batchSize);
  // defer this one to do it on the RecastModel at runtime
  //  iteratedModel.init_communicators(pl_iter, subspaceSamples);
}


void EfficientSubspaceMethod::derived_set_communicators(ParLevLIter pl_iter)
{
  miPLIndex = methodPCIter->mi_parallel_level_index(pl_iter);
  iteratedModel.set_communicators(pl_iter, maxEvalConcurrency);//initialSamples
}


void EfficientSubspaceMethod::derived_free_communicators(ParLevLIter pl_iter)
{
  iteratedModel.free_communicators(pl_iter, batchSize);
  iteratedModel.free_communicators(pl_iter, maxEvalConcurrency);//initialSamples
  // defer this one to do it on the RecastModel at runtime
  //  iteratedModel.free_communicators(pl_iter, subspaceSamples);
}


void EfficientSubspaceMethod::
expand_basis(bool& mach_svtol_met, bool& user_svtol_met)
{
  // determine number of points to add
  unsigned int diff_samples = calculate_fullspace_samples();
  totalSamples += diff_samples;
  totalEvals += diff_samples;

  if (outputLevel >= NORMAL_OUTPUT)
    Cout << "ESM: Adding " << diff_samples << " full-space samples."
	 << std::endl;

  // evaluate samples with fullSpaceSampler
  generate_fullspace_samples(diff_samples);

  // add the generated points to the matrices
  append_sample_matrices(diff_samples);

  // factor the derivative matrix and estimate the information content
  compute_svd(mach_svtol_met, user_svtol_met);
}


unsigned int EfficientSubspaceMethod::calculate_fullspace_samples()
{
  // on first pass, use initial build samples
  if (currIter == 1)
    return initialSamples;

  unsigned int diff_samples;

  // for now only supporting growth by fixed batch size
  diff_samples = batchSize;

  // TODO: allow options for fullspace sample growth based on
  // fixed user-supplied factor, logarithmic growth

  // factor by which to increase number samples at each minor iteration
  // therefore also each major iteration (> 1.0) (not used for now)
  //double sample_growth_factor = 2.0;
  //double desired_total_samples = sample_growth_factor*totalSamples;
  //diff_samples =
  //  (unsigned int) std::ceil(desired_total_samples) - totalSamples;

  // limit the increment by maximum evaluations control
  diff_samples = std::min(diff_samples, maxFunctionEvals - totalEvals);

  return diff_samples;
}


void EfficientSubspaceMethod::
generate_fullspace_samples(unsigned int diff_samples)
{
  // Rank-revealing phase requires derivatives (for now)
  ActiveSet dace_set = fullSpaceSampler.active_set(); // copy
  // TODO: always gradients only; no functions
  //       analysis_driver needs to parse active_set
  unsigned short request_value = 3;
  dace_set.request_values(request_value);
  fullSpaceSampler.active_set(dace_set);

  // Generate the samples.  Have to adjust the base number of samples
  // with sampling_reference() since the number of samples may go down
  // from intialSamples to batchSamples
  fullSpaceSampler.sampling_reference(diff_samples);
  fullSpaceSampler.sampling_reset(diff_samples, true, false);

  // select the right comm for current model evaluation concurrency
  // (initial build or batch update)
  ParLevLIter pl_iter = methodPCIter->mi_parallel_level_iterator(miPLIndex);
  if (currIter == 1)
    iteratedModel.set_communicators(pl_iter, initialSamples);
  else
    iteratedModel.set_communicators(pl_iter, batchSize);
  // and generate the additional samples
  fullSpaceSampler.run();//(pl_iter);
}


void EfficientSubspaceMethod::
append_sample_matrices(unsigned int diff_samples)
{
  // extract into a matrix
  // all_samples vs. all_variables
  const RealMatrix& all_vars = fullSpaceSampler.all_samples();
  const IntResponseMap& all_responses = fullSpaceSampler.all_responses();

  // TODO: could easily filter NaN/Inf responses and omit
  if (outputLevel >= DEBUG_OUTPUT) {
    Cout << "\nESM: Iteration " << currIter << ". DACE iterator returned "
	     << all_responses.size() << " samples. (Expected diff_samples = "
	     << diff_samples << ".)" << std::endl;
  }

  int sample_insert_point = varsMatrix.numCols();
  derivativeMatrix.reshape(numContinuousVars,
			   totalSamples*numFunctions);
  varsMatrix.reshape(numContinuousVars, totalSamples);

  unsigned int diff_sample_ind = 0;
  IntRespMCIter resp_it  = all_responses.begin();
  IntRespMCIter resp_end = all_responses.end();

  for ( ; resp_it != resp_end ; ++resp_it, ++diff_sample_ind) {
    // the absolute sample number to insert into
    unsigned int sample_ind = sample_insert_point + diff_sample_ind;
    // matrix of num_variables x num_functions
    const RealMatrix& resp_matrix = resp_it->second.function_gradients();
    for (unsigned int fn_ind = 0; fn_ind < numFunctions; ++fn_ind) {
      unsigned int col_ind = sample_ind * numFunctions + fn_ind;
      for (unsigned int var_ind = 0; var_ind < numContinuousVars; ++var_ind) {
	derivativeMatrix(var_ind, col_ind) = resp_matrix(var_ind, fn_ind);
      }
    }
    for (unsigned int var_ind = 0; var_ind < numContinuousVars; ++var_ind) {
      varsMatrix(var_ind, sample_ind) = all_vars(var_ind, diff_sample_ind);
    }
  }

  if (outputLevel >= DEBUG_OUTPUT) {
    Cout << "\nESM: Iteration " << currIter
	 << ". Compiled derivative matrix is:\n";
    derivativeMatrix.print(Cout);
    Cout << std::endl;
  }
}


void EfficientSubspaceMethod::
compute_svd(bool& mach_svtol_met, bool& user_svtol_met)
{
  // Want eigenvalues of derivMatrix*derivMatrix^T, so perform SVD of
  // derivMatrix and square them

    RealVector singular_values;
    RealMatrix V_transpose;
    // BMA: This overwrites the derivativeMatrix with left singular
    // vectors, probably not what we want...
    svd(derivativeMatrix, singular_values, V_transpose);

////////////////////////////////////////////////////////////////////////////////
    // BEGIN: Computational kernel to compute Bing Li's criterion for the
    // eigenvalue gap

    // TODO: Analyze whether we need to worry about this
    if(singular_values.length() == 0)
    {
      Cerr << "No computed singular values available!" << std::endl;
      abort_handler(-1);
    }

    int num_vars = derivativeMatrix.numRows();
    // TODO: Analyze whether we need this check and can have differing numbers
    // of singular values returned
    if(num_vars != singular_values.length())
    {
      Cerr << "Number of computed singular_values does not match the dimension "
              "of the space of gradient samples! Logic not currently supported!"
           << std::endl;
      abort_handler(-1);
    }

    // Stores Bing Li's criterion
    std::vector<RealMatrix::scalarType> bing_li_criterion(num_vars, 0);

    // Compute part 1 of criterion: relative energy in next eigenvalue in the
    // spectrum

    RealMatrix::scalarType eigen_sum = singular_values[0] * singular_values[0];

    for(size_t i = 0; i < num_vars - 1; ++i)
    {
      RealMatrix::scalarType eigen_val = singular_values[i+1] *
                                         singular_values[i+1];
      bing_li_criterion[i] = eigen_val;
      eigen_sum += eigen_val;
    }

    for(size_t i = 0; i < num_vars; ++i)
      bing_li_criterion[i] /= eigen_sum;

    // Compute part 2 of criterion: bootstrapped determinant metric

    RealMatrix bootstrapped_sample(num_vars, derivativeMatrix.numCols());
    RealVector sample_sing_vals;
    RealMatrix sample_sing_vectors;

    Teuchos::LAPACK<RealMatrix::ordinalType, RealMatrix::scalarType> lapack;

    std::vector<RealMatrix::scalarType> bootstrapped_determinants(num_vars);

    BootstrapSampler<RealMatrix> bootstrap_sampler(derivativeMatrix,
                                                   numFunctions);

    // TODO: Get number of bootstrap samples from problem desc database
    size_t num_replicates = 100;
    for (size_t i = 0; i < num_replicates; ++i)
    {
      bootstrap_sampler(bootstrapped_sample);

      svd(bootstrapped_sample, sample_sing_vals, sample_sing_vectors);

      // Overwrite bootstrap replicate with singular matrix product
      bootstrapped_sample.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, 1.0,
                                   V_transpose, sample_sing_vectors, 0.0);

      for(size_t j = 0; j < num_vars; ++j)
      {
        size_t num_sing_vec = j + 1;

        RealMatrix submatrix(Teuchos::Copy, bootstrapped_sample, num_sing_vec,
                             num_sing_vec);

        // Get determinant from LU decomposition

        Teuchos::SerialDenseVector<RealMatrix::ordinalType,
                                   RealMatrix::ordinalType> pivot(num_sing_vec);
        RealMatrix::ordinalType info;

        lapack.GETRF(num_sing_vec, num_sing_vec, submatrix.values(),
                     num_sing_vec, pivot.values(), &info);

        RealMatrix::scalarType det = 1.0;
        for (size_t i = 0; i < j; ++i)
        {
          det *= submatrix(i,i);
        }

        bootstrapped_determinants[j] += std::abs(det);
      }
    }

    RealMatrix::scalarType det_sum = 0.0;
    for (size_t i = 0; i < num_vars; ++i)
    {
      bootstrapped_determinants[i] = 1 - bootstrapped_determinants[i] /
        static_cast<RealMatrix::scalarType>(num_replicates);
      det_sum += bootstrapped_determinants[i];
    }

    for (size_t i = 0; i < num_vars; ++i)
    {
      bing_li_criterion[i] += bootstrapped_determinants[i] / det_sum;
    }

    // Cutoff is minimum of the criterion
    reducedRank = 0;
    RealMatrix::scalarType criterion_min =
      bing_li_criterion[reducedRank];
    for (size_t i = 1; i < num_vars; ++i)
    {
      if(criterion_min > bing_li_criterion[i])
      {
        criterion_min = bing_li_criterion[i];
        reducedRank = i;
      }
    }

    // END: Computational kernel to compute Bing Li's criterion
////////////////////////////////////////////////////////////////////////////////

    // TODO: if a reducedRank met the tolerance, but we added more
    // samples to meet construction error, need to allow this bigger
    // than svdtol and not truncate
    double inf_norm = derivativeMatrix.normInf();

    // See Golub and VanLoan discussion of numerical rank; use
    // tolerance applied to the sup norm...

    // if first time, iterate until we meet user tolerance
    // otherwise continue until we're numerically rank deficient
    double user_svtol = inf_norm * userSVTol;
    double mach_svtol = inf_norm * std::numeric_limits<Real>::epsilon();

    int num_singular_values = singular_values.length();
    reducedRank = num_singular_values;
    double sv_small = 1., sv_large = DBL_MAX;
    for (unsigned int i=0; i<num_singular_values; ++i) {

      // we assume the user tolerance is looser than machine; could improve
      // this logic
      sv_small = singular_values[i];
      sv_large = singular_values[0];
      svRatio = singular_values[i] / singular_values[0];

      if (svRatio < mach_svtol) {
	mach_svtol_met = true;
	reducedRank = i;
	if (svRatio < user_svtol)
	  user_svtol_met = true;
	break;
      }
      else if (svRatio < user_svtol) {
	user_svtol_met = true;
	reducedRank = i;
	break;
      }
    }

    if (outputLevel >= DEBUG_OUTPUT) {
      Cout << "\nESM: Iteration " << currIter << ": singular values are [ ";
      for (unsigned int i=0; i<num_singular_values; ++i)
	Cout << singular_values[i] << " ";
      Cout << "]" << std::endl;
    }
}


void EfficientSubspaceMethod::print_svd_stats()
{
  if (outputLevel >= NORMAL_OUTPUT)
    Cout << "\n --- ESM Iteration " << currIter << " Statistics --- "
	 << "\n  small/large singular value: " << svRatio
	 << "\n  build samples: " << totalSamples
	 << "\n  estimated reduced rank: " << reducedRank
	 << std::endl;
}


/** This function is experimental and needs to be carefully reviewed
    and cleaned up */
void EfficientSubspaceMethod::
assess_reconstruction(bool& recon_tol_met)
{
  // calculate upper bound on number of verification samples
  // allow maxFunctionEvals to be exceeded
  // TODO: likely shouldn't do more than N-reducedRank samples?
  // TODO: Could generate additional random samples in the perp space
  unsigned int wilks_max_samples = 100; // maximum based on Wilk's formula
  unsigned int verif_samples = std::min(wilks_max_samples, totalSamples);
  if (verif_samples > maxFunctionEvals - totalEvals) {
    verif_samples = maxFunctionEvals - totalEvals;
    Cout << "\nESM: Warning: " << verif_samples << " verification samples "
	 << "desired, but budget only permits " << verif_samples << std::endl;
  }
  totalEvals += verif_samples;

  // Calculate perturbations of the samples from nominal and project
  // orthogonal to reducedBasis U

  // TODO: Relies on normal distribution for now
  const RealVector& nominal_vars =
    iteratedModel.aleatory_distribution_parameters().normal_means();

  // Find vectors orthogonal to each initial perturbation
  // for test problem, nominal is 0.5 for all x
  // TODO: use actual nominal values

  RealMatrix perturbations(numContinuousVars, verif_samples);
  for (int i=0; i<numContinuousVars; ++i)
    for (int j=0; j<verif_samples; ++j)
      perturbations(i,j) = varsMatrix(i,j) - nominal_vars(i);


  // want to compute perturbations delta_x_perp in the nullspace of
  // U^T, i.e.,
  //   delta_x_perp = (I - U*U^T)*delta_x = delta_x - U*U^T*delta_x

  // compute U^T * delta_x via
  // u_trans_delta_x = 1.0 * U^T * delta_x + 0.0 * u_trans_delta_x
  RealMatrix u_trans_delta_x(reducedRank, verif_samples); // init to zero

  Real alpha = 1.0;
  Real beta = 0.0;
  u_trans_delta_x.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, alpha,
			   reducedBasis, perturbations, beta);

  // compute -1.0 * ( U * u_trans_delta_x ) + 1.0 * delta_x
  alpha = -1.0;
  beta = 1.0;
  perturbations.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha,
			 reducedBasis, u_trans_delta_x, beta);

  // the matrix perturbations now contains the delta_x_perp
  // NEED to evaluate at nomimal + delta_x_perp

  // evaluate the model at the sample points in perturbations
  // for now in serial
  // TODO: initialize to nominal
  RealMatrix perp_points(numContinuousVars, verif_samples, false);
  for (int i=0; i<numContinuousVars; ++i)
    for (int j=0; j<verif_samples; ++j)
      perp_points(i,j) = nominal_vars(i);
  perp_points += perturbations;

  // TODO: restrict perp_points to the variable bounds

  // Reconstruction phase requires only functions
  // TODO: get off the model?
//   ActiveSet recon_set = dace_iterator.active_set(); // copy
//   unsigned short request_value = 1;
//   recon_set.request_values(request_value);
  activeSet.request_values(1);
  if (outputLevel >= NORMAL_OUTPUT)
    Cout << "\nESM: Evaluating at nominal, and at " << verif_samples
	 << " points orthogonal to the subspace" << std::endl;

  // evaluate model at nominal values
  iteratedModel.continuous_variables(nominal_vars);
  iteratedModel.compute_response(activeSet);
  const RealVector& ynominal =
    iteratedModel.current_response().function_values();

  // TODO: use array of RealVectors here

  // TODO: convenience function to evaluate a matrix of samples
  // TODO: asynch as well
  RealMatrix Kmat(verif_samples, numFunctions);
  for (int j=0; j<verif_samples; ++j) {
    iteratedModel.continuous_variables(getCol(Teuchos::View, perp_points, j));
    iteratedModel.compute_response(activeSet);
    // compute y(perp) - y(nominal)
    RealVector deviation(iteratedModel.current_response().function_values());
    deviation -= ynominal;
    Teuchos::setCol(deviation, j, Kmat);
  }

  recon_tol_met = true;
  RealVector recon_error(numFunctions);
  for (int i=0; i<numFunctions; ++i) {

    RealVector Kvec = Teuchos::getCol(Teuchos::View, Kmat, i);

    // two-norm of recon error for now
    double recon_error = std::sqrt(Kvec.dot(Kvec)/verif_samples);

    // standard deviation of Kvec
    double K_mu = std::accumulate(&Kvec[0], &Kvec[0]+verif_samples, 0.0) /
      (Real) verif_samples;
    double K_sigma = 0.0;
    // Could decimate Kvec with (K-mu)^2; opt for simplicity
    for (unsigned int j=0; j<verif_samples; ++j) {
      K_sigma += (Kvec(j) - K_mu)*(Kvec(j) - K_mu);
    }
    K_sigma = std::sqrt(K_sigma / (Real) (verif_samples-1));

    Cout << "\nESM: Reconstruction statistics over " << verif_samples
	 << " samples for function " << i << ":\n  K_sigma = " << K_sigma
	 << "\n  K_mu = " << K_mu
	 << "\n  L2 (recon) error = " << recon_error << std::endl;

    // if any function violates, continue process
    if (recon_error > nullspaceTol)
      recon_tol_met = false;
  }
}


/** This function is experimental and needs to be reviewed and cleaned
    up.  In particular the translation of the correlations from full
    to reduced space is likely wrong.  Transformation may be correct
    for covariance, but likely not correlations. */
void EfficientSubspaceMethod::reduced_space_uq()
{
  // Create a recast model such that LHS can sample over the reduced
  // variable space, but evaluate the original model back in the full
  // space
  size_t num_recast_vars = reducedRank;

  // doesn't appear to be used; for now, indicate that each submodel
  // variable depends on all of the recast (reduced) variables, since
  // in general, they will
  Sizet2DArray vars_map_indices(numContinuousVars,
				SizetArray(num_recast_vars));
  for (size_t i=0; i<numContinuousVars; ++i)
    for (size_t recasti=0; recasti<num_recast_vars; ++recasti)
      vars_map_indices[i][recasti] = recasti;

  bool nonlinear_vars_map = false; // vars map is linear

  // TODO: set map needed;

  // must be initialized! TODO: handle constraints?
  Sizet2DArray primary_resp_map_indices(numFunctions,
					SizetArray(numFunctions, 0));
  for (size_t i=0; i<numFunctions; ++i)
    primary_resp_map_indices[i][i] = i;

  // empty as no secondary mapping; could init to identity map
  Sizet2DArray secondary_resp_map_indices;

  size_t recast_secondary_offset = 0;

  BoolDequeArray nonlinear_resp_map(numFunctions,
				    BoolDeque(numFunctions, false));

  // primary and secondary resp_maps are NULL

  Model vars_transform_model;

  // TODO: this is hardwiring the number of aleatory uncertain
  // aleatory uncertain (normal) variables in the reduced space
  SizetArray recast_vars_comps_total(16, 0);
  recast_vars_comps_total[TOTAL_CAUV] = reducedRank;
  BitArray all_relax_di, all_relax_dr; // default: empty; no discrete relaxation
  const Response& curr_resp = iteratedModel.current_response();
  short recast_resp_order = 1; // recast resp order to be same as original resp
  if (!curr_resp.function_gradients().empty()) recast_resp_order |= 2;
  if (!curr_resp.function_hessians().empty())  recast_resp_order |= 4;

  vars_transform_model.assign_rep(
    new RecastModel(iteratedModel, vars_map_indices, recast_vars_comps_total,
		    all_relax_di, all_relax_dr, nonlinear_vars_map, map_xi_to_x,
		    NULL, primary_resp_map_indices,  secondary_resp_map_indices,
		    recast_secondary_offset, recast_resp_order,
		    nonlinear_resp_map, NULL, NULL), false);

  ParLevLIter pl_iter = methodPCIter->mi_parallel_level_iterator(miPLIndex);
  vars_transform_model.init_communicators(pl_iter, subspaceSamples);

  // convert the normal distributionsto the reduced space and set in the
  // reduced model
  uncertain_vars_to_subspace(iteratedModel, vars_transform_model);

  // Perform UQ on it
  Iterator reduced_space_sampler;

  // convenience function from NonD
  // TODO: additional UQ methods, allow seed
  unsigned short sample_type = SUBMETHOD_RANDOM;

  int mc_seed = (seedSpec) ? seedSpec : generate_system_seed();
  if (seedSpec)
    Cout << "ESM: subspace seed (user-specified) = " << mc_seed << std::endl;
  else
    Cout << "ESM: subspace seed (system-generated) = " << mc_seed << std::endl;

  // might want true for multiple calls...
  bool vary_pattern = false;
  NonD::construct_lhs(reduced_space_sampler, vars_transform_model, sample_type,
		      subspaceSamples, mc_seed, String(), vary_pattern);

  bool all_data = true;
  bool gen_stats = true;
  reduced_space_sampler.sampling_reset(subspaceSamples, all_data, gen_stats);

  reduced_space_sampler.sub_iterator_flag(false);
  reduced_space_sampler.run();//(pl_iter);

  // reduced space UQ results
  Cout << " --- ESM: Results of reduced-space UQ --- \n";
  reduced_space_sampler.print_results(Cout);

  vars_transform_model.free_communicators(pl_iter, subspaceSamples);

}


/** Convert the user-specified normal random variables to the
    appropriate reduced space variables, based on the orthogonal
    transformation.

    TODO: Generalize to convert other random variable types
*/
/// transform and set the distribution parameters in the reduced model
void EfficientSubspaceMethod::
uncertain_vars_to_subspace(Model& native_model,
			   Model& vars_transform_model)
{
  const Pecos::AleatoryDistParams& native_params =
    native_model.aleatory_distribution_parameters();

  // native space characterization
  const RealVector& mu_x = native_params.normal_means();
  const RealVector& sd_x = native_params.normal_std_deviations();
  const RealSymMatrix& correl_x = native_params.uncertain_correlations();

  if (outputLevel >= DEBUG_OUTPUT)
    Cout << "\ncorrel_x = \n" << correl_x;

  bool native_correl = correl_x.empty() ? false : true;
  if (native_correl && correl_x.numRows() != numContinuousVars) {
    Cerr << "Wrong correlation size";
    abort_handler(-1);
  }

  // reduced space characterization: mean mu, std dev sd
  RealVector mu_xi(reducedRank), sd_xi(reducedRank);


  // xi_mu = reducedBasis^T * x_mu
  int m = reducedBasis.numRows();
  int n = reducedBasis.numCols();
  double alpha = 1.0;
  double beta = 0.0;

  int incx = 1;
  int incy = 1;

  // y <-- alpha*A*x + beta*y
  // mu_xi <-- 1.0 * reducedBasis^T * mu_x + 0.0 * mu_xi
  Teuchos::BLAS<int, Real> teuchos_blas;
  teuchos_blas.GEMV(Teuchos::TRANS, m, n, alpha, reducedBasis.values(), m,
		    mu_x.values(), incx, beta, mu_xi.values(), incy);

  // convert the correlations C_x to variance V_x
  // V_x <-- diag(sd_x) * C_x * diag(sd_x)
  // not using symmetric so we can multiply() below
  RealMatrix V_x(reducedBasis.numRows(), reducedBasis.numRows(), false);
  if (native_correl) {
    for (int row=0; row<reducedBasis.numRows(); ++row)
      for (int col=0; col<reducedBasis.numRows(); ++col)
	V_x(row, col) = sd_x(row)*correl_x(row,col)*sd_x(col);
  }
  else {
    V_x = 0.0;
    for (int row=0; row<reducedBasis.numRows(); ++row)
      V_x(row, row) = sd_x(row)*sd_x(row);
  }


  if (outputLevel >= DEBUG_OUTPUT) {
    Cout << "\nreducedBasis = \n" << reducedBasis;
    Cout << "\nV_x =\n" << V_x;
  }

  // compute V_xi = U^T * V_x * U
  alpha = 1.0;
  beta = 0.0;
  RealMatrix UTVx(n, m, false);
  UTVx.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, alpha, reducedBasis, V_x, beta);
  RealMatrix V_xi(reducedRank, reducedRank, false);
  V_xi.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, UTVx, reducedBasis, beta);

  if (outputLevel >= DEBUG_OUTPUT)
    Cout << "\nV_xi = \n" << V_xi;

  // compute the standard deviations in reduced space
  for (int i=0; i<reducedRank; ++i)
    sd_xi = std::sqrt(V_xi(i,i));

  // update the reduced space model
  Pecos::AleatoryDistParams& reduced_dist_params =
    vars_transform_model.aleatory_distribution_parameters();

  reduced_dist_params.normal_means(mu_xi);
  reduced_dist_params.normal_std_deviations(sd_xi);


  // compute the correlations in reduced space
  // TODO: fix symmetric access to not loop over whole matrix
  //  if (native_correl) {

  // Unless the native correl was alpha*I, the reduced variables will
  // be correlated in general, so always set the correltions
  RealSymMatrix correl_xi(reducedRank, false);
  for (int row=0; row<reducedRank; ++row)
    for (int col=0; col<reducedRank; ++col)
      correl_xi(row, col) = V_xi(row,col)/sd_xi(row)/sd_xi(col);

  if (outputLevel >= DEBUG_OUTPUT)
    Cout << "\ncorrel_xi = \n" << correl_xi;

  reduced_dist_params.uncertain_correlations(correl_xi);
    //  }

//   RealVector x_reduced(reducedRank);
//   x_reduced = 1.0;

//   Cout << "x_reduced is\n";
//   Cout << x_reduced << std::endl;

//   vars_transform_model.continuous_variables(x_reduced);
//   vars_transform_model.compute_response();

}



/**
  Perform the variables mapping from recast reduced dimension
  variables xi to original model x variables via linear transformation.
  Maps only continuous variables.
*/
void EfficientSubspaceMethod::
map_xi_to_x(const Variables& recast_xi_vars, Variables& sub_model_x_vars)
{
  Teuchos::BLAS<int, Real> teuchos_blas;

  const RealVector& xi = recast_xi_vars.continuous_variables();
  // TODO: does this yield a view or a copy?
  //RealVector x = sub_model_x_vars.continuous_variables();
  RealVector x;
  copy_data(sub_model_x_vars.continuous_variables(), x);

  //  Calculate x = reducedBasis*xi via matvec directly into x cv in submodel
  //void 	GEMV (ETransp trans, const OrdinalType m, const OrdinalType n, const alpha_type alpha, const A_type *A, const OrdinalType lda, const x_type *x, const OrdinalType incx, const beta_type beta, ScalarType *y, const OrdinalType incy) const
  // 	Performs the matrix-std::vector operation: y <- alpha*A*x+beta*y or y <- alpha*A'*x+beta*y where A is a general m by n matrix.
  const RealMatrix& reduced_basis = esmInstance->reducedBasis;
  int m = reduced_basis.numRows();
  int n = reduced_basis.numCols();

  double alpha = 1.0;
  double beta = 0.0;

  int incx = 1;
  int incy = 1;

  teuchos_blas.GEMV(Teuchos::NO_TRANS, m, n, alpha, reduced_basis.values(), m,
		    xi.values(), incx, beta, x.values(), incy);

  if (esmInstance->outputLevel >= DEBUG_OUTPUT) {
    Cout << "Recast vars are\n";
    Cout << recast_xi_vars << std::endl;
    //    xi.print(Cout);

    Cout << "Submodel vars are\n";
    Cout << sub_model_x_vars << std::endl;
    //x.print(Cout);
  }

  sub_model_x_vars.continuous_variables(x);

}


} // namespace Dakota
