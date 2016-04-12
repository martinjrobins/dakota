/*  _______________________________________________________________________

    PECOS: Parallel Environment for Creation Of Stochastics
    Copyright (c) 2011, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Pecos directory.
    _______________________________________________________________________ */

//- Class:        InterpPolyApproximation
//- Description:  Class for Lagrange Interpolation Polynomial Approximation
//-               
//- Owner:        Mike Eldred

#ifndef INTERP_POLY_APPROXIMATION_HPP
#define INTERP_POLY_APPROXIMATION_HPP

#include "PolynomialApproximation.hpp"
#include "BasisPolynomial.hpp"

namespace Pecos {


/// Derived approximation class for interpolation polynomials (global
/// approximation).

/** The InterpPolyApproximation class provides a global approximation
    based on interpolation polynomials.  It is used primarily for
    stochastic collocation approaches to uncertainty quantification. */

class InterpPolyApproximation: public PolynomialApproximation
{
public:

  //
  //- Heading: Constructor and destructor
  //

  /// default constructor
  InterpPolyApproximation(const SharedBasisApproxData& shared_data);
  /// destructor
  ~InterpPolyApproximation();

  //
  //- Heading: member functions
  //

  //
  //- Heading: Virtual function redefinitions
  //

  /// compute the coefficients for the expansion of multivariate
  /// interpolation polynomials
  void compute_coefficients();

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  int min_coefficients() const;

  /// update the coefficients for the expansion of multivariate Lagrange
  /// interpolation polynomials
  void increment_coefficients();
  /// restore the coefficients to their previous state prior to last increment
  void decrement_coefficients();
  /// restore the coefficients to a previously incremented state as
  /// identified by the current increment to the Smolyak multi index
  void restore_coefficients();
  /// finalize the coefficients by applying all previously evaluated increments
  void finalize_coefficients();

  /// size expansionType{1,2}Coeffs and expansionType1CoeffGrads
  void allocate_arrays();

  /// computes component (main and interaction) Sobol' indices
  void compute_component_sobol();
  /// computes total Sobol' indices
  void compute_total_sobol();

  /// compute numerical moments to order 4
  void compute_moments();
  /// compute numerical moments in all-variables mode to order 2
  void compute_moments(const RealVector& x);

  //
  //- Heading: New virtual functions
  //

  /// derived portion of allocate_arrays()
  virtual void allocate_expansion_coefficients() = 0;
  /// derived portion of compute_coefficients()
  virtual void compute_expansion_coefficients() = 0;
  /// increment expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  /// within increment_coefficients()
  virtual void increment_expansion_coefficients() = 0;
  /// decrement expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  /// within decrement_coefficients()
  virtual void decrement_expansion_coefficients() = 0;
  /// restore expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  /// within restore_coefficients()
  virtual void restore_expansion_coefficients() = 0;
  /// finalize expansion{Type1Coeffs,Type2Coeffs,Type1CoeffGrads}
  /// within finalize_coefficients()
  virtual void finalize_expansion_coefficients() = 0;

  /// compute moments of response using numerical integration
  virtual void integrate_response_moments(size_t num_moments) = 0;
  /// compute moments of expansion using numerical integration
  virtual void integrate_expansion_moments(size_t num_moments) = 0;

  virtual void compute_total_sobol_indices() = 0;
  virtual void compute_partial_variance(const BitArray& set_value);

  //
  //- Heading: Convenience functions
  //

  //
  //- Heading: Data
  //

  /// the partialVariances of subset functions f_alpha
  RealVector partialVariance;

private:

  //
  //- Heading: Convenience functions
  //

  /// recursively identifies constituent subsets that are children of
  /// a parent set
  void proper_subsets(const BitArray& parent_set, BitArraySet& children);

  //
  //- Heading: Data
  //
};


inline InterpPolyApproximation::
InterpPolyApproximation(const SharedBasisApproxData& shared_data):
  PolynomialApproximation(shared_data)
{ }


inline InterpPolyApproximation::~InterpPolyApproximation()
{ }


inline void InterpPolyApproximation::compute_moments()
{
  // standard variables mode supports four moments using the collocation rules
  // as integration rules
  integrate_response_moments(4);

  // do this second so that clearing any existing rules does not cause rework
  //if (expConfigOptions.outputLevel >= VERBOSE_OUTPUT)
    integrate_expansion_moments(4);
}


inline void InterpPolyApproximation::compute_moments(const RealVector& x)
{
  // all variables mode only supports first two moments
  mean(x); variance(x);
  //standardize_moments(numericalMoments);
  //integrate_expansion_moments(4, x);

  // Note: it would be feasible to implement an all_variables version of
  // integrate_expansion_moments() by evaluating the combined expansion at
  // {design/epistemic=initialPtU,aleatory=Gauss points}
  // --> can't do this for integrate_response_moments() (lacking reqd resp data)
  // --> would require generation of new TPQ/SSG grid only over aleatory vars
  // --> could possibly retire redundant all_vars functions
}

} // namespace Pecos

#endif
