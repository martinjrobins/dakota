/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014 Sandia Corporation.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

//- Class:        TestDriverInterface
//- Description:  Direct interfaces to test drivers and "simple" linked
//-               applications that don't require separate setup and tear-down
//- Owner:        Mike Eldred, Brian Adams
//- Version: $Id$

#ifndef MARTINJROBINS_DRIVER_INTERFACE_H
#define MARTINJROBINS_DRIVER_INTERFACE_H

#include "DirectApplicInterface.hpp"

namespace Dakota {

/** Specialization of DirectApplicInterface to embed algebraic test function
    drivers directly in Dakota */
class MartinJRobinsDriverInterface: public DirectApplicInterface
{
public:

  //
  //- Heading: Constructor and destructor
  //

  MartinJRobinsDriverInterface(const ProblemDescDB& problem_db); ///< constructor
  ~MartinJRobinsDriverInterface();                               ///< destructor

protected:

  //
  //- Heading: Virtual function redefinitions
  //

  /// execute an analysis code portion of a direct evaluation invocation
  virtual int derived_map_ac(const Dakota::String& ac_name);

private:

  //
  //- Heading: Simulators and test functions
  //

  int e_surface_driver();
  int e_solution_driver();

};

} // namespace Dakota

#endif  // MARTINJROBINS_DRIVER_INTERFACE_HPP
