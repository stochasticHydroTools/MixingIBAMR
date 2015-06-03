// Filename: AdvDiffOperator.h
// Created on 02 Feb 2013 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_AdvDiffOperator
#define included_AdvDiffOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/ConvectiveOperator.h>
#include <ibamr/app_namespaces.h>

// IBTK INCLUDES
#include <ibtk/LaplaceOperator.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

/*!
 * \brief Class AdvDiffOperator is a concrete LaplaceOperator implementing a
 * linear advection-diffusion operator A of the form A*c = L*c + N*c, in which L
 * is a Laplace-type operator, and N is a convective-type operator.
 */
class AdvDiffOperator
    : public LaplaceOperator
{
public:
    /*!
     * \brief Constructor for class AdvDiffOperator initializes the operator
     * coefficients and boundary conditions to default values.
     */
    AdvDiffOperator(
        const string& object_name,
        Pointer<Variable<NDIM> > u_var,
        Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
        Pointer<LaplaceOperator> laplace_op,
        Pointer<ConvectiveOperator> convective_op,
        TimeSteppingType convective_time_stepping_type,
        bool homogeneous_bc=true);

    /*!
     * \brief Destructor.
     */
    ~AdvDiffOperator();

    /*!
     * \brief Set whether the operator should use homogeneous boundary conditions.
     */
    void
    setHomogeneousBc(
        bool homogeneous_bc);

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    void
    setSolutionTime(
        double solution_time);

    /*!
     * \brief Set the current time interval.
     */
    void
    setTimeInterval(
        double current_time,
        double new_time);

    /*!
     * \brief Set the HierarchyMathOps object used by the operator.
     */
    void
    setHierarchyMathOps(
        Pointer<HierarchyMathOps> hier_math_ops);

    /*!
     * \brief Set the PoissonSpecifications object used to specify
     * the coefficients for the Laplace operator.
     */
    void
    setPoissonSpecifications(
        const PoissonSpecifications& poisson_spec);

    /*!
     * \brief Set the RobinBcCoefStrategy object used to specify physical
     * boundary conditions.
     *
     * \note \a bc_coef may be NULL.  In this case, default boundary conditions
     * (as supplied to the class constructor) are employed.
     *
     * \param bc_coef  Pointer to an object that can set the Robin boundary condition coefficients
     */
    void
    setPhysicalBcCoef(
        RobinBcCoefStrategy<NDIM>* bc_coef);

    /*!
     * \brief Set the RobinBcCoefStrategy objects used to specify physical
     * boundary conditions.
     *
     * \note Any of the elements of \a bc_coefs may be NULL.  In this case,
     * default boundary conditions (as supplied to the class constructor) are
     * employed for that data depth.
     *
     * \param bc_coefs  Vector of pointers to objects that can set the Robin boundary condition coefficients
     */
    void
    setPhysicalBcCoefs(
        const vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * \name Linear operator functionality.
     */
    //\{

    /*!
     * \brief Compute y=Ax.
     *
     * Before calling this function, the form of the vectors x and y should be
     * set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * Thus, the user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must
     * have enough ghost cells for computation of Ax).
     *
     * \note In general, the vectors x and y \em cannot be the same.
     *
     * Upon return from this function, the y vector will contain the result of
     * the application of A to x.
     *
     * initializeOperatorState must be called prior to any calls to
     * applyOperator.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y output: y=Ax
     */
    void
    apply(
        SAMRAIVectorReal<NDIM,double>& x,
        SAMRAIVectorReal<NDIM,double>& y);

    /*!
     * \brief Compute hierarchy-dependent data required for computing y=Ax (and
     * y=A'x).
     *
     * \param in input vector
     * \param out output vector
     *
     * \see KrylovLinearSolver::initializeSolverState
     */
    void
    initializeOperatorState(
        const SAMRAIVectorReal<NDIM,double>& in,
        const SAMRAIVectorReal<NDIM,double>& out);

    /*!
     * \brief Remove all hierarchy-dependent data computed by
     * initializeOperatorState().
     *
     * Remove all hierarchy-dependent data set by initializeOperatorState().  It
     * is safe to call deallocateOperatorState() even if the state is already
     * deallocated.
     *
     * \see initializeOperatorState
     * \see KrylovLinearSolver::deallocateSolverState
     */
    void
    deallocateOperatorState();

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffOperator(
        const AdvDiffOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffOperator&
    operator=(
        const AdvDiffOperator& that);

    Pointer<Variable<NDIM> > d_u_var;
    Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_solver;
    Pointer<LaplaceOperator> d_laplace_op;
    Pointer<ConvectiveOperator> d_convective_op;
    TimeSteppingType d_convective_time_stepping_type;
};

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/AdvDiffOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffOperator
