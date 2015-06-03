// Filename: AdvDiffOperator.C
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

#include "AdvDiffOperator.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <SAMRAIVectorReal.h>

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffOperator::AdvDiffOperator(
    const string& object_name,
    Pointer<Variable<NDIM> > u_var,
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_solver,
    Pointer<LaplaceOperator> laplace_op,
    Pointer<ConvectiveOperator> convective_op,
    TimeSteppingType convective_time_stepping_type,
    const bool homogeneous_bc)
    : LaplaceOperator(object_name, homogeneous_bc),
      d_u_var(u_var),
      d_adv_diff_solver(adv_diff_solver),
      d_laplace_op(laplace_op),
      d_convective_op(convective_op),
      d_convective_time_stepping_type(convective_time_stepping_type)
{
    // Setup the operator to use default scalar-valued boundary conditions.
    setHomogeneousBc(homogeneous_bc);
    setPhysicalBcCoef(NULL);
    return;
}// AdvDiffOperator()

AdvDiffOperator::~AdvDiffOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
}// ~AdvDiffOperator()

void
AdvDiffOperator::setHomogeneousBc(
    bool homogeneous_bc)
{
    LaplaceOperator::setHomogeneousBc(homogeneous_bc);
    d_laplace_op->setHomogeneousBc(getHomogeneousBc());
    d_convective_op->setHomogeneousBc(getHomogeneousBc());
    return;
}// setHomogeneousBc

void
AdvDiffOperator::setSolutionTime(
    double solution_time)
{
    LaplaceOperator::setSolutionTime(solution_time);
    d_laplace_op->setSolutionTime(getSolutionTime());
    d_convective_op->setSolutionTime(getSolutionTime());
    return;
}// setSolutionTime

void
AdvDiffOperator::setTimeInterval(
    double current_time,
    double new_time)
{
    LaplaceOperator::setTimeInterval(current_time, new_time);
    d_laplace_op->setTimeInterval(getTimeInterval().first, getTimeInterval().second);
    d_convective_op->setTimeInterval(getTimeInterval().first, getTimeInterval().second);
    return;
}// setTimeInterval

void
AdvDiffOperator::setHierarchyMathOps(
    Pointer<HierarchyMathOps> hier_math_ops)
{
    LaplaceOperator::setHierarchyMathOps(hier_math_ops);
    d_laplace_op->setHierarchyMathOps(getHierarchyMathOps());
    d_convective_op->setHierarchyMathOps(getHierarchyMathOps());
    return;
}// setHierarchyMathOps

void
AdvDiffOperator::setPoissonSpecifications(
    const PoissonSpecifications& poisson_spec)
{
    LaplaceOperator::setPoissonSpecifications(poisson_spec);
    d_laplace_op->setPoissonSpecifications(getPoissonSpecifications());
    return;
}// setPoissonSpecifications

void
AdvDiffOperator::setPhysicalBcCoef(
    RobinBcCoefStrategy<NDIM>* bc_coef)
{
    LaplaceOperator::setPhysicalBcCoef(bc_coef);
    d_laplace_op->setPhysicalBcCoefs(getPhysicalBcCoefs());
    return;
}// setPhysicalBcCoef

void
AdvDiffOperator::setPhysicalBcCoefs(
    const vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
    LaplaceOperator::setPhysicalBcCoefs(bc_coefs);
    d_laplace_op->setPhysicalBcCoefs(getPhysicalBcCoefs());
    return;
}// setPhysicalBcCoefs

void
AdvDiffOperator::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    bool evaluating_rhs = MathUtilities<double>::equalEps(getSolutionTime(), getTimeInterval().first);
    double convective_wgt = 0.0;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int u_idx = -1;
    switch (d_convective_time_stepping_type)
    {
        case FORWARD_EULER:
            convective_wgt = evaluating_rhs ? -1.0 :  0.0;
            u_idx = var_db->mapVariableAndContextToIndex(d_u_var, d_adv_diff_solver->getCurrentContext());
            break;
        case BACKWARD_EULER:
            convective_wgt = evaluating_rhs ?  0.0 : +1.0;
            u_idx = var_db->mapVariableAndContextToIndex(d_u_var, d_adv_diff_solver->getNewContext());
            break;
        case MIDPOINT_RULE:
            if (d_adv_diff_solver->getCurrentCycleNumber() > 0)
            {
                TBOX_ERROR(d_object_name << "::apply():\n"
                       << "  current implementation does not work correctly for MIDPOINT_RULE time stepping and NUM_CYCLES > 1\n");
            }
            u_idx = var_db->mapVariableAndContextToIndex(d_u_var, d_adv_diff_solver->getScratchContext());
            convective_wgt = evaluating_rhs ? -0.5 : +0.5;
            break;
        case TRAPEZOIDAL_RULE:
            convective_wgt = evaluating_rhs ? -0.5 : +0.5;
            u_idx = var_db->mapVariableAndContextToIndex(d_u_var, evaluating_rhs ? d_adv_diff_solver->getCurrentContext() : d_adv_diff_solver->getNewContext());
            break;
        default:
            TBOX_ERROR(d_object_name << "::apply():\n"
                       << "  unsupported convective time stepping type: " << IBAMR::enum_to_string<TimeSteppingType>(d_convective_time_stepping_type) << " \n"
                       << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    if (convective_wgt != 0.0)
    {
        d_convective_op->setAdvectionVelocity(u_idx);
        d_convective_op->apply(x,y);
        y.scale(convective_wgt, Pointer<SAMRAIVectorReal<NDIM,double> >(&y,false));
        d_laplace_op->applyAdd(x, y, y);
    }
    else
    {
        d_laplace_op->apply(x, y);
    }
    return;
}// apply

void
AdvDiffOperator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    d_laplace_op->initializeOperatorState(in,out);
    d_convective_op->initializeOperatorState(in,out);
    return;
}// initializeOperatorState

void
AdvDiffOperator::deallocateOperatorState()
{
    d_laplace_op->deallocateOperatorState();
    d_convective_op->deallocateOperatorState();
    return;
}// deallocateOperatorState

//////////////////////////////////////////////////////////////////////////////
