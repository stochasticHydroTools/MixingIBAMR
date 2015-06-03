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

// Config files
#include <IBAMR_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petsc.h>

// Headers for basic SAMRAI objects
#include <PatchLevel.h>
#include <VariableDatabase.h>
#include <tbox/Database.h>
#include <tbox/InputDatabase.h>
#include <tbox/InputManager.h>
#include <tbox/MathUtilities.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/RestartManager.h>
#include <tbox/SAMRAIManager.h>
#include <tbox/SAMRAI_MPI.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// Headers for major algorithm/data structure objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <StandardTagAndInitialize.h>
#include <VisItDataWriter.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffConvectiveOperatorManager.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/AdvDiffStochasticForcing.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/INSStaggeredStochasticForcing.h>
#include <ibamr/RNG.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/CopyToRootSchedule.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include "AdvDiffOperator.h"
#include "BoussinesqForcing.h"

#ifndef DISABLE_ANALYSIS_CODE
extern "C"
{
#include "HydroGrid.h"
}
#endif

namespace  // private namespace
{

// These variables are declared as global variables to facilitate their use in
// callback functions.
Pointer<CellVariable<NDIM,double> > C_var = NULL;
RobinBcCoefStrategy<NDIM>* C_bc_coef = NULL;
Pointer<FaceVariable<NDIM,double> > u_adv_var = NULL;
Pointer<FaceVariable<NDIM,double> > u_s_var = NULL;

// Callback function to update the fluid solver convergence criterion.
void
update_fluid_solver_convergence_criterion(
    const double /*current_time*/,
    const double /*new_time*/,
    const int /*num_cycles*/,
    void* ctx)
{
    INSStaggeredHierarchyIntegrator* fluid_solver = static_cast<INSStaggeredHierarchyIntegrator*>(ctx);
    Pointer<PETScKrylovLinearSolver> stokes_petsc_solver = fluid_solver->getStokesSolver();
    if (stokes_petsc_solver)
    {
        IBAMR_DO_ONCE(pout << "Modifying Stokes    solver convergence criteria to use the minimum initial residual norm for convergence criterion.\n";);
        KSP stokes_ksp = stokes_petsc_solver->getPETScKSP();
        KSPDefaultConvergedSetUMIRNorm(stokes_ksp);
    }
    else
    {
        IBAMR_DO_ONCE(pout << "Warning: Convergence Criteria for Stokes    solver is not robust. It is recommended to run with a PETSc Stokes solver.\n");
    }
    return;
}// update_fluid_solver_convergence_criterion

// Callback function to update the advection-diffucion solver convergence
// criterion.
void
update_adv_diff_solver_convergence_criterion(
    const double /*current_time*/,
    const double /*new_time*/,
    const int /*num_cycles*/,
    void* ctx)
{
    if (C_var)
    {
        AdvDiffHierarchyIntegrator* adv_diff_solver = static_cast<AdvDiffHierarchyIntegrator*>(ctx);
        Pointer<PETScKrylovLinearSolver> adv_diff_petsc_solver = adv_diff_solver->getHelmholtzSolver(C_var);
        if (adv_diff_petsc_solver)
        {
            IBAMR_DO_ONCE(pout << "Modifying Helmholtz solver convergence criteria to use the minimum initial residual norm for convergence criterion.\n");
            KSP adv_diff_ksp = adv_diff_petsc_solver->getPETScKSP();
            KSPDefaultConvergedSetUMIRNorm(adv_diff_ksp);
        }
        else
        {
            IBAMR_DO_ONCE(pout << "Warning: Convergence Criteria for Helmholtz solver is not robust. It is recommended to run with a PETSc Helmholtz solver.\n");
        }
    }
    return;
}// update_adv_diff_solver_convergence_criterion

}

/************************************************************************
 * For each run, the input filename and restart information (if         *
 * needed) must be given on the command line.  For non-restarted case,  *
 * command line is:                                                     *
 *                                                                      *
 *    executable <input file name>                                      *
 *                                                                      *
 * For restarted run, command line is:                                  *
 *                                                                      *
 *    executable <input file name> <restart directory> <restart number> *
 *                                                                      *
 ************************************************************************
 */

int
main(
    int argc,
    char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::setMaxNumberPatchDataEntries(2048);
    SAMRAIManager::startup();

    {// cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Get various application-specific options set in the input file.
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
        const int postprocess_interval = main_db->getIntegerWithDefault("postprocess_interval",0);
        const bool postprocess_data = (postprocess_interval > 0);
        const int snapshot_interval = main_db->getIntegerWithDefault("snapshot_interval",0);
        const int iteration_skip = main_db->getIntegerWithDefault("iteration_skip",0); // How many "equilibration" steps to skip
        const int check_interval = main_db->getIntegerWithDefault("check_interval",0);
        const bool check_data = (check_interval > 0);

        // A flag to control whether two separate HydroGrid analysis are carried
        // out; one for the actual physical domain, and another for the
        // projection along y.
        const bool project_y = main_db->getBoolWithDefault("project_y_axis", false);
        const int project_y_axis = (project_y ? 1 : 0) ;

        // Create major algorithm and data objects which comprise the
        // application.  Each object will be initialized either from input data
        // or restart files, or a combination of both.  Refer to each class
        // constructor for details.  For more information on the composition of
        // objects for this application, see comments at top of file.
        Pointer<INSStaggeredHierarchyIntegrator> fluid_solver = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_solver = new AdvDiffSemiImplicitHierarchyIntegrator(
            "AdvDiffSemiImplicitHierarchyIntegrator", app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>(
            "PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", fluid_solver, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>(
            "LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

        // Create stochastic forcing function specification object.
        Pointer<CartGridFunctionSet> forcing_fcns = new CartGridFunctionSet("forcing_fcns");
        forcing_fcns->addFunction(new INSStaggeredStochasticForcing("INSStaggeredStochasticForcing", app_initializer->getComponentDatabase("INSStaggeredStochasticForcing"), fluid_solver));
        fluid_solver->registerBodyForceFunction(forcing_fcns);

        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction("u_init", input_db->getDatabase("VelocityInitialConditions"), grid_geometry);
        fluid_solver->registerVelocityInitialConditions(u_init);

        Pointer<CartGridFunction> c_init = new muParserCartGridFunction(
            "c_init", app_initializer->getComponentDatabase("ConcentrationInitialConditions"), grid_geometry);

        // Create boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        const bool periodic_domain = periodic_shift.min() != 0;
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_domain)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();
                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();
                u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, input_db->getDatabase(bc_coefs_db_name), grid_geometry);
            }
            fluid_solver->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create a transported quantity and register it with the solver.
        if (input_db->getBoolWithDefault("INCLUDE_CONCENTRATION", true))
        {
            // Create the concentration variable.
            C_var = new CellVariable<NDIM,double>("C",1);
            adv_diff_solver->registerTransportedQuantity(C_var);
            adv_diff_solver->setDiffusionCoefficient(C_var, input_db->getDouble("KAPPA"));
            adv_diff_solver->setInitialConditions(C_var,c_init);

            // Setup the advection velocity.
            const bool include_convective_transport = input_db->getBoolWithDefault("INCLUDE_CONVECTIVE_TRANSPORT", true);
            if (include_convective_transport)
            {
                u_adv_var = new FaceVariable<NDIM,double>("u_adv");
                adv_diff_solver->registerAdvectionVelocity(u_adv_var);
                adv_diff_solver->setAdvectionVelocity(C_var, u_adv_var);
            }

            // Create stochastic forcing function specification objects for the
            // concentration.
            Pointer<CellVariable<NDIM,double> > G_var = new CellVariable<NDIM,double>("G",1);
            adv_diff_solver->registerSourceTerm(G_var);
            adv_diff_solver->setSourceTermFunction(
                G_var, new AdvDiffStochasticForcing("AdvDiffStochasticForcing", app_initializer->getComponentDatabase("AdvDiffStochasticForcing"), C_var, adv_diff_solver));

            // Create boundary condition specification objects (when necessary).
            if (!periodic_domain)
            {
                C_bc_coef = new muParserRobinBcCoefs("C_bc_coef", input_db->getDatabase("ConcentrationBcCoefs"), grid_geometry);
                adv_diff_solver->setPhysicalBcCoef(C_var, C_bc_coef);
            }

            // Create a transported quantity and register it with the solver.
            if (input_db->getBoolWithDefault("INCLUDE_SORET_TRANSPORT", true))
            {
                // Create the velocity variable.
                u_s_var = new FaceVariable<NDIM,double>("u_s");
                adv_diff_solver->registerAdvectionVelocity(u_s_var);
                adv_diff_solver->setAdvectionVelocityFunction(u_s_var, new muParserCartGridFunction("u_s_fcn", input_db->getDatabase("SoretVelocityFunction"), grid_geometry));
                Pointer<KrylovLinearSolver> helmholtz_solver = adv_diff_solver->getHelmholtzSolver(C_var);
                if (!helmholtz_solver)
                {
                    TBOX_ERROR("need to use a Krylov solver for the advection-diffusion solver when Soret terms are enabled.\n");
                }
                Pointer<LaplaceOperator> helmholtz_op = helmholtz_solver->getOperator();
                helmholtz_op->setPhysicalBcCoef(C_bc_coef);
                Pointer<MemoryDatabase> convective_op_db = new MemoryDatabase("");
                convective_op_db->putString("outflow_bdry_extrap_type", "NONE");  // use physical BCs at all boundaries
                Pointer<ConvectiveOperator> convective_op = AdvDiffConvectiveOperatorManager::getManager()->allocateOperator(
                    input_db->getString("SORET_CONVECTIVE_OP"),
                    "soret_convective_op",
                    C_var,
                    convective_op_db,
                    IBAMR::string_to_enum<ConvectiveDifferencingType>(input_db->getString("SORET_DIFFERENCE_FORM")),
                    vector<RobinBcCoefStrategy<NDIM>*>(1, C_bc_coef));
                TimeSteppingType convective_time_stepping_type = IBAMR::string_to_enum<TimeSteppingType>(input_db->getString("SORET_TIME_STEPPING_TYPE"));
                Pointer<AdvDiffOperator> adv_diff_op = new AdvDiffOperator("soret_solver_op", u_s_var, adv_diff_solver, helmholtz_op, convective_op, convective_time_stepping_type, /*homogeneous_bc*/ false);
                helmholtz_solver->setOperator(adv_diff_op);
                adv_diff_solver->setHelmholtzRHSOperator(C_var, adv_diff_op);
            }

            if (C_var && input_db->keyExists("ForcingFunction"))
            {
                Pointer<CellVariable<NDIM,double> > F_var = new CellVariable<NDIM,double>("F",1);
                Pointer<CartGridFunction> F_fcn = new muParserCartGridFunction("F_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
                adv_diff_solver->registerSourceTerm(F_var);
                adv_diff_solver->setSourceTermFunction(F_var, F_fcn);
                adv_diff_solver->setSourceTerm(C_var, F_var);
            }

        }

        if (input_db->getBoolWithDefault("INCLUDE_BOUSSINESQ_FORCING", false))
        {
            forcing_fcns->addFunction(new BoussinesqForcing(C_var, adv_diff_solver, input_db->getDouble("BOUSSINESQ_FORCING_GAMMA")));
        }

        // Seed the random number generator.
        int seed = 0;
        if (input_db->keyExists("SEED"))
        {
            seed = input_db->getInteger("SEED");
        }
        else
        {
            TBOX_ERROR("Key data `seed' not found in input.");
        }
        RNG::parallel_seed(seed);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            fluid_solver->registerVisItDataWriter(visit_data_writer);
            adv_diff_solver->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        fluid_solver->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        if (C_var)
        {
            adv_diff_solver->initializeHierarchyIntegrator(patch_hierarchy, gridding_algorithm);
            // TODO: User code should not have to call these methods explicitly.
            adv_diff_solver->initializeLevelData(patch_hierarchy, /*level_number*/ 0, /*data_time*/ 0.0, /*can_be_refined*/ false, /*initial_time*/ true, /*old_level*/ NULL, /*allocate_data*/ true);
            adv_diff_solver->resetHierarchyConfiguration(patch_hierarchy, /*coarsest_level*/ 0, /*finest_level*/ patch_hierarchy->getFinestLevelNumber());
        }

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Initialize analysis code.
        double dt = fluid_solver->getMaximumTimeStepSize();
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        vector<int> patch_data_idxs;
        patch_data_idxs.push_back(var_db->mapVariableAndContextToIndex(fluid_solver->getVelocityVariable(), fluid_solver->getCurrentContext()));
        if (C_var) patch_data_idxs.push_back(var_db->mapVariableAndContextToIndex(C_var, adv_diff_solver->getCurrentContext()));
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(0);
        static const int ROOT_PROC = 0;
        CopyToRootSchedule copy_to_root_schedule(ROOT_PROC, level, patch_data_idxs);
        const vector<Pointer<PatchData<NDIM> > >& root_patch_data = copy_to_root_schedule.getRootPatchData();

#ifndef DISABLE_ANALYSIS_CODE
        string nmlfile_str = main_db->getStringWithDefault("nml_file", "hydroGridOptions.nml");
        const char* nmlfile = nmlfile_str.c_str();
        if (SAMRAI_MPI::getRank() == ROOT_PROC) setHydroInputFile_C(nmlfile);
#endif
        const Box<NDIM>& domain_box = grid_geometry->getPhysicalDomain()[0];
        int ilower0 = domain_box.lower()(0);
        int iupper0 = domain_box.upper()(0);
        int ilower1 = domain_box.lower()(1);
        int iupper1 = domain_box.upper()(1);
#if (NDIM == 3)
        int ilower2 = domain_box.lower()(2);
        int iupper2 = domain_box.upper()(2);
#endif
        int NX = iupper0 - ilower0 + 1;
        int NY = iupper1 - ilower1 + 1;

        const double* const x_lower = grid_geometry->getXLower();
        const double* const x_upper = grid_geometry->getXUpper();
        double Lx = x_upper[0] - x_lower[0];
        double Ly = x_upper[1] - x_lower[1];
#if (NDIM == 3)
        int NZ = iupper2 - ilower2 + 1;
        double Lz = x_upper[2] - x_lower[2];
#else
        int NZ = 1;
        double Lz = input_db->getDoubleWithDefault("THICKNESS",1.0);
        // NOTE: The thickness should also be used in the stochastic forcing
        // amplitude; this is done manually now.
#endif
        double dx = Lx/NX;
        double dy = Ly/NY;
        double dz = Lz/NZ;
        double dV = dx*dy*dz;
#ifndef DISABLE_ANALYSIS_CODE
        if (postprocess_data)
        {
            int nCells[3] = {NX,NY,NZ};
            double systemLength[3] = { Lx , Ly , Lz };
            double heatCapacity[1] = {1.0};
            double eq_variance = input_db->getDoubleWithDefault("VARIANCE",1.0);
            pout <<  "Grid: dx=" << dx  <<  " dy=" << dy <<  " dz=" << dz <<  " dt=" << dt << "\n";
            if (SAMRAI_MPI::getRank() == ROOT_PROC)
            {
                // We always project along the y axes here
                createHydroAnalysis_C(nCells, 1, NDIM, 1, systemLength, heatCapacity,
                                      dt*postprocess_interval, 0, 1.0/eq_variance, 1);
            }
        }
#endif

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Setup convergence criteria for the solvers.
        fluid_solver->registerPreprocessIntegrateHierarchyCallback(&update_fluid_solver_convergence_criterion, fluid_solver.getPointer());
        adv_diff_solver->registerPreprocessIntegrateHierarchyCallback(&update_adv_diff_solver_convergence_criterion, adv_diff_solver.getPointer());

        // Setup convergence criteria for the Helmholtz solver.

        // Prepare to compute momentum and variances.
        int n_checks = 0;
        double mom0, mom1, mom2, conc;
        double mean0=0.0, mean1=0.0, mean2=0.0, meanc=0.0;
        double var00, var01, var11;
        double var00_mean=0.0, var01_mean=0.0, var11_mean=0.0;

        // Write out initial visualization data.
        int iteration_num = fluid_solver->getIntegratorStep();
        double loop_time = fluid_solver->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            fluid_solver->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
        }

        // Main time step loop.
        double loop_time_end = fluid_solver->getEndTime();
        while (!MathUtilities<double>::equalEps(loop_time,loop_time_end) && fluid_solver->stepsRemaining())
        {
            iteration_num = fluid_solver->getIntegratorStep();
            loop_time = fluid_solver->getIntegratorTime();
            double dt = fluid_solver->getMaximumTimeStepSize();
            if (C_var) dt = min(dt, adv_diff_solver->getMaximumTimeStepSize());
#if 0
            pout << "step=" << iteration_num << " t=" << loop_time << " dt=" << dt << "\n";
            pout <<                                                    "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " <<  iteration_num << "\n";
            pout << "Simulation time is " << loop_time << " dt=" << dt << "\n";
#endif
            if (postprocess_data && (iteration_num>=iteration_skip) && ((iteration_num-iteration_skip)%postprocess_interval == 0))
            {
                copy_to_root_schedule.communicate();
#ifndef DISABLE_ANALYSIS_CODE
                if (SAMRAI_MPI::getRank() == ROOT_PROC)
                {
                    Pointer<SideData<NDIM,double> > u_data = root_patch_data[0];
                    int u_gcw = u_data->getGhostCellWidth().max();
                    TBOX_ASSERT(u_gcw == 0); // The analysis code assumes the velocity data has no ghost cells.

                    double* u0 = u_data->getPointer(0);
                    double* u1 = u_data->getPointer(1);
#if (NDIM == 3)
                    double* u2 = u_data->getPointer(2);
#endif
                    Pointer<CellData<NDIM,double> > C_data = (C_var ? root_patch_data[1] : Pointer<PatchData<NDIM> >(NULL));
                    int nScalars = (C_var ? C_data->getDepth() : 0);
                    int C_gcw = (C_var ? C_data->getGhostCellWidth().max() : 0);
                    double* C = (C_var ? C_data->getPointer() : NULL);
#if (NDIM == 2)
                    updateHydroAnalysisStaggered_C (u_gcw+1, nScalars, u0, u1, NULL, C_gcw, C);
#endif
#if (NDIM == 3)
                    updateHydroAnalysisStaggered_C (u_gcw+1, nScalars, u0, u1, u2, C_gcw, C);
#endif
                    if (project_y_axis && C_var)
                    {
                        // There is no density and concentration here (it is constant), so we pass C twice
                        // Passing -2 here means do analysis of both full data and projection but not save the full grid
                        projectHydroGrid_C(C, C, "", (iteration_num-iteration_skip)/postprocess_interval, -2);
                    }
                    if ((snapshot_interval > 0) && C_var)
                    {
                        if ((iteration_num-iteration_skip)%snapshot_interval == 0)
                        {
                            writeHydroGridMixture_C (C, C, "", (iteration_num-iteration_skip)/snapshot_interval);
                            if(!project_y_axis) writeToFiles_C((iteration_num-iteration_skip)/snapshot_interval); // Write statistics up to now
                        }
                    }
                }
#endif
            }

            // Here we take control over the time stepping scheme (instead of
            // relying on HierarchyIntegrator::advanceHierarchy()) in order to
            // control the interlacing of the various solves.
            const double current_time = loop_time;
            const double new_time = current_time + dt;
            const int num_cycles = 2;
            fluid_solver->preprocessIntegrateHierarchy(current_time, new_time, num_cycles);
            adv_diff_solver->preprocessIntegrateHierarchy(current_time, new_time, num_cycles);
           
            for (int cycle_num = 0; cycle_num < num_cycles; ++cycle_num)
            {
                // Solve for the fluid velocity.
                fluid_solver->integrateHierarchy(current_time, new_time, cycle_num);
            
                if (C_var)
                {
                    // Copy the updated velocity into the current, scratch, and
                    // new velocity variables used by the concentration solver.
                    // This only has the desired effect when we use
                    // MIDPOINT_RULE time stepping for the convective term.
                    TBOX_ASSERT(adv_diff_solver->getConvectiveTimeSteppingType(C_var) == MIDPOINT_RULE);
                    if (u_adv_var)
                    {
                        for (int level_num = 0; level_num <= patch_hierarchy->getFinestLevelNumber(); ++level_num)
                        {
                            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(level_num);
                            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                            {
                                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                                Pointer<FaceData<NDIM,double> > u_adv_current_data = patch->getPatchData(u_adv_var, adv_diff_solver->getCurrentContext());
                                Pointer<FaceData<NDIM,double> > u_adv_scratch_data = patch->getPatchData(u_adv_var, adv_diff_solver->getScratchContext());
                                Pointer<FaceData<NDIM,double> > u_adv_new_data = patch->getPatchData(u_adv_var, adv_diff_solver->getNewContext());
                                Pointer<SideData<NDIM,double> > u_new_data = patch->getPatchData(fluid_solver->getVelocityVariable(), fluid_solver->getNewContext());
                                for (int axis = 0; axis < NDIM; ++axis)
                                {
                                    for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch->getBox(),axis)); b; b++)
                                    {
                                        const Index<NDIM>& i = b();
                                        const SideIndex<NDIM> i_s(i, axis, SideIndex<NDIM>::Lower);
                                        const FaceIndex<NDIM> i_f(i, axis, FaceIndex<NDIM>::Lower);
                                        (*u_adv_current_data)(i_f) = (*u_new_data)(i_s);
                                        (*u_adv_scratch_data)(i_f) = (*u_new_data)(i_s);
                                        (*u_adv_new_data)(i_f) = (*u_new_data)(i_s);
                                    }
                                }
                            }
                        }
                    }

                    // Solve for the concentration, using the most recently
                    // computed velocity as the advection velocity.
                    adv_diff_solver->integrateHierarchy(current_time, new_time, cycle_num);
                }
            }
                        
            // Execute postprocessing operations.
            static const bool skip_synchronize_new_state_data = true;
            fluid_solver->postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data, num_cycles);
            fluid_solver->synchronizeHierarchyData(NEW_DATA);
            fluid_solver->resetTimeDependentHierarchyData(new_time);
            if (C_var)
            {
                adv_diff_solver->postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data, num_cycles);
                adv_diff_solver->synchronizeHierarchyData(NEW_DATA);
                adv_diff_solver->resetTimeDependentHierarchyData(new_time);
            }
            
            loop_time += dt;
#if 0
            pout <<                                                    "\n";
            pout << "At end       of timestep # " <<  iteration_num << "\n";
            pout << "Simulation time is " << loop_time              << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout <<                                                    "\n";
#endif
            // At specified intervals, write visualization and restart files,
            // print out timer data, and perform postprocessing.
            iteration_num += 1;
            const bool last_step = !fluid_solver->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num%viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                fluid_solver->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num%restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num%timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (check_data && (iteration_num>=iteration_skip) && ((iteration_num-iteration_skip)%check_interval == 0))
            {
                pout << "At timestep # " <<  iteration_num << " t=" << loop_time << " dt=" << dt << "\n";
                n_checks += 1;

                // Calculate moments.
                mom0 = 0.0; mom1 = 0.0; mom2 = 0.0;
                conc = 0.0;
                var00 = 0.0; var01 = 0.0; var11 = 0.0;
                copy_to_root_schedule.communicate();
                if (SAMRAI_MPI::getRank() == ROOT_PROC)
                {
                    Pointer<SideData<NDIM,double> > u_data = root_patch_data[0];
                    for (BoxIterator<NDIM> b(u_data->getBox()); b; b++)
                    {
                        const Index<NDIM>& i = b();
                        const SideIndex<NDIM> i_x_lower(i,0,0);
                        const SideIndex<NDIM> i_x_upper(i,0,1);
                        const SideIndex<NDIM> i_y_lower(i,1,0);
                        const SideIndex<NDIM> i_y_upper(i,1,1);
                        mom0  += (*u_data)(i_x_lower);
                        mom1  += (*u_data)(i_y_lower);
                        var00 += (*u_data)(i_x_lower)*(*u_data)(i_x_lower);
                        var11 += (*u_data)(i_y_lower)*(*u_data)(i_y_lower);
#if (NDIM == 3)
                        const SideIndex<NDIM> i_z_lower(i,2,0);
                        const SideIndex<NDIM> i_z_upper(i,2,1);
                        mom2 += (*u_data)(i_z_lower);
#endif
                    }
                    if (C_var)
                    {
                        Pointer<CellData<NDIM,double> > C_data = (C_var ? root_patch_data[1] : Pointer<PatchData<NDIM> >(NULL));
                        for (BoxIterator<NDIM> b(C_data->getBox()); b; b++)
                        {
                            const Index<NDIM>& i = b();
                            conc  += (*C_data)(i);
                            var01 += (*C_data)(i)*(*C_data)(i);
                        }
                    }
                }

                MPI_Bcast(& mom0, 1, MPI_DOUBLE, ROOT_PROC, SAMRAI_MPI::getCommunicator());
                MPI_Bcast(& mom1, 1, MPI_DOUBLE, ROOT_PROC, SAMRAI_MPI::getCommunicator());
                MPI_Bcast(& mom2, 1, MPI_DOUBLE, ROOT_PROC, SAMRAI_MPI::getCommunicator());
                MPI_Bcast(&var00, 1, MPI_DOUBLE, ROOT_PROC, SAMRAI_MPI::getCommunicator());
                MPI_Bcast(&var11, 1, MPI_DOUBLE, ROOT_PROC, SAMRAI_MPI::getCommunicator());
                MPI_Bcast(& conc, 1, MPI_DOUBLE, ROOT_PROC, SAMRAI_MPI::getCommunicator());
                MPI_Bcast(&var01, 1, MPI_DOUBLE, ROOT_PROC, SAMRAI_MPI::getCommunicator());

                mean0 += mom0;
                mean1 += mom1;
                mean2 += mom2;
                meanc += conc;
                var00_mean += var00;
                var11_mean += var11;
                var01_mean += var01;

                mom0  /= (NX*NY*NZ);
                mom1  /= (NX*NY*NZ);
                mom2  /= (NX*NY*NZ);
                conc  /= (NX*NY*NZ);
                var00 /= (NX*NY*NZ);
                var11 /= (NX*NY*NZ);
                var01 /= (NX*NY*NZ);

                var00 = (var00 - mom0*mom0) * dV;
                var11 = (var11 - mom1*mom1) * dV;
                var01 = (var01 - conc*conc) * dV;

                // Verify conservation:
                pout << iteration_num << " <px>=" << mom0 << " <py>=" << mom1 << " <pz>=" << mom2 << " <c>=" << conc << endl;

                // Monitor variances:
                pout << iteration_num << " <dvx^2>=" << var00 << " <dvy^2>=" << var11 << " <dc^2>=" << var01 << endl;
                pout << "t=" << loop_time << " <cfl_x>=" << sqrt(var00)*dt/dx << " <cfl_y>=" << sqrt(var11)*dt/dy << endl;

            }
        } // end while loop for stepping

        mean0 /= n_checks*(NX*NY*NZ);
        mean1 /= n_checks*(NX*NY*NZ);
        mean2 /= n_checks*(NX*NY*NZ);
        meanc /= n_checks*(NX*NY*NZ);
        var00_mean /= n_checks*(NX*NY*NZ);
        var11_mean /= n_checks*(NX*NY*NZ);
        var01_mean /= n_checks*(NX*NY*NZ);

        var00_mean = (var00_mean - mean0*mean0) * dV;
        var11_mean = (var11_mean - mean1*mean1) * dV;
        var01_mean = (var01_mean - meanc*meanc) * dV;

        pout << "Final <px>=" << mom0 << " <py>=" << mom1 << " <c>=" << meanc << endl;
        pout << "Final <dvx^2>=" << var00_mean << " <dvy^2>=" << var11_mean << " <dc^2>=" << var01_mean << endl;

        // Print out the results of the analysis code and clean up
        if (postprocess_data && (SAMRAI_MPI::getRank() == ROOT_PROC))
        {
#ifndef DISABLE_ANALYSIS_CODE
            if(!project_y_axis) writeToFiles_C(-1); // Write to files
            destroyHydroAnalysis_C();
#endif
        }

        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        delete C_bc_coef;

    }// cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
}// main
