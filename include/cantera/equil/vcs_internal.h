/**
 *  @file vcs_internal.h Internal declarations for the VCSnonideal package
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef _VCS_INTERNAL_H
#define _VCS_INTERNAL_H

#include "cantera/base/global.h"
namespace Cantera
{
//! define this Cantera function to replace printf
/*!
 * We can replace this with printf easily
 */
#define plogf writelogf

//! define this Cantera function to replace cout << endl;
/*!
 * We use this to place an endl in the log file, and ensure that the IO buffers
 * are flushed.
 */
#define plogendl() writelogendl()

//! Global hook for turning on and off time printing.
/*!
 * Default is to allow printing. But, you can assign this to zero globally to
 * turn off all time printing. This is helpful for test suite purposes where you
 * are interested in differences in text files.
 */
extern int vcs_timing_print_lvl;

// Forward references
class VCS_SPECIES_THERMO;
class VCS_PROB;

//! Class to keep track of time and iterations
/*!
 * class keeps all of the counters together.
 */
class VCS_COUNTERS
{
public:
    //! Total number of iterations in the main loop
    //! of vcs_TP() to solve for thermo equilibrium
    int T_Its;

    //! Current number of iterations in the main loop
    //! of vcs_TP() to solve for thermo equilibrium
    int Its;

    //! Total number of optimizations of the components basis set done
    int T_Basis_Opts;

    //! number of optimizations of the components basis set done
    int Basis_Opts;

    //! Current number of times the initial thermo equilibrium estimator has
    //! been called
    int T_Calls_Inest;

    //! Current number of calls to vcs_TP
    int T_Calls_vcs_TP;

    //! Current time spent in vcs_TP
    double T_Time_vcs_TP;

    //! Current time spent in vcs_TP
    double Time_vcs_TP;

    //! Total Time spent in basopt
    double T_Time_basopt;

    //! Current Time spent in basopt
    double Time_basopt;

    //! Time spent in initial estimator
    double T_Time_inest;

    //! Time spent in the vcs suite of programs
    double T_Time_vcs;
};

//! Returns the value of the gas constant in the units specified by parameter
/*!
 *  @param mu_units Specifies the units.
 *           -  VCS_UNITS_KCALMOL: kcal gmol-1 K-1
 *           -  VCS_UNITS_UNITLESS:  1.0 K-1
 *           -  VCS_UNITS_KJMOL:   kJ gmol-1 K-1
 *           -  VCS_UNITS_KELVIN:    1.0 K-1
 *           -  VCS_UNITS_MKS:   joules kmol-1 K-1 =  kg m2 s-2 kmol-1 K-1
 */
double vcsUtil_gasConstant(int mu_units);

//! Definition of the function pointer for the root finder
/*!
 *  see vcsUtil_root1d for a definition of how to use this.
 */
typedef double(*VCS_FUNC_PTR)(double xval, double Vtarget,
                              int varID, void* fptrPassthrough,
                              int* err);

//! determine the l2 norm of a vector of doubles
/*!
 * @param vec vector of doubles
 * @return  Returns the l2 norm of the vector
 */
double vcs_l2norm(const vector_fp vec);

//! Finds the location of the maximum component in a double vector
/*!
 * @param x pointer to a vector of doubles
 * @param xSize pointer to a vector of doubles used as a multiplier to x[]
 *              before making the decision. Ignored if set to NULL.
 * @param j lowest index to search from
 * @param n highest index to search from
 * @returns index of the greatest value on X(i) searched, j <= i < n
 */
size_t vcs_optMax(const double* x, const double* xSize, size_t j, size_t n);

//! Returns a const char string representing the type of the species given by
//! the first argument
/*!
 * @param speciesStatus  Species status integer representing the type
 *                       of the species.
 * @param length         Maximum length of the string to be returned.
 *                       Shorter values will yield abbreviated strings.
 *                       Defaults to a value of 100.
 */
const char* vcs_speciesType_string(int speciesStatus, int length = 100);

//! Print a string within a given space limit
/*!
 * This routine limits the amount of the string that will be printed to a
 * maximum of "space" characters. Printing is done to to Cantera's writelog()
 * function.
 *
 * @param str  String, which must be null terminated.
 * @param space   space limit for the printing.
 * @param alignment Alignment of string within the space:
 *                 -  0 centered
 *                 -  1 right aligned
 *                 -  2 left aligned
 */
void vcs_print_stringTrunc(const char* str, size_t space, int alignment);

//! Simple routine to check whether two doubles are equal up to roundoff error
/*!
 * Currently it's set to check for 10 digits of relative accuracy.
 *
 * @param d1 first double
 * @param d2 second double
 *
 * @returns true if the doubles are "equal" and false otherwise
 */
bool vcs_doubleEqual(double d1, double d2);
}

#endif