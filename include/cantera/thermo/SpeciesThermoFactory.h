/**
 *  @file SpeciesThermoFactory.h
 *    Header for factory functions to build instances of classes that manage the
 *    standard-state thermodynamic properties of a set of species
 *    (see \ref spthermo);
 */
// Copyright 2001  California Institute of Technology

#ifndef SPECIESTHERMO_FACTORY_H
#define SPECIESTHERMO_FACTORY_H

#include "SpeciesThermoInterpType.h"

namespace Cantera
{

class XML_Node;

//! Create a new SpeciesThermoInterpType object given a corresponding constant.
/*!
 *  @param type A constant specifying the type to be created
 *  @param tlow The lowest temperature at which the parameterization is valid
 *  @param thigh The highest temperature at which the parameterization is valid
 *  @param pref The reference pressure for the parameterization
 *  @param coeffs The array of coefficients for the parameterization
 *  @return       Returns the pointer to the newly allocated
 *                SpeciesThermoInterpType object
 */
SpeciesThermoInterpType* newSpeciesThermoInterpType(int type, double tlow,
    double thigh, double pref, const double* coeffs);

//! Create a new SpeciesThermoInterpType object given a string
/*!
 *  @param type String name for the species thermo type
 *  @param tlow The lowest temperature at which the parameterization is valid
 *  @param thigh The highest temperature at which the parameterization is valid
 *  @param pref The reference pressure for the parameterization
 *  @param coeffs The array of coefficients for the parameterization
 *  @return       Returns the pointer to the newly allocated
 *                SpeciesThermoInterpType object
 */
SpeciesThermoInterpType* newSpeciesThermoInterpType(const std::string& type,
    double tlow, double thigh, double pref, const double* coeffs);

//! Create a new SpeciesThermoInterpType object from XML_Node
/*!
 *  @param thermoNode 'thermo' XML_Node (child of the 'species' node) with child
 *      nodes representing parameterizations for one or more temperature ranges
 *  @return       Returns the pointer to the newly allocated
 *                SpeciesThermoInterpType object
 */
SpeciesThermoInterpType* newSpeciesThermoInterpType(const XML_Node& thermoNode);

}

#endif
