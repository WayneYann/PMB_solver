/**
 *  @file InterfaceKinetics.cpp
 */

// Copyright 2002  California Institute of Technology

#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/EdgeKinetics.h"
#include "cantera/kinetics/RateCoeffMgr.h"
#include "cantera/kinetics/ImplicitSurfChem.h"
#include "cantera/thermo/SurfPhase.h"

#include <cstdio>

using namespace std;

namespace Cantera
{

InterfaceKinetics::InterfaceKinetics(thermo_t* thermo) :
    m_redo_rates(false),
    m_surf(0),
    m_integrator(0),
    m_logp0(0.0),
    m_logc0(0.0),
    m_ROP_ok(false),
    m_temp(0.0),
    m_logtemp(0.0),
    m_finalized(false),
    m_has_coverage_dependence(false),
    m_has_electrochem_rxns(false),
    m_has_exchange_current_density_formulation(false),
    m_phaseExistsCheck(false),
    m_ioFlag(0)
{
    if (thermo != 0) {
        addPhase(*thermo);
    }
}

InterfaceKinetics::~InterfaceKinetics()
{
    delete m_integrator;
}

InterfaceKinetics::InterfaceKinetics(const InterfaceKinetics& right)
{
    // Call the assignment operator
    operator=(right);
}

InterfaceKinetics& InterfaceKinetics::operator=(const InterfaceKinetics& right)
{
    // Check for self assignment.
    if (this == &right) {
        return *this;
    }

    Kinetics::operator=(right);

    m_grt = right.m_grt;
    m_revindex = right.m_revindex;
    m_rates = right.m_rates;
    m_redo_rates = right.m_redo_rates;
    m_irrev = right.m_irrev;
    m_conc = right.m_conc;
    m_actConc = right.m_actConc;
    m_mu0 = right.m_mu0;
    m_mu = right.m_mu;
    m_mu0_Kc = right.m_mu0_Kc;
    m_phi = right.m_phi;
    m_pot = right.m_pot;
    deltaElectricEnergy_ = right.deltaElectricEnergy_;
    m_E = right.m_E;
    m_surf = right.m_surf; //DANGER - shallow copy
    m_integrator = right.m_integrator; //DANGER - shallow copy
    m_beta = right.m_beta;
    m_ctrxn = right.m_ctrxn;
    m_ctrxn_BVform = right.m_ctrxn_BVform;
    m_ctrxn_ecdf = right.m_ctrxn_ecdf;
    m_StandardConc = right.m_StandardConc;
    m_deltaG0 = right.m_deltaG0;
    m_deltaG = right.m_deltaG;
    m_ProdStanConcReac = right.m_ProdStanConcReac;
    m_logp0 = right.m_logp0;
    m_logc0 = right.m_logc0;
    m_ROP_ok = right.m_ROP_ok;
    m_temp = right.m_temp;
    m_logtemp = right.m_logtemp;
    m_finalized = right.m_finalized;
    m_has_coverage_dependence = right.m_has_coverage_dependence;
    m_has_electrochem_rxns = right.m_has_electrochem_rxns;
    m_has_exchange_current_density_formulation = right.m_has_exchange_current_density_formulation;
    m_phaseExistsCheck = right.m_phaseExistsCheck;
    m_phaseExists = right.m_phaseExists;
    m_phaseIsStable = right.m_phaseIsStable;
    m_rxnPhaseIsReactant = right.m_rxnPhaseIsReactant;
    m_rxnPhaseIsProduct = right.m_rxnPhaseIsProduct;
    m_ioFlag = right.m_ioFlag;

    return *this;
}

int InterfaceKinetics::type() const
{
    return cInterfaceKinetics;
}

Kinetics* InterfaceKinetics::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
{
    InterfaceKinetics* iK = new InterfaceKinetics(*this);
    iK->assignShallowPointers(tpVector);
    return iK;
}

void InterfaceKinetics::setElectricPotential(int n, doublereal V)
{
    thermo(n).setElectricPotential(V);
    m_redo_rates = true;
}

void InterfaceKinetics::_update_rates_T()
{
    // First task is update the electrical potentials from the Phases
    _update_rates_phi();
    if (m_has_coverage_dependence) {
        m_surf->getCoverages(m_actConc.data());
        m_rates.update_C(m_actConc.data());
        m_redo_rates = true;
    }

    // Go find the temperature from the surface
    doublereal T = thermo(surfacePhaseIndex()).temperature();
    m_redo_rates = true;
    if (T != m_temp || m_redo_rates) {
        m_logtemp = log(T);

        //  Calculate the forward rate constant by calling m_rates and store it in m_rfn[]
        m_rates.update(T, m_logtemp, m_rfn.data());
        applyStickingCorrection(m_rfn.data());

        // If we need to do conversions between exchange current density
        // formulation and regular formulation (either way) do it here.
        if (m_has_exchange_current_density_formulation) {
            convertExchangeCurrentDensityFormulation(m_rfn.data());
        }
        if (m_has_electrochem_rxns) {
            applyVoltageKfwdCorrection(m_rfn.data());
        }
        m_temp = T;
        updateKc();
        m_ROP_ok = false;
        m_redo_rates = false;
    }
}

void InterfaceKinetics::_update_rates_phi()
{
    // Store electric potentials for each phase in the array m_phi[].
    for (size_t n = 0; n < nPhases(); n++) {
        if (thermo(n).electricPotential() != m_phi[n]) {
            m_phi[n] = thermo(n).electricPotential();
            m_redo_rates = true;
        }
    }
}

void InterfaceKinetics::_update_rates_C()
{
    for (size_t n = 0; n < nPhases(); n++) {
        const ThermoPhase* tp = m_thermo[n];
        /*
         * We call the getActivityConcentrations function of each ThermoPhase
         * class that makes up this kinetics object to obtain the generalized
         * concentrations for species within that class. This is collected in
         * the vector m_conc. m_start[] are integer indices for that vector
         * denoting the start of the species for each phase.
         */
        tp->getActivityConcentrations(m_actConc.data() + m_start[n]);

        // Get regular concentrations too
        tp->getConcentrations(m_conc.data() + m_start[n]);
    }
    m_ROP_ok = false;
}

void InterfaceKinetics::getActivityConcentrations(doublereal* const conc)
{
    _update_rates_C();
    copy(m_actConc.begin(), m_actConc.end(), conc);
}

void InterfaceKinetics::updateKc()
{
    fill(m_rkcn.begin(), m_rkcn.end(), 0.0);

    if (m_revindex.size() > 0) {
        /*
         * Get the vector of standard state electrochemical potentials for
         * species in the Interfacial kinetics object and store it in m_mu0[]
         * and m_mu0_Kc[]
         */
        updateMu0();
        doublereal rrt = 1.0 / thermo(0).RT();

        // compute Delta mu^0 for all reversible reactions
        getRevReactionDelta(m_mu0_Kc.data(), m_rkcn.data());

        for (size_t i = 0; i < m_revindex.size(); i++) {
            size_t irxn = m_revindex[i];
            if (irxn == npos || irxn >= nReactions()) {
                throw CanteraError("InterfaceKinetics", "illegal value: irxn = {}", irxn);
            }
            // WARNING this may overflow HKM
            m_rkcn[irxn] = exp(m_rkcn[irxn]*rrt);
        }
        for (size_t i = 0; i != m_irrev.size(); ++i) {
            m_rkcn[ m_irrev[i] ] = 0.0;
        }
    }
}

void InterfaceKinetics::updateMu0()
{
    // First task is update the electrical potentials from the Phases
    _update_rates_phi();

    updateExchangeCurrentQuantities();
    size_t nsp, ik = 0;
    size_t np = nPhases();
    for (size_t n = 0; n < np; n++) {
        thermo(n).getStandardChemPotentials(m_mu0.data() + m_start[n]);
        nsp = thermo(n).nSpecies();
        for (size_t k = 0; k < nsp; k++) {
            m_mu0_Kc[ik] = m_mu0[ik] + Faraday * m_phi[n] * thermo(n).charge(k);
            m_mu0_Kc[ik] -= thermo(0).RT() * thermo(n).logStandardConc(k);
            ik++;
        }
    }
}

void InterfaceKinetics::checkPartialEquil()
{
    // First task is update the electrical potentials from the Phases
    _update_rates_phi();

    vector_fp dmu(nTotalSpecies(), 0.0);
    vector_fp rmu(std::max<size_t>(nReactions(), 1), 0.0);
    if (m_revindex.size() > 0) {
        cout << "T = " << thermo(0).temperature() << " " << thermo(0).RT() << endl;
        size_t nsp, ik=0;
        doublereal delta;
        for (size_t n = 0; n < nPhases(); n++) {
            thermo(n).getChemPotentials(dmu.data() + m_start[n]);
            nsp = thermo(n).nSpecies();
            for (size_t k = 0; k < nsp; k++) {
                delta = Faraday * m_phi[n] * thermo(n).charge(k);
                dmu[ik] += delta;
                ik++;
            }
        }

        // compute Delta mu^ for all reversible reactions
        getRevReactionDelta(dmu.data(), rmu.data());
        updateROP();
        for (size_t i = 0; i < m_revindex.size(); i++) {
            size_t irxn = m_revindex[i];
            writelog("Reaction {} {}\n",
                     reactionString(irxn), rmu[irxn]/thermo(0).RT());
            writelogf("%12.6e  %12.6e  %12.6e  %12.6e \n",
                      m_ropf[irxn], m_ropr[irxn], m_ropnet[irxn],
                      m_ropnet[irxn]/(m_ropf[irxn] + m_ropr[irxn]));
        }
    }
}

void InterfaceKinetics::getEquilibriumConstants(doublereal* kc)
{
    updateMu0();
    doublereal rrt = 1.0 / thermo(0).RT();
    std::fill(kc, kc + nReactions(), 0.0);
    getReactionDelta(m_mu0_Kc.data(), kc);
    for (size_t i = 0; i < nReactions(); i++) {
        kc[i] = exp(-kc[i]*rrt);
    }
}

void InterfaceKinetics::updateExchangeCurrentQuantities()
{
    // Calculate:
    //   - m_StandardConc[]
    //   - m_ProdStanConcReac[]
    //   - m_deltaG0[]
    //   - m_mu0[]

    // First collect vectors of the standard Gibbs free energies of the
    // species and the standard concentrations
    //   - m_mu0
    //   - m_StandardConc
    size_t ik = 0;

    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(m_mu0.data() + m_start[n]);
        size_t nsp = thermo(n).nSpecies();
        for (size_t k = 0; k < nsp; k++) {
            m_StandardConc[ik] = thermo(n).standardConcentration(k);
            ik++;
        }
    }

    getReactionDelta(m_mu0.data(), m_deltaG0.data());

    //  Calculate the product of the standard concentrations of the reactants
    for (size_t i = 0; i < nReactions(); i++) {
        m_ProdStanConcReac[i] = 1.0;
    }
    m_reactantStoich.multiply(m_StandardConc.data(), m_ProdStanConcReac.data());
}

void InterfaceKinetics::applyVoltageKfwdCorrection(doublereal* const kf)
{
    // Compute the electrical potential energy of each species
    size_t ik = 0;
    for (size_t n = 0; n < nPhases(); n++) {
        size_t nsp = thermo(n).nSpecies();
        for (size_t k = 0; k < nsp; k++) {
            m_pot[ik] = Faraday * thermo(n).charge(k) * m_phi[n];
            ik++;
        }
    }

    // Compute the change in electrical potential energy for each reaction. This
    // will only be non-zero if a potential difference is present.
    getReactionDelta(m_pot.data(), deltaElectricEnergy_.data());

    // Modify the reaction rates. Only modify those with a non-zero activation
    // energy. Below we decrease the activation energy below zero but in some
    // debug modes we print out a warning message about this.

    // NOTE, there is some discussion about this point. Should we decrease the
    // activation energy below zero? I don't think this has been decided in any
    // definitive way. The treatment below is numerically more stable, however.
    doublereal eamod;
    for (size_t i = 0; i < m_beta.size(); i++) {
        size_t irxn = m_ctrxn[i];

        // If we calculate the BV form directly, we don't add the voltage
        // correction to the forward reaction rate constants.
        if (m_ctrxn_BVform[i] == 0) {
            eamod = m_beta[i] * deltaElectricEnergy_[irxn];
            if (eamod != 0.0) {
                kf[irxn] *= exp(-eamod/thermo(0).RT());
            }
        }
    }
}

void InterfaceKinetics::convertExchangeCurrentDensityFormulation(doublereal* const kfwd)
{
    updateExchangeCurrentQuantities();
    // Loop over all reactions which are defined to have a voltage transfer
    // coefficient that affects the activity energy for the reaction
    for (size_t i = 0; i < m_ctrxn.size(); i++) {
        size_t irxn = m_ctrxn[i];

        // Determine whether the reaction rate constant is in an exchange
        // current density formulation format.
        int iECDFormulation = m_ctrxn_ecdf[i];
        if (iECDFormulation) {
            // If the BV form is to be converted into the normal form then we go
            // through this process. If it isn't to be converted, then we don't
            // go through this process.
            //
            // We need to have the straight chemical reaction rate constant to
            // come out of this calculation.
            if (m_ctrxn_BVform[i] == 0) {
                //  Calculate the term and modify the forward reaction
                double tmp = exp(- m_beta[i] * m_deltaG0[irxn] / thermo(0).RT());
                double tmp2 = m_ProdStanConcReac[irxn];
                tmp *= 1.0 / tmp2 / Faraday;
                kfwd[irxn] *= tmp;
            }
            //  If BVform is nonzero we don't need to do anything.
        } else {
            // kfwd[] is the chemical reaction rate constant
            //
            // If we are to calculate the BV form directly, then we will do the
            // reverse. We will calculate the exchange current density
            // formulation here and substitute it.
            if (m_ctrxn_BVform[i] != 0) {
                // Calculate the term and modify the forward reaction rate
                // constant so that it's in the exchange current density
                // formulation format
                double tmp = exp(m_beta[i] * m_deltaG0[irxn] * thermo(0).RT());
                double tmp2 = m_ProdStanConcReac[irxn];
                tmp *= Faraday * tmp2;
                kfwd[irxn] *= tmp;
            }
        }
    }
}

void InterfaceKinetics::getFwdRateConstants(doublereal* kfwd)
{
    updateROP();

    // copy rate coefficients into kfwd
    copy(m_rfn.begin(), m_rfn.end(), kfwd);

    // multiply by perturbation factor
    multiply_each(kfwd, kfwd + nReactions(), m_perturb.begin());
}

void InterfaceKinetics::getRevRateConstants(doublereal* krev, bool doIrreversible)
{
    getFwdRateConstants(krev);
    if (doIrreversible) {
        getEquilibriumConstants(m_ropnet.data());
        for (size_t i = 0; i < nReactions(); i++) {
            krev[i] /= m_ropnet[i];
        }
    } else {
        multiply_each(krev, krev + nReactions(), m_rkcn.begin());
    }
}

void InterfaceKinetics::updateROP()
{
    // evaluate rate constants and equilibrium constants at temperature and phi
    // (electric potential)
    _update_rates_T();
    // get updated activities (rates updated below)
    _update_rates_C();

    if (m_ROP_ok) {
        return;
    }

    // Copy the reaction rate coefficients, m_rfn, into m_ropf
    m_ropf = m_rfn;

    // Multiply by the perturbation factor
    multiply_each(m_ropf.begin(), m_ropf.end(), m_perturb.begin());

    // Copy the forward rate constants to the reverse rate constants
    m_ropr = m_ropf;

    // For reverse rates computed from thermochemistry, multiply
    // the forward rates copied into m_ropr by the reciprocals of
    // the equilibrium constants
    multiply_each(m_ropr.begin(), m_ropr.end(), m_rkcn.begin());

    // multiply ropf by the activity concentration reaction orders to obtain
    // the forward rates of progress.
    m_reactantStoich.multiply(m_actConc.data(), m_ropf.data());

    // For reversible reactions, multiply ropr by the activity concentration
    // products
    m_revProductStoich.multiply(m_actConc.data(), m_ropr.data());

    // Fix up these calculations for cases where the above formalism doesn't hold
    double OCV = 0.0;
    for (size_t jrxn = 0; jrxn != nReactions(); ++jrxn) {
        if (reactionType(jrxn) == BUTLERVOLMER_RXN) {
            // OK, the reaction rate constant contains the current density rate
            // constant calculation the rxnstoich calculation contained the
            // dependence of the current density on the activity concentrations
            // We finish up with the ROP calculation
            //
            // Calculate the overpotential of the reaction
            double nStoichElectrons=1;
            getDeltaGibbs(0);
            if (nStoichElectrons != 0.0) {
                OCV = m_deltaG[jrxn]/Faraday/ nStoichElectrons;
            }
        }
    }

    for (size_t j = 0; j != nReactions(); ++j) {
        m_ropnet[j] = m_ropf[j] - m_ropr[j];
    }

    // For reactions involving multiple phases, we must check that the phase
    // being consumed actually exists. This is particularly important for phases
    // that are stoichiometric phases containing one species with a unity
    // activity
    if (m_phaseExistsCheck) {
        for (size_t j = 0; j != nReactions(); ++j) {
            if ((m_ropr[j] > m_ropf[j]) && (m_ropr[j] > 0.0)) {
                for (size_t p = 0; p < nPhases(); p++) {
                    if (m_rxnPhaseIsProduct[j][p] && !m_phaseExists[p]) {
                        m_ropnet[j] = 0.0;
                        m_ropr[j] = m_ropf[j];
                        if (m_ropf[j] > 0.0) {
                            for (size_t rp = 0; rp < nPhases(); rp++) {
                                if (m_rxnPhaseIsReactant[j][rp] && !m_phaseExists[rp]) {
                                    m_ropnet[j] = 0.0;
                                    m_ropr[j] = m_ropf[j] = 0.0;
                                }
                            }
                        }
                    }
                    if (m_rxnPhaseIsReactant[j][p] && !m_phaseIsStable[p]) {
                        m_ropnet[j] = 0.0;
                        m_ropr[j] = m_ropf[j];
                    }
                }
            } else if ((m_ropf[j] > m_ropr[j]) && (m_ropf[j] > 0.0)) {
                for (size_t p = 0; p < nPhases(); p++) {
                    if (m_rxnPhaseIsReactant[j][p] && !m_phaseExists[p]) {
                        m_ropnet[j] = 0.0;
                        m_ropf[j] = m_ropr[j];
                        if (m_ropf[j] > 0.0) {
                            for (size_t rp = 0; rp < nPhases(); rp++) {
                                if (m_rxnPhaseIsProduct[j][rp] && !m_phaseExists[rp]) {
                                    m_ropnet[j] = 0.0;
                                    m_ropf[j] = m_ropr[j] = 0.0;
                                }
                            }
                        }
                    }
                    if (m_rxnPhaseIsProduct[j][p] && !m_phaseIsStable[p]) {
                        m_ropnet[j] = 0.0;
                        m_ropf[j] = m_ropr[j];
                    }
                }
            }
        }
    }
    m_ROP_ok = true;
}

void InterfaceKinetics::getDeltaGibbs(doublereal* deltaG)
{
    // Get the chemical potentials of the species in the all of the phases used
    // in the kinetics mechanism
    for (size_t n = 0; n < nPhases(); n++) {
        m_thermo[n]->getChemPotentials(m_mu.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_mu.data(), m_deltaG.data());
    if (deltaG != 0 && (m_deltaG.data() != deltaG)) {
        for (size_t j = 0; j < nReactions(); ++j) {
            deltaG[j] = m_deltaG[j];
        }
    }
}

void InterfaceKinetics::getDeltaElectrochemPotentials(doublereal* deltaM)
{
    // Get the chemical potentials of the species
    size_t np = nPhases();
    for (size_t n = 0; n < np; n++) {
        thermo(n).getElectrochemPotentials(m_grt.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_grt.data(), deltaM);
}

void InterfaceKinetics::getDeltaEnthalpy(doublereal* deltaH)
{
    // Get the partial molar enthalpy of all species
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEnthalpies(m_grt.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(m_grt.data(), deltaH);
}

void InterfaceKinetics::getDeltaEntropy(doublereal* deltaS)
{
    // Get the partial molar entropy of all species in all of the phases
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getPartialMolarEntropies(m_grt.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(m_grt.data(), deltaS);
}

void InterfaceKinetics::getDeltaSSGibbs(doublereal* deltaGSS)
{
    // Get the standard state chemical potentials of the species. This is the
    // array of chemical potentials at unit activity We define these here as the
    // chemical potentials of the pure species at the temperature and pressure
    // of the solution.
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getStandardChemPotentials(m_mu0.data() + m_start[n]);
    }

    // Use the stoichiometric manager to find deltaG for each reaction.
    getReactionDelta(m_mu0.data(), deltaGSS);
}

void InterfaceKinetics::getDeltaSSEnthalpy(doublereal* deltaH)
{
    // Get the standard state enthalpies of the species. This is the array of
    // chemical potentials at unit activity We define these here as the
    // enthalpies of the pure species at the temperature and pressure of the
    // solution.
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getEnthalpy_RT(m_grt.data() + m_start[n]);
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= thermo(0).RT();
    }

    // Use the stoichiometric manager to find deltaH for each reaction.
    getReactionDelta(m_grt.data(), deltaH);
}

void InterfaceKinetics::getDeltaSSEntropy(doublereal* deltaS)
{
    // Get the standard state entropy of the species. We define these here as
    // the entropies of the pure species at the temperature and pressure of the
    // solution.
    for (size_t n = 0; n < nPhases(); n++) {
        thermo(n).getEntropy_R(m_grt.data() + m_start[n]);
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_grt[k] *= GasConstant;
    }

    // Use the stoichiometric manager to find deltaS for each reaction.
    getReactionDelta(m_grt.data(), deltaS);
}

bool InterfaceKinetics::addReaction(shared_ptr<Reaction> r_base)
{
    size_t i = nReactions();
    bool added = Kinetics::addReaction(r_base);
    if (!added) {
        return false;
    }

    InterfaceReaction& r = dynamic_cast<InterfaceReaction&>(*r_base);
    SurfaceArrhenius rate = buildSurfaceArrhenius(i, r);
    m_rates.install(i, rate);

    // Turn on the global flag indicating surface coverage dependence
    if (!r.coverage_deps.empty()) {
        m_has_coverage_dependence = true;
    }

    // Store activation energy
    m_E.push_back(rate.activationEnergy_R());

    ElectrochemicalReaction* re = dynamic_cast<ElectrochemicalReaction*>(&r);
    if (re) {
        m_has_electrochem_rxns = true;
        m_beta.push_back(re->beta);
        m_ctrxn.push_back(i);
        if (re->exchange_current_density_formulation) {
            m_has_exchange_current_density_formulation = true;
            m_ctrxn_ecdf.push_back(1);
        } else {
            m_ctrxn_ecdf.push_back(0);
        }
        m_ctrxn_resistivity_.push_back(re->film_resistivity);

        if (r.reaction_type == BUTLERVOLMER_NOACTIVITYCOEFFS_RXN ||
            r.reaction_type == BUTLERVOLMER_RXN ||
            r.reaction_type == SURFACEAFFINITY_RXN ||
            r.reaction_type == GLOBAL_RXN) {
            //   Specify alternative forms of the electrochemical reaction
            if (r.reaction_type == BUTLERVOLMER_RXN) {
                m_ctrxn_BVform.push_back(1);
            } else if (r.reaction_type == BUTLERVOLMER_NOACTIVITYCOEFFS_RXN) {
                m_ctrxn_BVform.push_back(2);
            } else {
                // set the default to be the normal forward / reverse calculation method
                m_ctrxn_BVform.push_back(0);
            }
            if (!r.orders.empty()) {
                vector_fp orders(nTotalSpecies(), 0.0);
                for (const auto& order : r.orders) {
                    orders[kineticsSpeciesIndex(order.first)] = order.second;
                }
            }
        } else {
            m_ctrxn_BVform.push_back(0);
            if (re->film_resistivity > 0.0) {
                throw CanteraError("InterfaceKinetics::addReaction()",
                                   "film resistivity set for elementary reaction");
            }
        }
    }

    if (r.reversible) {
        m_revindex.push_back(i);
    } else {
        m_irrev.push_back(i);
    }

    m_rxnPhaseIsReactant.emplace_back(nPhases(), false);
    m_rxnPhaseIsProduct.emplace_back(nPhases(), false);

    for (const auto& sp : r.reactants) {
        size_t k = kineticsSpeciesIndex(sp.first);
        size_t p = speciesPhaseIndex(k);
        m_rxnPhaseIsReactant[i][p] = true;
    }
    for (const auto& sp : r.products) {
        size_t k = kineticsSpeciesIndex(sp.first);
        size_t p = speciesPhaseIndex(k);
        m_rxnPhaseIsProduct[i][p] = true;
    }
    return true;
}

void InterfaceKinetics::modifyReaction(size_t i, shared_ptr<Reaction> r_base)
{
    Kinetics::modifyReaction(i, r_base);
    InterfaceReaction& r = dynamic_cast<InterfaceReaction&>(*r_base);
    SurfaceArrhenius rate = buildSurfaceArrhenius(npos, r);
    m_rates.replace(i, rate);

    // Invalidate cached data
    m_redo_rates = true;
    m_temp += 0.1;
}

SurfaceArrhenius InterfaceKinetics::buildSurfaceArrhenius(
    size_t i, InterfaceReaction& r)
{
    double A_rate = r.rate.preExponentialFactor();
    double b_rate = r.rate.temperatureExponent();

    if (r.is_sticking_coefficient) {
        // Identify the interface phase
        size_t iInterface = npos;
        size_t min_dim = 4;
        for (size_t n = 0; n < nPhases(); n++) {
            if (thermo(n).nDim() < min_dim) {
                iInterface = n;
                min_dim = thermo(n).nDim();
            }
        }

        b_rate += 0.5;
        std::string sticking_species = r.sticking_species;
        if (sticking_species == "") {
            // Identify the sticking species if not explicitly given
            bool foundStick = false;
            for (const auto& sp : r.reactants) {
                size_t iPhase = speciesPhaseIndex(kineticsSpeciesIndex(sp.first));
                if (iPhase != iInterface) {
                    // Non-interface species. There should be exactly one of these
                    if (foundStick) {
                        throw CanteraError("InterfaceKinetics::addReaction",
                            "Multiple non-interface species found"
                            "in sticking reaction: '" + r.equation() + "'");
                    }
                    foundStick = true;
                    sticking_species = sp.first;
                }
            }
            if (!foundStick) {
                throw CanteraError("InterfaceKinetics::addReaction",
                    "No non-interface species found"
                    "in sticking reaction: '" + r.equation() + "'");
            }
        }

        double surface_order = 0.0;
        // Adjust the A-factor
        for (const auto& sp : r.reactants) {
            size_t iPhase = speciesPhaseIndex(kineticsSpeciesIndex(sp.first));
            const ThermoPhase& p = thermo(iPhase);
            const ThermoPhase& surf = thermo(surfacePhaseIndex());
            size_t k = p.speciesIndex(sp.first);
            if (sp.first == sticking_species) {
                A_rate *= sqrt(GasConstant/(2*Pi*p.molecularWeight(k)));
            } else {
                // Non-sticking species. Convert from coverages used in the
                // sticking probability expression to the concentration units
                // used in the mass action rate expression. For surface phases,
                // the dependence on the site density is incorporated when the
                // rate constant is evaluated, since we don't assume that the
                // site density is known at this time.
                double order = getValue(r.orders, sp.first, sp.second);
                if (&p == &surf) {
                    A_rate *= pow(p.size(k), order);
                    surface_order += order;
                } else {
                    A_rate *= pow(p.standardConcentration(k), -order);
                }
            }
        }
        if (i != npos) {
            m_sticking_orders.emplace_back(i, surface_order);
        }
    }

    SurfaceArrhenius rate(A_rate, b_rate, r.rate.activationEnergy_R());

    // Set up coverage dependencies
    for (const auto& sp : r.coverage_deps) {
        size_t k = thermo(reactionPhaseIndex()).speciesIndex(sp.first);
        rate.addCoverageDependence(k, sp.second.a, sp.second.m, sp.second.E);
    }
    return rate;
}

void InterfaceKinetics::setIOFlag(int ioFlag)
{
    m_ioFlag = ioFlag;
    if (m_integrator) {
        m_integrator->setIOFlag(ioFlag);
    }
}

void InterfaceKinetics::addPhase(thermo_t& thermo)
{
    Kinetics::addPhase(thermo);
    m_phaseExists.push_back(true);
    m_phaseIsStable.push_back(true);
}

void InterfaceKinetics::init()
{
    m_kk = 0;
    for (size_t n = 0; n < nPhases(); n++) {
        m_kk += thermo(n).nSpecies();
    }
    m_actConc.resize(m_kk);
    m_conc.resize(m_kk);
    m_mu0.resize(m_kk);
    m_mu.resize(m_kk);
    m_mu0_Kc.resize(m_kk);
    m_grt.resize(m_kk);
    m_pot.resize(m_kk, 0.0);
    m_phi.resize(nPhases(), 0.0);
}

void InterfaceKinetics::finalize()
{
    Kinetics::finalize();
    deltaElectricEnergy_.resize(nReactions());
    size_t ks = reactionPhaseIndex();
    if (ks == npos) throw CanteraError("InterfaceKinetics::finalize",
                                           "no surface phase is present.");

    // Check to see that the interface routine has a dimension of 2
    m_surf = (SurfPhase*)&thermo(ks);
    if (m_surf->nDim() != 2) {
        throw CanteraError("InterfaceKinetics::finalize",
                           "expected interface dimension = 2, but got dimension = {}",
                           m_surf->nDim());
    }
    m_StandardConc.resize(m_kk, 0.0);
    m_deltaG0.resize(nReactions(), 0.0);
    m_deltaG.resize(nReactions(), 0.0);
    m_ProdStanConcReac.resize(nReactions(), 0.0);

    if (m_thermo.size() != m_phaseExists.size()) {
        throw CanteraError("InterfaceKinetics::finalize", "internal error");
    }
    m_finalized = true;
}

doublereal InterfaceKinetics::electrochem_beta(size_t irxn) const
{
    for (size_t i = 0; i < m_ctrxn.size(); i++) {
        if (m_ctrxn[i] == irxn) {
            return m_beta[i];
        }
    }
    return 0.0;
}

bool InterfaceKinetics::ready() const
{
    return m_finalized;
}

void InterfaceKinetics::advanceCoverages(doublereal tstep)
{
    if (m_integrator == 0) {
        vector<InterfaceKinetics*> k;
        k.push_back(this);
        m_integrator = new ImplicitSurfChem(k);
        m_integrator->initialize();
    }
    m_integrator->integrate(0.0, tstep);
    delete m_integrator;
    m_integrator = 0;
}

void InterfaceKinetics::solvePseudoSteadyStateProblem(
    int ifuncOverride, doublereal timeScaleOverride)
{
    // create our own solver object
    if (m_integrator == 0) {
        vector<InterfaceKinetics*> k;
        k.push_back(this);
        m_integrator = new ImplicitSurfChem(k);
        m_integrator->initialize();
    }
    m_integrator->setIOFlag(m_ioFlag);
    // New direct method to go here
    m_integrator->solvePseudoSteadyStateProblem(ifuncOverride, timeScaleOverride);
}

void InterfaceKinetics::setPhaseExistence(const size_t iphase, const int exists)
{
    if (iphase >= m_thermo.size()) {
        throw CanteraError("InterfaceKinetics:setPhaseExistence", "out of bounds");
    }
    if (exists) {
        if (!m_phaseExists[iphase]) {
            m_phaseExistsCheck--;
            m_phaseExistsCheck = std::max(m_phaseExistsCheck, 0);
            m_phaseExists[iphase] = true;
        }
        m_phaseIsStable[iphase] = true;
    } else {
        if (m_phaseExists[iphase]) {
            m_phaseExistsCheck++;
            m_phaseExists[iphase] = false;
        }
        m_phaseIsStable[iphase] = false;
    }
}

int InterfaceKinetics::phaseExistence(const size_t iphase) const
{
    if (iphase >= m_thermo.size()) {
        throw CanteraError("InterfaceKinetics:phaseExistence()", "out of bounds");
    }
    return m_phaseExists[iphase];
}

int InterfaceKinetics::phaseStability(const size_t iphase) const
{
    if (iphase >= m_thermo.size()) {
        throw CanteraError("InterfaceKinetics:phaseStability()", "out of bounds");
    }
    return m_phaseIsStable[iphase];
}

void InterfaceKinetics::setPhaseStability(const size_t iphase, const int isStable)
{
    if (iphase >= m_thermo.size()) {
        throw CanteraError("InterfaceKinetics:setPhaseStability", "out of bounds");
    }
    if (isStable) {
        m_phaseIsStable[iphase] = true;
    } else {
        m_phaseIsStable[iphase] = false;
    }
}

void InterfaceKinetics::determineFwdOrdersBV(ElectrochemicalReaction& r, vector_fp& fwdFullOrders)
{
    // Start out with the full ROP orders vector.
    // This vector will have the BV exchange current density orders in it.
    fwdFullOrders.assign(nTotalSpecies(), 0.0);
    for (const auto& order : r.orders) {
        fwdFullOrders[kineticsSpeciesIndex(order.first)] = order.second;
    }

    //   forward and reverse beta values
    double betaf = r.beta;

    // Loop over the reactants doing away with the BV terms.
    // This should leave the reactant terms only, even if they are non-mass action.
    for (const auto& sp : r.reactants) {
        size_t k = kineticsSpeciesIndex(sp.first);
        fwdFullOrders[k] += betaf * sp.second;
        // just to make sure roundoff doesn't leave a term that should be zero (haven't checked this out yet)
        if (abs(fwdFullOrders[k]) < 0.00001) {
            fwdFullOrders[k] = 0.0;
        }
    }

    // Loop over the products doing away with the BV terms.
    // This should leave the reactant terms only, even if they are non-mass action.
    for (const auto& sp : r.products) {
        size_t k = kineticsSpeciesIndex(sp.first);
        fwdFullOrders[k] -= betaf * sp.second;
        // just to make sure roundoff doesn't leave a term that should be zero (haven't checked this out yet)
        if (abs(fwdFullOrders[k]) < 0.00001) {
            fwdFullOrders[k] = 0.0;
        }
    }
}

void InterfaceKinetics::applyStickingCorrection(double* kf)
{
    if (m_sticking_orders.empty()) {
        return;
    }

    static const int cacheId = m_cache.getId();
    CachedArray cached = m_cache.getArray(cacheId);
    vector_fp& factors = cached.value;

    SurfPhase& surf = dynamic_cast<SurfPhase&>(thermo(reactionPhaseIndex()));
    double n0 = surf.siteDensity();
    if (!cached.validate(n0)) {
        factors.resize(m_sticking_orders.size());
        for (size_t n = 0; n < m_sticking_orders.size(); n++) {
            factors[n] = pow(n0, -m_sticking_orders[n].second);
        }
    }

    for (size_t n = 0; n < m_sticking_orders.size(); n++) {
        kf[m_sticking_orders[n].first] *= factors[n];
    }
}


void EdgeKinetics::finalize()
{
    // Note we can't call the Interface::finalize() routine because we need to
    // check for a dimension of 1 below. Therefore, we have to malloc room in
    // arrays that would normally be handled by the
    // InterfaceKinetics::finalize() call.
    Kinetics::finalize();

    size_t safe_reaction_size = std::max<size_t>(nReactions(), 1);
    deltaElectricEnergy_.resize(safe_reaction_size);
    size_t ks = reactionPhaseIndex();
    if (ks == npos) throw CanteraError("EdgeKinetics::finalize",
                                           "no surface phase is present.");

    // Check to see edge phase has a dimension of 1
    m_surf = (SurfPhase*)&thermo(ks);
    if (m_surf->nDim() != 1) {
        throw CanteraError("EdgeKinetics::finalize",
                           "expected interface dimension = 1, but got dimension = {}",
                           m_surf->nDim());
    }
    m_StandardConc.resize(m_kk, 0.0);
    m_deltaG0.resize(safe_reaction_size, 0.0);
    m_deltaG.resize(safe_reaction_size, 0.0);

    m_ProdStanConcReac.resize(safe_reaction_size, 0.0);

    if (m_thermo.size() != m_phaseExists.size()) {
        throw CanteraError("InterfaceKinetics::finalize", "internal error");
    }

    // Guarantee that these arrays can be converted to double* even in the
    // special case where there are no reactions defined.
    if (!nReactions()) {
        m_perturb.resize(1, 1.0);
        m_ropf.resize(1, 0.0);
        m_ropr.resize(1, 0.0);
        m_ropnet.resize(1, 0.0);
        m_rkcn.resize(1, 0.0);
    }
    m_finalized = true;
}

}
