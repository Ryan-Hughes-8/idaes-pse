#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

from pyomo.environ import units as units
from idaes.core.util.constants import Constants as const
from pyomo.environ import (
    # Constraint,
    # Var,
    Param,
    # value,
    # Set,
    exp,
    log,
    # PositiveReals,
    # NonPositiveReals,
    # TransformationFactory,
    # units,
    # Block,
)


def add_Langmuir_Freundlich_parameters(blk):
    """
    Method for adding parameters of the Langmuir isotherm model.
    """

    blk.LF_qm = Param(
        blk.isotherm_components,
        initialize={"CO2": 5.0, "N2": 0.0},
        units=units.mol / units.kg,
        doc="Langmuir-Freundlich saturation capacity",
    )
    blk.LF_b0 = Param(
        blk.isotherm_components,
        initialize={"CO2": 1e-10, "N2": 0.0},
        units=units.Pa**-1,
        doc="Langmuir-Freundlich b0 (pre-exponential)",
    )
    blk.LF_E = Param(
        blk.isotherm_components,
        initialize={"CO2": -35.0, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="Langmuir-Freundlich E",
    )
    blk.LF_nu = Param(
        blk.isotherm_components,
        initialize={"CO2": 0.8, "N2": 0.0},
        units=units.dimensionless,
        doc="Langmuir-Freundlich nu",
    )


def Langmuir_Freundlich_isotherm(blk, i, pressure, temperature):
    """
    Method to add isotherm for components.
    Isotherm equation: Langmuir-Freundlich

    NOTE: CO2 is considered as the only adsorbing component

    Keyword Arguments:
        i : component
        pressure : partial pressure of components
        temperature : temperature

    """

    T = temperature
    p = {}
    loading = {}

    for j in blk.isotherm_components:
        p[j] = units.convert(pressure[j], to_units=units.Pa)

    if i == "CO2":

        affinity = blk.LF_b0[i] * exp(
            units.convert(-blk.LF_E[i], to_units=units.J / units.mol)
            / const.gas_constant
            / T
        )

        loading[i] = (
            blk.LF_qm[i]
            * affinity
            * p[i] ** blk.LF_nu[i]
            / (1 + affinity * p[i] ** blk.LF_nu[i])
        )

    elif i == "N2":
        # no adsorption is assumed of N2 in this adsorbent
        loading[i] = 1e-10 * units.mol / units.kg

    return loading[i]


def add_Henry_parameters(blk):
    """
    Method for adding parameters of the Langmuir isotherm model.
    """

    blk.Henry_b0 = Param(
        blk.isotherm_components,
        initialize={"CO2": 1e-10, "N2": 0.0},
        units=units.Pa**-1,
        doc="Henry b0 (pre-exponential)",
    )
    blk.Henry_E = Param(
        blk.isotherm_components,
        initialize={"CO2": -35.0, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="Henry E",
    )


def Henry_isotherm(blk, i, pressure, temperature):
    """
    Method to add isotherm for components.
    Isotherm equation: Henry

    NOTE: CO2 is considered as the only adsorbing component

    Keyword Arguments:
        i : component
        pressure : partial pressure of components
        temperature : temperature

    """

    T = temperature
    p = {}
    loading = {}

    for j in blk.isotherm_components:
        p[j] = units.convert(pressure[j], to_units=units.Pa)

    if i == "CO2":

        affinity = blk.Henry_b0[i] * exp(
            units.convert(-blk.Henry_E[i], to_units=units.J / units.mol)
            / const.gas_constant
            / T
        )

        loading[i] = affinity * p[i]

    elif i == "N2":
        # no adsorption is assumed of N2 in this adsorbent
        loading[i] = 1e-10 * units.mol / units.kg

    return loading[i]


def add_Langmuir_parameters(blk):
    """
    Method for adding parameters of the Langmuir isotherm model.
    """

    blk.Langmuir_qm = Param(
        blk.isotherm_components,
        initialize={"CO2": 3.0, "N2": 0.0},
        units=units.mol / units.kg,
        doc="Langmuir saturation capacity",
    )
    blk.Langmuir_b0 = Param(
        blk.isotherm_components,
        initialize={"CO2": 1e-10, "N2": 0.0},
        units=units.Pa**-1,
        doc="Langmuir b0 (pre-exponential)",
    )
    blk.Langmuir_E = Param(
        blk.isotherm_components,
        initialize={"CO2": -35.0, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="Langmuir E",
    )


def Langmuir_isotherm(blk, i, pressure, temperature):
    """
    Method to add isotherm for components.
    Isotherm equation: Langmuir

    NOTE: CO2 is considered as the only adsorbing component

    Keyword Arguments:
        i : component
        pressure : partial pressure of components
        temperature : temperature

    """

    T = temperature
    p = {}
    loading = {}

    for j in blk.isotherm_components:
        p[j] = units.convert(pressure[j], to_units=units.Pa)

    if i == "CO2":

        affinity = blk.Langmuir_b0[i] * exp(
            units.convert(-blk.Langmuir_E[i], to_units=units.J / units.mol)
            / const.gas_constant
            / T
        )

        loading[i] = blk.Langmuir_qm[i] * affinity * p[i] / (1 + affinity * p[i])

    elif i == "N2":
        # no adsorption is assumed of N2 in this adsorbent
        loading[i] = 1e-10 * units.mol / units.kg

    return loading[i]


def add_dual_site_Langmuir_parameters(blk):
    """
    Method for adding parameters of the dual site Langmuir isotherm model.
    """

    blk.temperature_ref = Param(
        initialize=298.15,
        units=units.K,
        doc="Reference temperature",
    )
    blk.saturation_capacity_site_b = Param(
        blk.isotherm_components,
        initialize={"CO2": 2.387, "N2": 0.0},
        units=units.mol / units.kg,
        doc="saturation capacity at site b",
    )
    blk.saturation_capacity_site_d = Param(
        blk.isotherm_components,
        initialize={"CO2": 3.2711, "N2": 0.0},
        units=units.mol / units.kg,
        doc="saturation capacity at site d",
    )
    blk.dual_site_langmuir_constant_pre_exp_b = Param(
        blk.isotherm_components,
        initialize={"CO2": 5.519e-7, "N2": 0.0},
        units=units.meter**3 / units.mol,
        doc="dual site langmuir constant for site b",
    )
    blk.dual_site_langmuir_constant_pre_exp_d = Param(
        blk.isotherm_components,
        initialize={"CO2": 5.187e-08, "N2": 0.0},
        units=units.meter**3 / units.mol,
        doc="dual site langmuir constant for site d",
    )
    blk.internal_energy_b = Param(
        blk.isotherm_components,
        initialize={"CO2": -35.06, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="internal energy of site b",
    )
    blk.internal_energy_d = Param(
        blk.isotherm_components,
        initialize={"CO2": -28.95, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="internal energy of site d",
    )


def dual_site_Langmuir_isotherm(blk, i, pressure, temperature):
    """
    Method to add isotherm for components.
    Isotherm equation: Dual site Langmuir (DSL) isotherm

    NOTE: CO2 is considered as the only adsorbing component

    Keyword Arguments:
        i : component
        pressure : partial pressure of components
        temperature : temperature

    """
    T = temperature
    p = {}
    c = {}
    loading = {}

    for j in blk.isotherm_components:
        p[j] = units.convert(pressure[j], to_units=units.bar)
        c[j] = units.convert(
            (p[j] / const.gas_constant / T), to_units=units.mol / units.meter**3
        )

    if i == "CO2":

        dual_site_langmuir_constant_b = blk.dual_site_langmuir_constant_pre_exp_b[
            i
        ] * exp(
            units.convert(
                -blk.internal_energy_b[i],
                to_units=units.J / units.mol,
            )
            / const.gas_constant
            / T
        )

        dual_site_langmuir_constant_d = blk.dual_site_langmuir_constant_pre_exp_d[
            i
        ] * exp(
            units.convert(
                -blk.internal_energy_d[i],
                to_units=units.J / units.mol,
            )
            / const.gas_constant
            / T
        )

        # portion of the isotherm for site b
        loading_b = (
            blk.saturation_capacity_site_b[i]
            * dual_site_langmuir_constant_b
            * c[i]
            / (1 + dual_site_langmuir_constant_b * c[i])
        )

        # portion of the isotherm for site d
        loading_d = (
            blk.saturation_capacity_site_d[i]
            * dual_site_langmuir_constant_d
            * c[i]
            / (1 + dual_site_langmuir_constant_d * c[i])
        )

        loading[i] = loading_b + loading_d  # [mol/kg]

    elif i == "N2":
        # no adsorption is assumed of N2
        loading[i] = 1e-10 * units.mol / units.kg

    return loading[i]


def add_extended_Sips_parameters(blk):
    """
    Method to add extended sips isotherm parameters. CO2 and N2 are adsorbed.

    Reference: Hefti, M.; Marx, D.; Joss, L.; Mazzotti, M. Adsorption
    Equilibrium of Binary Mixtures of Carbon Dioxide and Nitrogen on
    Zeolites ZSM-5 and 13X. Microporous Mesoporous Materials, 215, 2014.

    """

    blk.temperature_ref = Param(
        initialize=298.15,
        units=units.K,
        doc="Reference temperature",
    )
    blk.saturation_capacity_ref = Param(
        blk.isotherm_components,
        initialize={"CO2": 7.268, "N2": 4.051},
        units=units.mol / units.kg,
        doc="Saturation capacity at reference temperature",
    )
    blk.saturation_capacity_exponential_factor = Param(
        blk.isotherm_components,
        initialize={"CO2": -0.61684, "N2": 0.0},
        units=units.dimensionless,
        doc="Isotherm fitting parameter",
    )
    blk.affinity_parameter_preexponential_factor = Param(
        blk.isotherm_components,
        initialize={"CO2": 1.129e-4, "N2": 5.8470e-5},
        units=units.bar**-1,
        doc="Affinity parameter pre-exponential factor",
    )
    blk.affinity_parameter_characteristic_energy = Param(
        blk.isotherm_components,
        initialize={"CO2": 28.389, "N2": 18.4740},
        units=units.kJ / units.mol,
        doc="Characteristic energy for the affinity parameter",
    )
    blk.heterogeneity_parameter_ref = Param(
        blk.isotherm_components,
        initialize={"CO2": 0.42456, "N2": 0.98624},
        units=units.dimensionless,
        doc="Heterogeneity parameter at reference temperature",
    )
    blk.heterogeneity_parameter_alpha = Param(
        blk.isotherm_components,
        initialize={"CO2": 0.72378, "N2": 0.0},
        units=units.dimensionless,
        doc="Isotherm fitting parameter",
    )


def extended_Sips_isotherm(blk, i, pressure, temperature):
    """
    Method to add isotherm for components.
    Isotherm equation: Extended Sips

    Keyword Arguments:
        i : component
        pressure : partial pressure of components
        temperature : temperature

    """
    T = temperature

    saturation_capacity = {}
    affinity_parameter = {}
    heterogeneity_parameter = {}
    p = {}

    for j in blk.isotherm_components:

        p[j] = units.convert(pressure[j], to_units=units.bar)

        saturation_capacity[j] = blk.saturation_capacity_ref[j] * exp(
            blk.saturation_capacity_exponential_factor[j]
            * (T / blk.temperature_ref - 1)
        )
        affinity_parameter[j] = blk.affinity_parameter_preexponential_factor[j] * exp(
            units.convert(
                blk.affinity_parameter_characteristic_energy[j],
                to_units=units.J / units.mol,
            )
            / const.gas_constant
            / T
        )
        heterogeneity_parameter[j] = blk.heterogeneity_parameter_ref[
            j
        ] + blk.heterogeneity_parameter_alpha[j] * (T / blk.temperature_ref - 1)

    loading = {}

    for j in blk.isotherm_components:

        loading[j] = (
            saturation_capacity[j]
            * (affinity_parameter[j] * p[j]) ** heterogeneity_parameter[j]
            / (
                1
                + sum(
                    (affinity_parameter[k] * p[k]) ** heterogeneity_parameter[k]
                    for k in blk.isotherm_components
                )
            )
        )

    return loading[i]


def add_weighted_DSL_parameters(blk):
    """
    Method to add isotherm parameters for the weighed dual-site
    Langmuir isotherm model. Default values taken for mmen-Mg-MOF-74.

    Reference: Joss, L.; Hefti, M.; Bjelobrk, Z.; Mazzotti, M.
    Investigating the potential of phase-change materials for CO2
    capture. Faraday Disc, 192, 2016

    """

    blk.temperature_ref = Param(
        initialize=313.15,
        units=units.K,
        doc="Reference temperature",
    )
    blk.lower_saturtion_capacity = Param(
        blk.isotherm_components,
        initialize={"CO2": 0.146, "N2": 0.0},
        units=units.mol / units.kg,
        doc="Lower isotherm saturation capacity",
    )
    blk.upper_saturtion_capacity = Param(
        blk.isotherm_components,
        initialize={"CO2": 3.478, "N2": 0.0},
        units=units.mol / units.kg,
        doc="Upper isotherm saturation capacity",
    )
    blk.lower_affinity_preexponential_factor = Param(
        blk.isotherm_components,
        initialize={"CO2": 0.009, "N2": 0.0},
        units=units.bar**-1,
        doc="Pre-exponential factor for the lower isotherm affinity parameter",
    )
    blk.upper_affinity_preexponential_factor_1 = Param(
        blk.isotherm_components,
        initialize={"CO2": 9.00e-07, "N2": 0.0},
        units=units.bar**-1,
        doc="Pre-exponential factor for the upper isotherm site 1 affinity parameter",
    )
    blk.upper_affinity_preexponential_factor_2 = Param(
        blk.isotherm_components,
        initialize={"CO2": 5.00e-04, "N2": 0.0},
        units=units.mol / units.kg / units.bar,
        doc="Pre-exponential factor for the upper isotherm site 2 affinity parameter",
    )
    blk.lower_affinity_characteristic_energy = Param(
        blk.isotherm_components,
        initialize={"CO2": 31.0, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="Characteristic energy for the lower isotherm affinity parameter",
    )
    blk.upper_affinity_characteristic_energy_1 = Param(
        blk.isotherm_components,
        initialize={"CO2": 59.0, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="Characteristic energy for the upper isotherm site 1 affinity parameter",
    )
    blk.upper_affinity_characteristic_energy_2 = Param(
        blk.isotherm_components,
        initialize={"CO2": 18.0, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="Characteristic energy for the upper isotherm site 2 affinity parameter",
    )
    blk.step_width_preexponential_factor = Param(
        blk.isotherm_components,
        initialize={"CO2": 1.24e-01, "N2": 0.0},
        units=units.dimensionless,
        doc="Pre-exponential factor for step width",
    )
    blk.step_width_exponential_factor = Param(
        blk.isotherm_components,
        initialize={"CO2": 0.0, "N2": 0.0},
        units=units.dimensionless,
        doc="Exponential factor for step width",
    )
    blk.weighting_function_exponent = Param(
        blk.isotherm_components,
        initialize={"CO2": 4.00, "N2": 0.0},
        units=units.dimensionless,
        doc="Isotherm weighting function exponent",
    )
    blk.step_partial_pressure_ref = Param(
        blk.isotherm_components,
        initialize={"CO2": 0.5 * 1e-3, "N2": 0.0},
        units=units.bar,
        doc="Step partial pressure at reference temperature",
    )
    blk.step_enthalpy = Param(
        blk.isotherm_components,
        initialize={"CO2": -74.1, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="Enthalpy of phase transition",
    )


def weighted_DSL_isotherm(blk, i, pressure, temperature):
    """
    Method to add isotherm for components.
    Isotherm equation: Weighted dual site Langmuir (w-DSL) isotherm
    Adsorbent: mmen-Mg(dobpdc): mmen-Mg-MOF-74
                mmen = N,N"-dimethylethylenediamine
                dobpdc4^- = 4,4"-dioxido-3,3"-biphenyldicarboxylate

    NOTE: the affinity of the material towards CO2 is not reduced by
          neither N2 nor H2O. N2 adsorption is negligible in
          mmen-Mg(dobpdc). Therefore, CO2 is considered as the only
          adsorbing component

    Keyword Arguments:
        i : component
        pressure : partial pressure of components
        temperature : temperature

    """
    T = temperature
    p = {}
    loading = {}

    for j in blk.isotherm_components:
        p[j] = units.convert(pressure[j], to_units=units.bar)

    if i == "CO2":

        lower_affinity_parameter = blk.lower_affinity_preexponential_factor[i] * exp(
            units.convert(
                blk.lower_affinity_characteristic_energy[i],
                to_units=units.J / units.mol,
            )
            / const.gas_constant
            / T
        )
        upper_affinity_parameter_1 = blk.upper_affinity_preexponential_factor_1[
            i
        ] * exp(
            units.convert(
                blk.upper_affinity_characteristic_energy_1[i],
                to_units=units.J / units.mol,
            )
            / const.gas_constant
            / T
        )
        upper_affinity_parameter_2 = blk.upper_affinity_preexponential_factor_2[
            i
        ] * exp(
            units.convert(
                blk.upper_affinity_characteristic_energy_2[i],
                to_units=units.J / units.mol,
            )
            / const.gas_constant
            / T
        )

        # lower portion of the isotherm
        lower_loading = (
            blk.lower_saturtion_capacity[i]
            * lower_affinity_parameter
            * p[i]
            / (1 + lower_affinity_parameter * p[i])
        )

        # upper portion of the isotherm
        upper_loading = (
            blk.upper_saturtion_capacity[i]
            * upper_affinity_parameter_1
            * p[i]
            / (1 + upper_affinity_parameter_1 * p[i])
        ) + upper_affinity_parameter_2 * p[i]

        #  weighting function
        sigma = blk.step_width_preexponential_factor[i] * exp(
            blk.step_width_exponential_factor[i]
            * (1 / blk.temperature_ref - 1 / T)
            * units.K
        )

        pstep = (
            blk.step_partial_pressure_ref[i]
            / units.bar
            * exp(
                (-blk.step_enthalpy[i] / const.gas_constant * 1e3 * units.J / units.kJ)
                * (1 / blk.temperature_ref - 1 / T)
            )
        )

        w = (
            exp((log(p[i] / units.bar) - log(pstep)) / sigma)
            / (1 + exp((log(p[i] / units.bar) - log(pstep)) / sigma))
        ) ** blk.weighting_function_exponent[i]

        loading[i] = lower_loading * (1 - w) + upper_loading * w  # [mol/kg]

    elif i == "N2":
        # no adsorption is assumed of N2 in mmen-Mg(dobpdc)
        loading[i] = 1e-10 * units.mol / units.kg

    return loading[i]


def add_Toth_parameters(blk):
    """
    Method to add adsorbent related parameters to run fixed bed TSA model.
    This method is to add parameters for polystyrene functionalized
    with primary amine.

    Elfvinga, J.; Bajamundia, C.; Kauppinena, J.; Sainiob, T. Modelling
    of equilibrium working capacity of PSA, TSA and TVSA processes for
    CO2 adsorption under direct air capture conditions. Journal of CO2
    Utilization, 22, 2017.

    """

    blk.temperature_ref = Param(
        initialize=298.15,
        units=units.K,
        doc="Reference temperature",
    )
    blk.saturation_capacity_ref = Param(
        blk.isotherm_components,
        initialize={"CO2": 1.71, "N2": 0.0},
        units=units.mol / units.kg,
        doc="Saturation capacity at reference temperature",
    )
    blk.affinity_preexponential_factor = Param(
        blk.isotherm_components,
        initialize={"CO2": 1.13e5, "N2": 0.0},
        units=units.bar**-1,
        doc="Pre-exponential factor for the affinity parameter",
    )
    blk.toth_constant_ref = Param(
        blk.isotherm_components,
        initialize={"CO2": 0.265, "N2": 0.0},
        units=units.dimensionless,
        doc="Toth constant at reference temperature",
    )
    blk.toth_constant_alpha = Param(
        blk.isotherm_components,
        initialize={"CO2": 0.601, "N2": 0.0},
        units=units.dimensionless,
        doc="Isotherm parameter",
    )
    blk.saturation_capacity_exponential_factor = Param(
        blk.isotherm_components,
        initialize={"CO2": 4.53, "N2": 0.0},
        units=units.dimensionless,
        doc="Exponential factor for saturation capacity",
    )
    blk.affinity_characteristic_energy = Param(
        blk.isotherm_components,
        initialize={"CO2": 62.2, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="Characteristic energy for the affinity parameter",
    )


def Toth_isotherm(blk, i, pressure, temperature):
    """
    Method to add isotherm for components.
    Isotherm equation: Toth isotherm

    NOTE: the affinity of the material towards CO2 is not reduced by
          neither N2 nor H2O. N2 adsorption is negligible in
          this polystyrene functionalized with primary amine. Therefore,
          CO2 is considered as the only adsorbing component

    Keyword Arguments:
        i : component
        pressure : partial pressure of components
        temperature : temperature

    """

    T = temperature
    p = {}
    loading = {}

    for j in blk.isotherm_components:
        p[j] = units.convert(pressure[j], to_units=units.bar)

    if i == "CO2":

        saturation_capacity = blk.saturation_capacity_ref[i] * exp(
            blk.saturation_capacity_exponential_factor[i]
            * (1 - T / blk.temperature_ref)
        )
        affinity = blk.affinity_preexponential_factor[i] * exp(
            units.convert(
                blk.affinity_characteristic_energy[i], to_units=units.J / units.mol
            )
            / const.gas_constant
            / blk.temperature_ref
            * (blk.temperature_ref / T - 1)
        )
        toth_constant = blk.toth_constant_ref[i] + blk.toth_constant_alpha[i] * (
            1 - blk.temperature_ref / T
        )

        loading[i] = (
            saturation_capacity
            * affinity
            * p[i]
            / (1 + (affinity * p[i]) ** toth_constant) ** (1 / toth_constant)
        )

    elif i == "N2":
        # no adsorption is assumed of N2 in this adsorbent
        loading[i] = 1e-10 * units.mol / units.kg

    return loading[i]
