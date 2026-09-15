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

from enum import Enum
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


class IsothermModel(Enum):
    """
    Enum for supported isotherm models to use with custom adsorbent
    """

    Langmuir = 1
    dual_site_Langmuir = 2
    weighted_DSL = 3
    extended_Sips = 4
    Toth = 5
    Henry = 6
    Langmuir_Freundlich = 7
    Sips = 8
    constant = 9


def add_parameters_custom_isotherm(blk):
    """
    upper level function to add isotherm model parameters
    """
    # need to group components into sets based on shared isotherm models
    grouped = {}

    for comp, model in blk.config.isotherm_models.items():
        grouped.setdefault(model, []).append(comp)

    for model, comp_set in grouped.items():
        print(f"adding {model} for component(s) {comp_set}")
        if model == IsothermModel.Langmuir_Freundlich:
            add_Langmuir_Freundlich_parameters(blk, comp_set)
        elif model == IsothermModel.Henry:
            add_Henry_parameters(blk, comp_set)
        elif model == IsothermModel.Langmuir:
            add_Langmuir_parameters(blk, comp_set)
        elif model == IsothermModel.dual_site_Langmuir:
            add_dual_site_Langmuir_parameters(blk, comp_set)
        elif model == IsothermModel.extended_Sips:
            add_extended_Sips_parameters(blk, comp_set)
        elif model == IsothermModel.Sips:
            add_Sips_parameters(blk, comp_set)
        elif model == IsothermModel.weighted_DSL:
            add_weighted_DSL_parameters(blk, comp_set)
        elif model == IsothermModel.Toth:
            add_Toth_parameters(blk, comp_set)
        elif model == IsothermModel.constant:
            add_constant_parameters(blk, comp_set)


def custom_isotherm(blk, i, pressure, temperature):
    """
    upper level function to add isotherm models and their parameters to blk.
    """
    model = blk.config.isotherm_models[i]
    if model == IsothermModel.Langmuir_Freundlich:
        return Langmuir_Freundlich_isotherm(blk, i, pressure, temperature)
    elif model == IsothermModel.Henry:
        return Henry_isotherm(blk, i, pressure, temperature)
    elif model == IsothermModel.Langmuir:
        return Langmuir_isotherm(blk, i, pressure, temperature)
    elif model == IsothermModel.dual_site_Langmuir:
        return dual_site_Langmuir_isotherm(blk, i, pressure, temperature)
    elif model == IsothermModel.extended_Sips:
        return extended_Sips_isotherm(blk, i, pressure, temperature)
    elif model == IsothermModel.Sips:
        return Sips_isotherm(blk, i, pressure, temperature)
    elif model == IsothermModel.weighted_DSL:
        return weighted_DSL_isotherm(blk, i, pressure, temperature)
    elif model == IsothermModel.Toth:
        return Toth_isotherm(blk, i, pressure, temperature)
    elif model == IsothermModel.constant:
        return constant_isotherm(blk, i, pressure, temperature)


def add_Langmuir_Freundlich_parameters(blk, comp_set):
    """
    Method for adding parameters of the Langmuir isotherm model for specified component.
    """

    blk.LF_qsat0 = Param(
        comp_set,
        initialize=5,
        units=units.mol / units.kg,
        doc="Langmuir-Freundlich saturation capacity reference [mol/kg] or [mmol/g]",
    )
    blk.LF_chi = Param(
        comp_set,
        initialize=2,
        units=units.dimensionless,
        doc="Langmuir-Freundlich saturation capacity chi",
    )
    blk.LF_b0 = Param(
        comp_set,
        initialize=1e-10,
        units=units.Pa**-1,
        doc="Langmuir-Freundlich b0 (pre-exponential) bar^-1",
    )
    blk.LF_dH = Param(
        comp_set,
        initialize=-35.0 * 1e3,
        units=units.J / units.mol,
        doc="Langmuir-Freundlich dH [J/mol]",
    )
    blk.LF_nu0 = Param(
        comp_set,
        initialize=0.8,
        units=units.dimensionless,
        doc="Langmuir-Freundlich nu reference",
    )
    blk.LF_c = Param(
        comp_set,
        initialize=-0.3,
        units=units.dimensionless,
        doc="Langmuir-Freundlich nu c",
    )
    blk.LF_T0 = Param(
        comp_set,
        initialize=298.15,
        units=units.K,
        doc="Langmuir-Freundlich reference T [K]",
    )


def Langmuir_Freundlich_isotherm(blk, i, pressure, temperature):
    """
    Method to add isotherm for components.
    Isotherm equation: Langmuir-Freundlich

    Keyword Arguments:
        i : component
        pressure : dict containing partial pressure of components
        temperature : temperature

    """

    T = temperature
    p = units.convert(pressure[i], to_units=units.Pa)

    b = blk.LF_b0[i] * exp(-blk.LF_dH[i] / const.gas_constant / T)
    q_sat = blk.LF_qsat0[i] * exp(blk.LF_chi[i] * (1 - T / blk.LF_T0[i]))
    nu = blk.LF_nu0[i] + blk.LF_c[i] * (1 - blk.LF_T0[i] / T)

    loading = q_sat * b * p**nu / (1 + b * p**nu)

    return loading


def add_Henry_parameters(blk, comp_set):
    """
    Method for adding parameters of the Langmuir isotherm model.
    """

    blk.Henry_a0 = Param(
        comp_set,
        initialize=4.57e-9,
        units=units.mmol / units.g / units.Pa,
        doc="Henry b0 (pre-exponential) [mmol/g/bar] or [mol/kg/bar]",
    )
    blk.Henry_dH = Param(
        comp_set,
        initialize=-2.12e4,
        units=units.J / units.mol,
        doc="Henry dh [J/mol]",
    )


def Henry_isotherm(blk, i, pressure, temperature):
    """
    Method to add isotherm for components.
    Isotherm equation: Henry

    Keyword Arguments:
        i : component
        pressure : dict containing partial pressure of components
        temperature : temperature

    """

    T = temperature
    p = units.convert(pressure[i], to_units=units.Pa)

    a = blk.Henry_a0[i] * exp(-blk.Henry_dH[i] / const.gas_constant / T)

    loading = a * p

    return loading


def add_Langmuir_parameters(blk, comp_set):
    """
    Method for adding parameters of the Langmuir isotherm model.
    """

    blk.Langmuir_q_sat0 = Param(
        comp_set,
        initialize=3.0,
        units=units.mol / units.kg,
        doc="Langmuir saturation capacity [mol/kg] or [mmol/g]",
    )
    blk.Langmuir_b0 = Param(
        comp_set,
        initialize=1e-10,
        units=units.Pa**-1,
        doc="Langmuir b0 (pre-exponential)",
    )
    blk.Langmuir_dH = Param(
        comp_set,
        initialize=-3.5e4,
        units=units.J / units.mol,
        doc="Langmuir E",
    )
    blk.Langmuir_chi = Param(
        comp_set,
        initialize=0.69,
        units=units.dimensionless,
        doc="Langmuir chi",
    )
    blk.Langmuir_T0 = Param(
        comp_set,
        initialize=298.15,
        units=units.K,
        doc="Langmuir reference T",
    )


def Langmuir_isotherm(blk, i, pressure, temperature):
    """
    Method to add isotherm for components.
    Isotherm equation: Langmuir

    Keyword Arguments:
        i : component
        pressure : partial pressure of components
        temperature : temperature

    """

    T = temperature
    p = units.convert(pressure[i], to_units=units.Pa)

    b = blk.Langmuir_b0[i] * exp(-blk.Langmuir_dH[i] / const.gas_constant / T)
    q = blk.Langmuir_q_sat0[i] * exp(blk.Langmuir_chi[i] * (1 - T / blk.Langmuir_T0[i]))

    loading = q * b * p / (1 + b * p)

    return loading


def add_dual_site_Langmuir_parameters(blk, comp_set):
    """
    Method for adding parameters of the dual site Langmuir isotherm model.
    """

    blk.temperature_ref = Param(
        initialize=298.15,
        units=units.K,
        doc="Reference temperature",
    )
    blk.saturation_capacity_site_b = Param(
        comp_set,
        initialize={"CO2": 2.387, "N2": 0.0},
        units=units.mol / units.kg,
        doc="saturation capacity at site b",
    )
    blk.saturation_capacity_site_d = Param(
        comp_set,
        initialize={"CO2": 3.2711, "N2": 0.0},
        units=units.mol / units.kg,
        doc="saturation capacity at site d",
    )
    blk.dual_site_langmuir_constant_pre_exp_b = Param(
        comp_set,
        initialize={"CO2": 5.519e-7, "N2": 0.0},
        units=units.meter**3 / units.mol,
        doc="dual site langmuir constant for site b",
    )
    blk.dual_site_langmuir_constant_pre_exp_d = Param(
        comp_set,
        initialize={"CO2": 5.187e-08, "N2": 0.0},
        units=units.meter**3 / units.mol,
        doc="dual site langmuir constant for site d",
    )
    blk.internal_energy_b = Param(
        comp_set,
        initialize={"CO2": -35.06, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="internal energy of site b",
    )
    blk.internal_energy_d = Param(
        comp_set,
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


def add_extended_Sips_parameters(blk, comp_set):
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
        comp_set,
        initialize={"CO2": 7.268, "N2": 4.051},
        units=units.mol / units.kg,
        doc="Saturation capacity at reference temperature",
    )
    blk.saturation_capacity_exponential_factor = Param(
        comp_set,
        initialize={"CO2": -0.61684, "N2": 0.0},
        units=units.dimensionless,
        doc="Isotherm fitting parameter",
    )
    blk.affinity_parameter_preexponential_factor = Param(
        comp_set,
        initialize={"CO2": 1.129e-4, "N2": 5.8470e-5},
        units=units.bar**-1,
        doc="Affinity parameter pre-exponential factor",
    )
    blk.affinity_parameter_characteristic_energy = Param(
        comp_set,
        initialize={"CO2": 28.389, "N2": 18.4740},
        units=units.kJ / units.mol,
        doc="Characteristic energy for the affinity parameter",
    )
    blk.heterogeneity_parameter_ref = Param(
        comp_set,
        initialize={"CO2": 0.42456, "N2": 0.98624},
        units=units.dimensionless,
        doc="Heterogeneity parameter at reference temperature",
    )
    blk.heterogeneity_parameter_alpha = Param(
        comp_set,
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


def add_Sips_parameters(blk, comp_set):
    """
    Method to add sips isotherm parameters. CO2 and N2 are adsorbed.

    """

    blk.Sips_T0 = Param(
        comp_set,
        initialize=298.15,
        units=units.K,
        doc="Reference temperature",
    )
    blk.Sips_qsat0 = Param(
        comp_set,
        initialize=12.76,
        units=units.mol / units.kg,
        doc="Sips reference capacity [mol/kg]",
    )
    blk.Sips_chi = Param(
        comp_set,
        initialize=-0.1994,
        units=units.dimensionless,
        doc="Sips chi",
    )
    blk.Sips_b0 = Param(
        comp_set,
        initialize=3.396e-10,
        units=units.Pa**-1,
        doc="Sips b0 [Pa^-1]",
    )
    blk.Sips_dH = Param(
        comp_set,
        initialize=-2.358e4,
        units=units.J / units.mol,
        doc="Sips energy [J/mol]",
    )
    blk.Sips_nu0 = Param(
        comp_set,
        initialize=0.9357,
        units=units.dimensionless,
        doc="Sips nu reference",
    )
    blk.Sips_c = Param(
        comp_set,
        initialize=0.1314,
        units=units.dimensionless,
        doc="Sips parameter c",
    )


def Sips_isotherm(blk, i, pressure, temperature):
    """
    Method to add isotherm for components.
    Isotherm equation: Extended Sips

    Keyword Arguments:
        i : component
        pressure : partial pressure of components
        temperature : temperature

    """
    T = temperature
    p = units.convert(pressure[i], to_units=units.Pa)

    b = blk.Sips_b0[i] * exp(-blk.Sips_dH[i] / const.gas_constant / T)
    q_sat = blk.Sips_qsat0[i] * exp(blk.Sips_chi[i] * (1 - T / blk.Sips_T0[i]))
    nu = blk.Sips_nu0[i] + blk.Sips_c[i] * (1 - blk.Sips_T0[i] / T)

    b_p = b * p
    loading = q_sat * b_p ** (1 / nu) / (1 + b_p ** (1 / nu))

    return loading


def add_weighted_DSL_parameters(blk, comp_set):
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
        comp_set,
        initialize={"CO2": 0.146, "N2": 0.0},
        units=units.mol / units.kg,
        doc="Lower isotherm saturation capacity",
    )
    blk.upper_saturtion_capacity = Param(
        comp_set,
        initialize={"CO2": 3.478, "N2": 0.0},
        units=units.mol / units.kg,
        doc="Upper isotherm saturation capacity",
    )
    blk.lower_affinity_preexponential_factor = Param(
        comp_set,
        initialize={"CO2": 0.009, "N2": 0.0},
        units=units.bar**-1,
        doc="Pre-exponential factor for the lower isotherm affinity parameter",
    )
    blk.upper_affinity_preexponential_factor_1 = Param(
        comp_set,
        initialize={"CO2": 9.00e-07, "N2": 0.0},
        units=units.bar**-1,
        doc="Pre-exponential factor for the upper isotherm site 1 affinity parameter",
    )
    blk.upper_affinity_preexponential_factor_2 = Param(
        comp_set,
        initialize={"CO2": 5.00e-04, "N2": 0.0},
        units=units.mol / units.kg / units.bar,
        doc="Pre-exponential factor for the upper isotherm site 2 affinity parameter",
    )
    blk.lower_affinity_characteristic_energy = Param(
        comp_set,
        initialize={"CO2": 31.0, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="Characteristic energy for the lower isotherm affinity parameter",
    )
    blk.upper_affinity_characteristic_energy_1 = Param(
        comp_set,
        initialize={"CO2": 59.0, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="Characteristic energy for the upper isotherm site 1 affinity parameter",
    )
    blk.upper_affinity_characteristic_energy_2 = Param(
        comp_set,
        initialize={"CO2": 18.0, "N2": 0.0},
        units=units.kJ / units.mol,
        doc="Characteristic energy for the upper isotherm site 2 affinity parameter",
    )
    blk.step_width_preexponential_factor = Param(
        comp_set,
        initialize={"CO2": 1.24e-01, "N2": 0.0},
        units=units.dimensionless,
        doc="Pre-exponential factor for step width",
    )
    blk.step_width_exponential_factor = Param(
        comp_set,
        initialize={"CO2": 0.0, "N2": 0.0},
        units=units.dimensionless,
        doc="Exponential factor for step width",
    )
    blk.weighting_function_exponent = Param(
        comp_set,
        initialize={"CO2": 4.00, "N2": 0.0},
        units=units.dimensionless,
        doc="Isotherm weighting function exponent",
    )
    blk.step_partial_pressure_ref = Param(
        comp_set,
        initialize={"CO2": 0.5 * 1e-3, "N2": 0.0},
        units=units.bar,
        doc="Step partial pressure at reference temperature",
    )
    blk.step_enthalpy = Param(
        comp_set,
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


def add_Toth_parameters(blk, comp_set):
    """
    Method to add adsorbent related parameters to run fixed bed TSA model.
    This method is to add parameters for polystyrene functionalized
    with primary amine.

    Elfvinga, J.; Bajamundia, C.; Kauppinena, J.; Sainiob, T. Modelling
    of equilibrium working capacity of PSA, TSA and TVSA processes for
    CO2 adsorption under direct air capture conditions. Journal of CO2
    Utilization, 22, 2017.

    """

    blk.Toth_T0 = Param(
        comp_set,
        initialize=298.15,
        units=units.K,
        doc="Toth Reference temperature [K]",
    )
    blk.Toth_q_sat0 = Param(
        comp_set,
        initialize=1.71,
        units=units.mol / units.kg,
        doc="Toth saturation capacity at reference temperature [mol/kg] or [mmol/g]",
    )
    blk.Toth_b0 = Param(
        comp_set,
        initialize=1e-10,
        units=units.Pa**-1,
        doc="Toth pre-exponential factor for the affinity parameter [Pa^-1]",
    )
    blk.Toth_nu0 = Param(
        comp_set,
        initialize=0.75,
        units=units.dimensionless,
        doc="Toth constant at reference temperature",
    )
    blk.Toth_c = Param(
        comp_set,
        initialize=0.601,
        units=units.dimensionless,
        doc="Toth nu temperature dependant parameter",
    )
    blk.Toth_chi = Param(
        comp_set,
        initialize=4.53,
        units=units.dimensionless,
        doc="Exponential factor for saturation capacity",
    )
    blk.Toth_dH = Param(
        comp_set,
        initialize=4.5e4,
        units=units.J / units.mol,
        doc="Toth dH [J/mol]",
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
    p = units.convert(pressure[i], to_units=units.Pa)

    q_sat = blk.Toth_q_sat0[i] * exp(blk.Toth_chi[i] * (1 - T / blk.Toth_T0[i]))
    b = blk.Toth_b0[i] * exp(-blk.Toth_dH[i] / const.gas_constant / T)
    nu = blk.Toth_nu0[i] + blk.Toth_c[i] * (1 - blk.Toth_T0[i] / T)

    loading = q_sat * b * p / (1 + (b * p) ** nu) ** (1 / nu)

    return loading


def add_constant_parameters(blk, comp_set):
    """
    no parameters needed
    """
    pass


def constant_isotherm(blk, i, pressure, temperature, value=1e-10):
    """
    only need to return value. Not dependent on comp, temp, or pressure.
    """

    return value * units.mol / units.kg
