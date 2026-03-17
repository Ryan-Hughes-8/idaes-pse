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
    # log,
    # PositiveReals,
    # NonPositiveReals,
    # TransformationFactory,
    # units,
    # Block,
)


def add_dual_site_Langmuir_parameters(blk):
    """
    Method for adding parameters of the dual site Langmuir isotherm model.
    """

    blk.dh_ads = Param(
        blk.isotherm_components,
        initialize={"CO2": -37000, "N2": 0},
        units=units.J / units.mol,
        doc="Heat of adsorption",
    )
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
