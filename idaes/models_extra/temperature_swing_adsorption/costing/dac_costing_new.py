#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
""" """
import json

# import pytest
import pandas as pd
from pyomo.environ import (
    Block,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Param,
    units,
    value,
    Var,
    Expression,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.config import ConfigValue

from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom


from idaes.models.properties import iapws95


from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)
import pyomo.environ as pyo

vacuum_pump_params = {
    "9": {
        "B": {
            "vac_2.9_eq1": {
                "Account Name": "Vacuum pump quote for 2.9 psia",
                "Exponent": 0.45,
                "Process Parameter": "CO2 Flowrate",
                "BEC": 47.8e6 / 1e3,  # 46e6/1e3,  # BEC in $1000 dollars
                "BEC_units": "K$2018",
                "Eng Fee": 0.2,
                "Process Contingency": 0.18,
                "Project Contingency": 0.2,
                "RP Value": 77.00,
                "Units": "MW",
            },
            "vac_6.5_eq1": {
                "Account Name": "Vacuum pump quote for 6.5 psia",
                "Exponent": 0.45,
                "Process Parameter": "Absorber volume",
                "BEC": 18.71e3,  # 18e6 / 1e3,
                "BEC_units": "K$2018",
                "Eng Fee": 0.2,
                "Process Contingency": 0.18,
                "Project Contingency": 0.2,
                "RP Value": 63.00,
                "Units": "MW",
            },
        }
    }
}


def get_dac_costing_data(case):
    if case == "electric_boiler":
        fname = "dac_eb_costing_data.json"
    elif case == "NGCC":
        fname = "dac_retrofit_costing_data.json"
    else:
        print("Invalid case")

    with open(fname, "r") as f:
        costing_data = json.load(f)

    return costing_data


def build_dac_costing(dac, dll_per_brick=100, replacement_time=3.5):
    # capital costs
    dac_cost_params = get_dac_costing_data(dac.config.steam_source)

    dac.costing = QGESSCosting()
    CE_index_year = "2018"
    CE_index_units = getattr(units, "MUSD_" + CE_index_year)

    # add a var to account for multiple DAC units
    dac.number_of_units = Var(initialize=1, bounds=(0, 100))
    dac.number_of_units.fix(1)

    if dac.config.steam_source == "electric_boiler":
        dac.raw_water_system = UnitModelBlock()
        dac.raw_water_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=dac.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["3.2", "3.4", "9.5", "14.6"],
                "scaled_param": dac.raw_water_withdrawal[0] * dac.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": dac_cost_params,
            },
        )

        dac.steam_system = UnitModelBlock()
        dac.steam_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=dac.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["3.1", "3.3", "3.5"],
                "scaled_param": dac.steam_flow_mass[0] * dac.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": dac_cost_params,
            },
        )

        dac.cooling_tower = UnitModelBlock()
        dac.cooling_tower.costing = UnitModelCostingBlock(
            flowsheet_costing_block=dac.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["9.1"],
                "scaled_param": dac.cooling_tower_duty[0] * dac.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": dac_cost_params,
            },
        )

        dac.water_discharge_system = UnitModelBlock()
        dac.water_discharge_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=dac.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["3.7"],
                "scaled_param": dac.process_water_discharge[0] * dac.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": dac_cost_params,
            },
        )

        dac.cooling_water_system = UnitModelBlock()
        dac.cooling_water_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=dac.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["9.2", "9.3", "9.4", "9.6", "9.7", "14.5"],
                "scaled_param": dac.circulating_water_flowrate[0] * dac.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": dac_cost_params,
            },
        )

        dac.electric_systems = UnitModelBlock()
        dac.electric_systems.costing = UnitModelCostingBlock(
            flowsheet_costing_block=dac.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": [
                    "11.1",
                    "11.2",
                    "11.3",
                    "11.4",
                    "11.5",
                    "11.6",
                    "11.7",
                    "11.8",
                    "11.9",
                    "12.4",
                    "12.5",
                    "12.6",
                    "12.7",
                    "12.8",
                    "12.9",
                    "13.1",
                    "13.2",
                    "13.3",
                    "14.4",
                    "14.7",
                    "14.8",
                    "14.9",
                    "14.10",
                ],
                "scaled_param": dac.auxiliary_load[0] * dac.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": dac_cost_params,
            },
        )

        # Electric Boiler 15.9
        dac.electric_boiler = UnitModelBlock()
        dac.electric_boiler.costing = UnitModelCostingBlock(
            flowsheet_costing_block=dac.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["15.9"],
                "scaled_param": dac.steam_flow_mass[0],
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": dac_cost_params,
            },
        )

    elif dac.config.steam_source == "NGCC":
        dac.steam_flow_system = UnitModelBlock()
        dac.steam_flow_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=dac.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["8.4"],
                "scaled_param": dac.steam_flow_mass[0] * dac.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": dac_cost_params,
            },
        )

        dac.electric_systems = UnitModelBlock()
        dac.electric_systems.costing = UnitModelCostingBlock(
            flowsheet_costing_block=dac.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["11.2", "11.3", "11.4", "11.5", "11.6"],
                "scaled_param": dac.auxiliary_load[0] * dac.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": dac_cost_params,
            },
        )

    dac.vessels = UnitModelBlock()
    dac.vessels.costing = UnitModelCostingBlock(
        flowsheet_costing_block=dac.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.1"],
            "scaled_param": units.convert(dac.bed_volume, to_units=units.ft**3),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": dac_cost_params,
        },
    )

    purge_pressure = dac.config.purge_pressure
    # 15.2 - DAC CO2 Compression & Drying
    dac.product_compression = UnitModelBlock()
    if purge_pressure == 1:
        dac.product_compression.costing = UnitModelCostingBlock(
            flowsheet_costing_block=dac.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["15.2"],
                "scaled_param": dac.compressor_power[0] * dac.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": dac_cost_params,
            },
        )

    elif purge_pressure == 0.2:
        dac.product_compression.costing = UnitModelCostingBlock(
            flowsheet_costing_block=dac.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["vac_2.9_eq1"],
                "scaled_param": dac.compressor_power[0] * dac.number_of_units,
                "tech": 9,
                "ccs": "B",
                "additional_costing_params": vacuum_pump_params,
            },
        )

    elif purge_pressure == 0.5:
        dac.product_compression.costing = UnitModelCostingBlock(
            flowsheet_costing_block=dac.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["vac_6.5_eq1"],
                "scaled_param": dac.compressor_power[0] * dac.number_of_units,
                "units": "MW",
                "tech": 9,
                "ccs": "B",
                "additional_costing_params": vacuum_pump_params,
            },
        )

    # 15.3 - DAC CO2 Compressor Aftercooler
    dac.compressor_aftercooler = UnitModelBlock()
    dac.compressor_aftercooler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=dac.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.3"],
            "scaled_param": dac.compressor_aftercooler_heat_duty[0]
            * dac.number_of_units,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": dac_cost_params,
        },
    )

    # 15.4 - DAC System Air Handling Duct and Dampers (1 system per 2 beds)
    dac.duct_dampers = UnitModelBlock()
    dac.duct_dampers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=dac.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.4"],
            "scaled_param": (dac.air_flow_mass[0] / dac.num_beds_total * 2),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": dac_cost_params,
        },
    )

    # 15.5 - DAC System Air Handling Fans (1 system per 2 beds)
    dac.feed_fans = UnitModelBlock()
    dac.feed_fans.costing = UnitModelCostingBlock(
        flowsheet_costing_block=dac.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.5"],
            "scaled_param": dac.fans.work_mechanical[0] / dac.num_beds_total * 2,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": dac_cost_params,
        },
    )

    # 15.6 - DAC Desorption Process Gas Handling System (pure CO2 gas)
    dac.desorption_gas_handling = UnitModelBlock()
    dac.desorption_gas_handling.costing = UnitModelCostingBlock(
        flowsheet_costing_block=dac.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.6"],
            "scaled_param": dac.co2_flow_mass[0],
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": dac_cost_params,
        },
    )

    # 15.7 - DAC Steam Distribution System
    dac.steam_distribution = UnitModelBlock()
    dac.steam_distribution.costing = UnitModelCostingBlock(
        flowsheet_costing_block=dac.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.7"],
            "scaled_param": dac.steam_flow_mass[0],
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": dac_cost_params,
        },
    )

    # 15.8 - DAC System Controls Equipment
    dac.controls_equipment = UnitModelBlock()
    dac.controls_equipment.costing = UnitModelCostingBlock(
        flowsheet_costing_block=dac.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.8"],
            "scaled_param": dac.auxiliary_load[0],
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": dac_cost_params,
        },
    )

    # # Electric Boiler 15.9
    # dac.electric_boiler = UnitModelBlock()
    # dac.electric_boiler.costing = UnitModelCostingBlock(
    #     flowsheet_costing_block=dac.costing,
    #     costing_method=QGESSCostingData.get_PP_costing,
    #     costing_method_arguments={
    #         "cost_accounts": ["15.9"],
    #         "scaled_param": dac.steam_flow_mass[0],
    #         "tech": 8,
    #         "ccs": "B",
    #         "additional_costing_params": dac_cost_params,
    #     },
    # )

    # we need a custom method of calculating the total TPC to account for distributed vs centralized systems
    distributed_systems = [
        dac.vessels.costing,
        dac.duct_dampers.costing,
        dac.feed_fans.costing,
        dac.desorption_gas_handling.costing,
        dac.controls_equipment.costing,
    ]

    centralized_TPCs = []

    for b in dac.costing._registered_unit_costing:
        if b not in distributed_systems:
            for key in b.total_plant_cost.keys():
                centralized_TPCs.append(b.total_plant_cost[key])

    dac.costing.total_TPC = Var(
        initialize=100,
        bounds=(0, 1e4),
        doc="total TPC in $MM",
    )

    @dac.costing.Constraint()
    def total_TPC_eq(b):
        return (
            b.total_TPC
            == sum(centralized_TPCs)
            + (
                dac.vessels.costing.total_plant_cost["15.1"] / 120 * dac.num_beds_total
                + dac.duct_dampers.costing.total_plant_cost["15.4"]
                * dac.num_beds_total
                / 2
                + dac.feed_fans.costing.total_plant_cost["15.5"]
                * dac.num_beds_total
                / 2
                + dac.desorption_gas_handling.costing.total_plant_cost["15.6"]
                + dac.controls_equipment.costing.total_plant_cost["15.8"]
            )
            * dac.number_of_units
        )

    # resorces to be costed
    resources = [
        "water",
        "water_treatment_chemicals",
        "aux_power",
        "sorbent",
        "waste_sorbent",
        "boiler_feed_water",
    ]
    if dac.config.steam_source == "NGCC":
        resources.append("IP_steam")

    # resource consumption rates for variable costing
    capacity_factor = 0.85

    @dac.costing.Expression(dac.time)
    def water_use(b, t):
        # water use is lineraly scaled from NETL reference
        ref_water = 60792.407 * units.gal / units.day
        ref_air = 1629629 * units.kmol / units.hr
        return dac.air_flow_mol[t] * (ref_water / ref_air) * dac.number_of_units

    @dac.costing.Expression(dac.time)
    def water_treatment_chems(b, t):
        # treatment chemical use is lineraly scaled from NETL reference
        ref_chem = 0.1811 * units.ton / units.day
        ref_air = 1629629 * units.kmol / units.hr
        return dac.air_flow_mol[t] * (ref_chem / ref_air) * dac.number_of_units

    @dac.costing.Expression(dac.time)
    def energy_purchased(b, t):  # in kWh/day
        hr_per_day = 24 * units.hr / units.day
        return dac.auxiliary_load[t] * hr_per_day * dac.number_of_units

    @dac.costing.Expression(dac.time)
    def sorbent_rate(b, t):
        return (
            units.convert(dac.bed_volume, to_units=units.ft**3)
            * (1 - dac.bed_voidage)
            * dac.num_beds_total
            / dac.config.sorbent_lifespan
            / 365
            / units.day
        )

    @dac.costing.Expression(dac.time)
    def bfw_rate(b, t):
        hr_per_day = 24 * units.hr / units.day
        return dac.BFW_makeup[t] * hr_per_day

    @dac.costing.Expression(dac.time)
    def steam_rate(b, t):
        hr_per_day = 24 * units.hr / units.day
        return (
            units.convert(dac.steam_flow_mass[t], to_units=units.kg / units.hr)
            * hr_per_day
        )

    # vars for resource consumption rates
    rates = [
        dac.costing.water_use,
        dac.costing.water_treatment_chems,
        dac.costing.energy_purchased,
        dac.costing.sorbent_rate,
        dac.costing.sorbent_rate,
        dac.costing.bfw_rate,
    ]
    if dac.config.steam_source == "NGCC":
        rates.append(dac.costing.steam_rate)

    # resource prices
    prices = {
        "sorbent": 4 * units.USD_2018 / units.ft**3,  # 4 or 100
        "aux_power": 0.06 * units.USD_2018 / units.kWh,
        "waste_sorbent": 0.86 * units.USD_2018 / units.ft**3,
        "IP_steam": 0.00733 * units.USD_2018 / units.kg,
        "boiler_feed_water": 2.45
        / 1000
        * units.USD_2018
        / units.kg,  # Turton et al., 2012
    }

    @dac.costing.Expression()
    def land_cost1(b):
        return (
            156000 * (dac.num_beds_total / 120) ** (0.78)
        ) * 1e-6  # scaled to Millions

    @dac.costing.Expression()
    def tonne_CO2_capture(b):
        # return dac.co2_flow_mass[0] / 2204.62 * 8760 * dac.number_of_units
        return (
            units.convert(dac.co2_flow_mass[0], to_units=units.tonne / units.year)
            * dac.number_of_units
        )

    dac.costing.build_process_costs(
        net_power=None,
        # arguments related to fixed OM costs
        total_plant_cost=True,
        labor_rate=38.50,
        labor_burden=30,
        operators_per_shift=8,
        tech=6,
        fixed_OM=True,
        # arguments related owners costs
        variable_OM=True,
        capacity_factor=capacity_factor,
        land_cost=dac.costing.land_cost1,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel=None,
        waste=None,
        chemicals=None,
        tonne_CO2_capture=dac.costing.tonne_CO2_capture,
    )

    # we want the number of operators to depend on the number of units
    dac.costing.operators_per_shift_var = Var(initialize=8, bounds=(0, 200))

    @dac.costing.Constraint()
    def operator_eqn(b):
        return b.operators_per_shift_var == 8 * dac.number_of_units

    dac.costing.annual_labor_cost_rule.deactivate()

    @dac.costing.Constraint()
    def annual_labor_cost_rule_new(c):
        return c.annual_operating_labor_cost == units.convert(
            (
                c.operators_per_shift_var
                * c.labor_rate
                * (1 + c.labor_burden / 100)
                * 8760
                * units.hr
            ),
            CE_index_units,
        )

    # # brick replacement variable cost
    # dac.costing.other_variable_costs.unfix()

    # @dac.costing.Constraint()
    # def other_var_costs_eqn(b):
    #     return (
    #         b.other_variable_costs[0]
    #         == dac.number_of_units
    #         * dll_per_brick
    #         * dac.number_of_panels[0]
    #         * dac.bricks_per_panel
    #         * 1e-6
    #         / replacement_time
    #     )

    dac.costing.costing_initialization()


"""
Centralized Systems
- Water withdrawal and pretreating
- Boiler feedwater system
- Cooling tower and cooling water circulation
- Waste water treatment
- Electric equipment (Switchyard, transformers, ..)
- Instrumentation and control equipment
- Site preparation and facilities
- Buildings
- Electric boiler
- CO2 compressor/vacuum pump

Distributed Systems
- DAC units (costed per brick)
- Air ducts
- Fans
- Desorption gas handling
- Steam distribution
- DAC specific controls

The distributed systems need to be multiplied by the number of units when added
to the total plant cost.
"""
