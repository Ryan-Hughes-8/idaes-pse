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
import os
import json
import textwrap
from sys import stdout
from pandas import DataFrame

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
from pyomo.common.fileutils import this_file_dir
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.config import ConfigValue

from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog
from idaes.core.util.tables import stream_table_dataframe_to_string

from idaes.models.properties import iapws95


from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)
import pyomo.environ as pyo

directory = this_file_dir()
_log = idaeslog.getLogger(__name__)

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
        fname = "costing_params_dac_electric_boiler.json"
    elif case == "retrofit_ngcc":
        fname = "costing_params_dac_retrofit_ngcc.json"
    else:
        print("Invalid case")

    # load custom costing parameters
    with open(os.path.join(directory, fname), "r") as f:
        costing_params = json.load(f)

    return costing_params


def get_dac_costing(unit, costing_case):
    # capital costs
    costing_params = get_dac_costing_data(costing_case)

    unit.costing = QGESSCosting()
    CE_index_year = "2018"
    CE_index_units = getattr(units, "MUSD_" + CE_index_year)

    # add a var to account for multiple DAC units
    unit.number_of_units = Var(initialize=1, bounds=(0, 100))
    unit.number_of_units.fix(1)

    if costing_case == "electric_boiler":
        unit.raw_water_system = UnitModelBlock()
        unit.raw_water_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["3.2", "3.4", "9.5", "14.6"],
                "scaled_param": unit.raw_water_withdrawal[0] * unit.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        unit.steam_system = UnitModelBlock()
        unit.steam_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["3.1", "3.3", "3.5"],
                "scaled_param": unit.steam_flow_mass[0] * unit.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        unit.cooling_tower = UnitModelBlock()
        unit.cooling_tower.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["9.1"],
                "scaled_param": unit.cooling_tower_duty[0] * unit.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        unit.water_discharge_system = UnitModelBlock()
        unit.water_discharge_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["3.7"],
                "scaled_param": unit.process_water_discharge[0] * unit.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        unit.cooling_water_system = UnitModelBlock()
        unit.cooling_water_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["9.2", "9.3", "9.4", "9.6", "9.7", "14.5"],
                "scaled_param": unit.circulating_water_flowrate[0]
                * unit.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        unit.electric_systems = UnitModelBlock()
        unit.electric_systems.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.costing,
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
                "scaled_param": unit.auxiliary_load[0] * unit.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        # Electric Boiler 15.9
        unit.electric_boiler = UnitModelBlock()
        unit.electric_boiler.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["15.9"],
                "scaled_param": unit.steam_flow_mass[0],
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

    elif costing_case == "retrofit_NGCC":
        unit.steam_flow_system = UnitModelBlock()
        unit.steam_flow_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["8.4"],
                "scaled_param": unit.steam_flow_mass[0] * unit.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        unit.electric_systems = UnitModelBlock()
        unit.electric_systems.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["11.2", "11.3", "11.4", "11.5", "11.6"],
                "scaled_param": unit.auxiliary_load[0] * unit.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

    unit.vessels = UnitModelBlock()
    unit.vessels.costing = UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.1"],
            "scaled_param": units.convert(unit.bed_volume, to_units=units.ft**3),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    purge_pressure = unit.config.purge_pressure
    # 15.2 - DAC CO2 Compression & Drying
    unit.product_compression = UnitModelBlock()
    if purge_pressure == 1:
        unit.product_compression.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["15.2"],
                "scaled_param": unit.compressor_power[0] * unit.number_of_units,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

    elif purge_pressure == 0.2:
        unit.product_compression.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["vac_2.9_eq1"],
                "scaled_param": unit.compressor_power[0] * unit.number_of_units,
                "tech": 9,
                "ccs": "B",
                "additional_costing_params": vacuum_pump_params,
            },
        )

    elif purge_pressure == 0.5:
        unit.product_compression.costing = UnitModelCostingBlock(
            flowsheet_costing_block=unit.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["vac_6.5_eq1"],
                "scaled_param": unit.compressor_power[0] * unit.number_of_units,
                "units": "MW",
                "tech": 9,
                "ccs": "B",
                "additional_costing_params": vacuum_pump_params,
            },
        )

    # 15.3 - DAC CO2 Compressor Aftercooler
    unit.compressor_aftercooler = UnitModelBlock()
    unit.compressor_aftercooler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.3"],
            "scaled_param": unit.compressor_aftercooler_heat_duty[0]
            * unit.number_of_units,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.4 - DAC System Air Handling Duct and Dampers (1 system per 2 beds)
    unit.duct_dampers = UnitModelBlock()
    unit.duct_dampers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.4"],
            "scaled_param": (unit.air_flow_mass[0] / unit.num_beds_total * 2),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.5 - DAC System Air Handling Fans (1 system per 2 beds)
    unit.feed_fans = UnitModelBlock()
    unit.feed_fans.costing = UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.5"],
            "scaled_param": unit.fans.work_mechanical[0] / unit.num_beds_total * 2,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.6 - DAC Desorption Process Gas Handling System (pure CO2 gas)
    unit.desorption_gas_handling = UnitModelBlock()
    unit.desorption_gas_handling.costing = UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.6"],
            "scaled_param": unit.co2_flow_mass[0],
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.7 - DAC Steam Distribution System
    unit.steam_distribution = UnitModelBlock()
    unit.steam_distribution.costing = UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.7"],
            "scaled_param": unit.steam_flow_mass[0],
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.8 - DAC System Controls Equipment
    unit.controls_equipment = UnitModelBlock()
    unit.controls_equipment.costing = UnitModelCostingBlock(
        flowsheet_costing_block=unit.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.8"],
            "scaled_param": unit.auxiliary_load[0],
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # # Electric Boiler 15.9
    # unit.electric_boiler = UnitModelBlock()
    # unit.electric_boiler.costing = UnitModelCostingBlock(
    #     flowsheet_costing_block=unit.costing,
    #     costing_method=QGESSCostingData.get_PP_costing,
    #     costing_method_arguments={
    #         "cost_accounts": ["15.9"],
    #         "scaled_param": unit.steam_flow_mass[0],
    #         "tech": 8,
    #         "ccs": "B",
    #         "additional_costing_params": costing_params,
    #     },
    # )

    # we need a custom method of calculating the total TPC to account for distributed vs centralized systems
    distributed_systems = [
        unit.vessels.costing,
        unit.duct_dampers.costing,
        unit.feed_fans.costing,
        unit.desorption_gas_handling.costing,
        unit.controls_equipment.costing,
    ]

    centralized_TPCs = []

    for b in unit.costing._registered_unit_costing:
        if b not in distributed_systems:
            for key in b.total_plant_cost.keys():
                centralized_TPCs.append(b.total_plant_cost[key])

    unit.costing.total_TPC = Var(
        initialize=100,
        bounds=(0, 1e4),
        doc="total TPC in $MM",
    )

    @unit.costing.Constraint()
    def total_TPC_eq(b):
        return (
            b.total_TPC
            == sum(centralized_TPCs)
            + (
                unit.vessels.costing.total_plant_cost["15.1"]
                / 120
                * unit.num_beds_total
                + unit.duct_dampers.costing.total_plant_cost["15.4"]
                * unit.num_beds_total
                / 2
                + unit.feed_fans.costing.total_plant_cost["15.5"]
                * unit.num_beds_total
                / 2
                + unit.desorption_gas_handling.costing.total_plant_cost["15.6"]
                + unit.controls_equipment.costing.total_plant_cost["15.8"]
            )
            * unit.number_of_units
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
    if costing_case == "retrofit_NGCC":
        resources.append("IP_steam")

    # resource consumption rates for variable costing
    capacity_factor = 0.85

    @unit.costing.Expression(unit.time)
    def water_use(b, t):
        # water use is lineraly scaled from NETL reference
        ref_water = 60792.407 * units.gal / units.day
        ref_air = 1629629 * units.kmol / units.hr
        return unit.air_flow_mol[t] * (ref_water / ref_air) * unit.number_of_units

    @unit.costing.Expression(unit.time)
    def water_treatment_chems(b, t):
        # treatment chemical use is lineraly scaled from NETL reference
        ref_chem = 0.1811 * units.ton / units.day
        ref_air = 1629629 * units.kmol / units.hr
        return unit.air_flow_mol[t] * (ref_chem / ref_air) * unit.number_of_units

    @unit.costing.Expression(unit.time)
    def energy_purchased(b, t):  # in kWh/day
        hr_per_day = 24 * units.hr / units.day
        return unit.auxiliary_load[t] * hr_per_day * unit.number_of_units

    @unit.costing.Expression(unit.time)
    def sorbent_rate(b, t):
        return (
            units.convert(unit.bed_volume, to_units=units.ft**3)
            * (1 - unit.bed_voidage)
            * unit.num_beds_total
            / unit.config.sorbent_lifespan
            / 365
            / units.day
        )

    @unit.costing.Expression(unit.time)
    def bfw_rate(b, t):
        hr_per_day = 24 * units.hr / units.day
        return unit.BFW_makeup[t] * hr_per_day

    @unit.costing.Expression(unit.time)
    def steam_rate(b, t):
        hr_per_day = 24 * units.hr / units.day
        return (
            units.convert(unit.steam_flow_mass[t], to_units=units.kg / units.hr)
            * hr_per_day
        )

    # vars for resource consumption rates
    rates = [
        unit.costing.water_use,
        unit.costing.water_treatment_chems,
        unit.costing.energy_purchased,
        unit.costing.sorbent_rate,
        unit.costing.sorbent_rate,
        unit.costing.bfw_rate,
    ]
    if costing_case == "retrofit_NGCC":
        rates.append(unit.costing.steam_rate)

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

    @unit.costing.Expression()
    def land_cost1(b):
        return (
            156000 * (unit.num_beds_total / 120) ** (0.78)
        ) * 1e-6  # scaled to Millions

    @unit.costing.Expression()
    def tonne_CO2_capture(b):
        # return unit.co2_flow_mass[0] / 2204.62 * 8760 * unit.number_of_units
        return (
            units.convert(unit.co2_flow_mass[0], to_units=units.tonne / units.year)
            * unit.number_of_units
        )

    unit.costing.build_process_costs(
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
        land_cost=unit.costing.land_cost1,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel=None,
        waste=None,
        chemicals=None,
        tonne_CO2_capture=unit.costing.tonne_CO2_capture,
    )

    # we want the number of operators to depend on the number of units
    unit.costing.operators_per_shift_var = Var(initialize=8, bounds=(0, 200))

    @unit.costing.Constraint()
    def operator_eqn(b):
        return b.operators_per_shift_var == 8 * unit.number_of_units

    unit.costing.annual_labor_cost_rule.deactivate()

    @unit.costing.Constraint()
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
    # unit.costing.other_variable_costs.unfix()

    # @unit.costing.Constraint()
    # def other_var_costs_eqn(b):
    #     return (
    #         b.other_variable_costs[0]
    #         == unit.number_of_units
    #         * dll_per_brick
    #         * unit.number_of_panels[0]
    #         * unit.bricks_per_panel
    #         * 1e-6
    #         / replacement_time
    #     )

    unit.costing.costing_initialization()


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


def print_dac_costing(tsa):
    fs = tsa.flowsheet()

    TPC_list = {}
    for o in fs.component_objects(descend_into=True):
        # look for costing blocks
        if hasattr(o, "costing") and hasattr(o.costing, "total_plant_cost"):
            for k in o.costing.total_plant_cost.keys():
                if k not in ["15.1", "15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k]
                if k in ["15.1"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] / 120 * tsa.number_beds
                if k in ["15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] * tsa.number_beds / 2

    for i, k in TPC_list.items():
        print(i, value(k))


def _var_dict_costing(tsa):

    # get flowsheet
    fs = tsa.flowsheet()

    # create dir with costing summary
    var_dict = {}

    var_dict["Annualized capital cost of dac unit [$MM/year]"] = value(
        fs.costing.annualized_cost
    )
    var_dict["Fixed O&M cost of dac unit [$MM/year]"] = value(
        fs.costing.total_fixed_OM_cost
    )
    var_dict["Variable O&M cost of dac unit [$MM/year]"] = value(
        fs.costing.total_variable_OM_cost[0]
    )
    var_dict["Total annualized cost of dac unit [$MM/year]"] = value(
        fs.costing.annualized_cost
        + fs.costing.total_fixed_OM_cost
        + fs.costing.total_variable_OM_cost[0] * fs.costing.capacity_factor
    )
    var_dict["Capture cost [$/tonne CO2]"] = value(fs.costing.cost_of_capture * 1e6)

    if hasattr(fs, "emissions_electric_boiler"):
        var_dict["Electric Boiler Emissions [mol/s]"] = value(
            fs.emissions_electric_boiler
        )

    if hasattr(fs, "emissions_electric_boiler_pv"):
        var_dict["Electric Boiler Emissions, PV electricity grid [mol/s]"] = value(
            fs.emissions_electric_boiler_pv
        )

    return var_dict


def dac_costing_summary(tsa, export=False):

    fs = tsa.flowsheet()

    if not hasattr(fs, "vessels"):
        raise ConfigurationError(f"{tsa.name} does not have any costing block.")

    var_dict = _var_dict_costing(tsa)

    summary_dir = {}
    summary_dir["Value"] = {}
    summary_dir["pos"] = {}

    count = 1
    for k, v in var_dict.items():
        summary_dir["Value"][k] = value(v)
        summary_dir["pos"][k] = count
        count += 1

    df = DataFrame.from_dict(summary_dir, orient="columns")
    del df["pos"]
    if export:
        df.to_csv(f"{tsa.local_name}_summary_costing.csv")

    print("\n" + "=" * 84)
    print(f"summary costing {tsa.local_name}")
    print("-" * 84)
    stdout.write(textwrap.indent(stream_table_dataframe_to_string(df), " " * 4))
    print("\n" + "=" * 84 + "\n")
