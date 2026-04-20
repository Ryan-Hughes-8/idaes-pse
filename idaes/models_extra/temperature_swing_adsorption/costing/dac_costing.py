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
from idaes.core.util.math import smooth_max

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


def _nonneg(expr, unit):
    """
    Return a smooth nonnegative version of expr with the given units.
    """
    return smooth_max(0, units.convert(expr, to_units=unit) / unit) * unit


def _warn_if_negative(expr, name):
    try:
        v = value(expr)
    except Exception:
        return
    if v < 0:
        _log.warning(
            "%s is negative (%.6g). Clamping to zero for costing stability.",
            name,
            v,
        )


def get_dac_costing(unit, costing_case):
    # capital costs
    costing_params = get_dac_costing_data(costing_case)

    # get flowsheet
    fs = unit.flowsheet()

    fs.costing = QGESSCosting()

    # reference parameters for accounts ==============================
    # === accounts present for both costing cases ===
    # compressor auxiliary load - from surrogates
    # TODO: this compressor power is just from surrogates, for vacuum support, need to model with compressor unit model
    _pcal_dimless = (
        0.0012
        * units.convert(unit.flow_mol_in_total, to_units=units.kmol / units.hr)
        * units.hr
        / units.kmol
        - 2.2798
    )
    _warn_if_negative(_pcal_dimless, "product_compressor_auxiliary_load surrogate")
    product_compressor_auxiliary_load = smooth_max(0, _pcal_dimless) * units.kW  # [kW]
    # compressor aftercooler heat exchanger duty - from surrogates
    _cahd_dimless = (
        2e-6
        * units.convert(unit.flow_mol_in_total, to_units=units.kmol / units.hr)
        * units.hr
        / units.kmol
        - 7e-8
    )
    _warn_if_negative(_cahd_dimless, "compressor_aftercooler_heat_duty surrogate")
    compressor_aftercooler_heat_duty = (
        smooth_max(0, _cahd_dimless) * units.MBtu / units.hr
    )  # [MMBtu/hr]

    auxiliary_load_2_beds = _nonneg(
        2
        * units.convert(unit.compressor.unit.work_mechanical[0], to_units=units.kW)
        / unit.number_beds,
        units.kW,
    )  # [kW]

    # CO2 product mass flow rate - from model
    CO2_product_mass_flow = units.convert(
        unit.mw["CO2"] * unit.flow_mol_co2_rich_stream[0, "CO2"],
        to_units=units.lb / units.hr,
    )  # [lb/hr]

    if costing_case == "electric_boiler":
        # raw water withdrawal flow rate - from surrogates
        # TODO: check which expression is correct (CDR has this as the correlation for ngcc)
        raw_water_withdrawal = (
            (1.0496e-03 * CO2_product_mass_flow * units.hr / units.lb + 5.3359e01)
            * units.gal
            / units.min
        )  # [gpm]
        # process water discharge flow rate - from surrogates
        # TODO: check which expression is correct (CDR has this as the correlation for ngcc)
        process_water_discharge = (
            (5.4007e-04 * CO2_product_mass_flow * units.hr / units.lb + 1.2000e01)
            * units.gal
            / units.min
        )  # [gpm]
        # cooling tower heat duty - from surrogates
        # TODO: check which expression is correct (CDR has this as the correlation for ngcc)
        cooling_tower_duty = (
            (3.0797e-04 * CO2_product_mass_flow * units.hr / units.lb + 2.5000e01)
            * units.MBtu
            / units.hr
        )  # [MMBtu/hr]
        # circulating water flow_rate - from surrogates
        # TODO: replace with 100*cooling_tower_duty*units_adjustment
        circulating_water_flow_rate = (
            (3.0797e-02 * CO2_product_mass_flow * units.hr / units.lb + 2.5000e03)
            * units.gal
            / units.min
        )  # [gpm]
        # boiler auxiliary load - from surrogates
        boiler_auxiliary_load = _nonneg(
            0.000335999
            * units.convert(unit.flow_mass_steam, to_units=units.lb / units.hr)
            * units.hr
            / units.lb
            * 1e3
            * units.kW,
            units.kW,
        )  # convert MW to kW
        # total auxiliary load
        # calculate with fans work + compressor for CO2 pure + boiler aux load
        total_auxiliary_load = _nonneg(
            product_compressor_auxiliary_load
            + boiler_auxiliary_load
            + units.convert(unit.compressor.unit.work_mechanical[0], to_units=units.kW),
            units.kW,
        )

        unit.raw_water_system = UnitModelBlock()
        unit.raw_water_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["3.2", "3.4", "9.5", "14.6"],
                "scaled_param": raw_water_withdrawal,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        unit.steam_system = UnitModelBlock()
        unit.steam_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["3.1", "3.3", "3.5"],
                "scaled_param": units.convert(
                    unit.flow_mass_steam, to_units=units.lb / units.hr
                ),
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        unit.cooling_tower = UnitModelBlock()
        unit.cooling_tower.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["9.1"],
                "scaled_param": cooling_tower_duty,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        unit.water_discharge_system = UnitModelBlock()
        unit.water_discharge_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["3.7"],
                "scaled_param": process_water_discharge,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        unit.cooling_water_system = UnitModelBlock()
        unit.cooling_water_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["9.2", "9.3", "9.4", "9.6", "9.7", "14.5"],
                "scaled_param": circulating_water_flow_rate,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        unit.electric_systems = UnitModelBlock()
        unit.electric_systems.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
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
                "scaled_param": total_auxiliary_load,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        # Electric Boiler 15.9
        unit.electric_boiler = UnitModelBlock()
        unit.electric_boiler.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["15.9"],
                "scaled_param": units.convert(
                    unit.flow_mass_steam, to_units=units.lb / units.hr
                ),
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

    elif costing_case == "retrofit_ngcc":
        # calculate with fans work + compressor for CO2 pure
        total_auxiliary_load = _nonneg(
            product_compressor_auxiliary_load
            + units.convert(unit.compressor.unit.work_mechanical[0], to_units=units.kW),
            units.kW,
        )

        unit.steam_flow_system = UnitModelBlock()
        unit.steam_flow_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["8.4"],
                "scaled_param": units.convert(
                    unit.flow_mass_steam, to_units=units.lb / units.hr
                ),
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        unit.electric_systems = UnitModelBlock()
        unit.electric_systems.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["11.2", "11.3", "11.4", "11.5", "11.6"],
                "scaled_param": total_auxiliary_load,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )
    else:
        print("Invalid Case")

    fs.vessels = UnitModelBlock()
    fs.vessels.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.1"],
            "scaled_param": units.convert(unit.bed_volume, to_units=units.ft**3),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # purge_pressure = unit.config.purge_pressure
    # 15.2 - DAC CO2 Compression & Drying
    fs.product_compression = UnitModelBlock()
    # if purge_pressure == 1:
    fs.product_compression.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.2"],
            "scaled_param": product_compressor_auxiliary_load,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # elif purge_pressure == 0.2:
    #     unit.product_compression.costing = UnitModelCostingBlock(
    #         flowsheet_costing_block=fs.costing,
    #         costing_method=QGESSCostingData.get_PP_costing,
    #         costing_method_arguments={
    #             "cost_accounts": ["vac_2.9_eq1"],
    #             "scaled_param": unit.compressor_power[0],
    #             "tech": 9,
    #             "ccs": "B",
    #             "additional_costing_params": vacuum_pump_params,
    #         },
    #     )

    # elif purge_pressure == 0.5:
    #     unit.product_compression.costing = UnitModelCostingBlock(
    #         flowsheet_costing_block=fs.costing,
    #         costing_method=QGESSCostingData.get_PP_costing,
    #         costing_method_arguments={
    #             "cost_accounts": ["vac_6.5_eq1"],
    #             "scaled_param": unit.compressor_power[0],
    #             "units": "MW",
    #             "tech": 9,
    #             "ccs": "B",
    #             "additional_costing_params": vacuum_pump_params,
    #         },
    #     )

    # 15.3 - DAC CO2 Compressor Aftercooler
    fs.compressor_aftercooler = UnitModelBlock()
    fs.compressor_aftercooler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.3"],
            "scaled_param": compressor_aftercooler_heat_duty,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.4 - DAC System Air Handling Duct and Dampers (1 system per 2 beds)
    # TODO: double check the 1 unit per two beds assumption here
    fs.duct_dampers = UnitModelBlock()
    fs.duct_dampers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.4"],
            "scaled_param": 2
            * units.convert(unit.flow_mass_in_total_bed, to_units=units.lb / units.hr),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.5 - DAC System Air Handling Fans (1 system per 2 beds)
    fs.feed_fans = UnitModelBlock()
    fs.feed_fans.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.5"],
            "scaled_param": auxiliary_load_2_beds,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.6 - DAC Desorption Process Gas Handling System (pure CO2 gas)
    fs.desorption_gas_handling = UnitModelBlock()
    fs.desorption_gas_handling.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.6"],
            "scaled_param": CO2_product_mass_flow,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.7 - DAC Steam Distribution System
    fs.steam_distribution = UnitModelBlock()
    fs.steam_distribution.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.7"],
            "scaled_param": units.convert(
                unit.flow_mass_steam, to_units=units.lb / units.hr
            ),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.8 - DAC System Controls Equipment
    fs.controls_equipment = UnitModelBlock()
    fs.controls_equipment.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.8"],
            "scaled_param": total_auxiliary_load,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # we need a custom method of calculating the total TPC to account for distributed vs centralized systems
    distributed_systems = [
        fs.vessels.costing,
        fs.duct_dampers.costing,
        fs.feed_fans.costing,
        fs.desorption_gas_handling.costing,
        fs.controls_equipment.costing,
    ]

    centralized_TPCs = []

    for b in fs.costing._registered_unit_costing:
        if b not in distributed_systems:
            for key in b.total_plant_cost.keys():
                centralized_TPCs.append(b.total_plant_cost[key])

    fs.costing.total_TPC = Var(
        initialize=100,
        bounds=(0, 1e4),
        doc="total TPC in $MM",
    )

    @fs.costing.Constraint()
    def total_TPC_eq(b):
        return b.total_TPC == sum(centralized_TPCs) + (
            fs.vessels.costing.total_plant_cost["15.1"] / 120 * unit.number_beds
            + fs.duct_dampers.costing.total_plant_cost["15.4"] * unit.number_beds / 2
            + fs.feed_fans.costing.total_plant_cost["15.5"] * unit.number_beds / 2
            + fs.desorption_gas_handling.costing.total_plant_cost["15.6"]
            + fs.controls_equipment.costing.total_plant_cost["15.8"]
        )

    # resource consumption rates for variable costing
    capacity_factor = 0.85
    sorbent_lifespan = 0.5

    @fs.costing.Expression(fs.time)
    def water_use(b, t):
        # water use is lineraly scaled from NETL reference
        ref_water = 60792.407 * units.gal / units.day
        ref_air = 1629629 * units.kmol / units.hr
        return units.convert(unit.flow_mol_in_total, to_units=units.kmol / units.hr) * (
            ref_water / ref_air
        )

    @fs.costing.Expression(fs.time)
    def water_treatment_chems(b, t):
        # treatment chemical use is lineraly scaled from NETL reference
        ref_chem = 0.1811 * units.ton / units.day
        ref_air = 1629629 * units.kmol / units.hr
        return units.convert(unit.flow_mol_in_total, to_units=units.kmol / units.hr) * (
            ref_chem / ref_air
        )

    @fs.costing.Expression(fs.time)
    def energy_purchased(b, t):  # in kWh/day
        hr_per_day = 24 * units.hr / units.day
        return total_auxiliary_load * hr_per_day

    @fs.costing.Expression(fs.time)
    def sorbent_rate(b, t):
        return (
            units.convert(unit.bed_volume, to_units=units.ft**3)
            * (1 - unit.bed_voidage)
            * unit.number_beds
            / sorbent_lifespan
            / 365
            / units.day
        )

    # TODO: add this back in
    # @fs.costing.Expression(fs.time)
    # def bfw_rate(b, t):
    #     hr_per_day = 24 * units.hr / units.day
    #     return unit.BFW_makeup[t] * hr_per_day

    @fs.costing.Expression(fs.time)
    def steam_rate(b, t):
        return units.convert(unit.flow_mass_steam, to_units=units.kg / units.day)

    fs.costing.net_power = Var(fs.time, initialize=690, units=units.MW)
    fs.costing.net_power.fix()

    # resorces to be costed
    resources = [
        "water",
        "water_treatment_chemicals",
        "aux_power",
        "sorbent",
        "waste_sorbent",
        # "boiler_feed_water", #TODO: add this back in
    ]
    if costing_case == "retrofit_NGCC":
        resources.append("IP_steam")

    # vars for resource consumption rates
    rates = [
        fs.costing.water_use,
        fs.costing.water_treatment_chems,
        fs.costing.energy_purchased,
        fs.costing.sorbent_rate,
        fs.costing.sorbent_rate,
        # fs.costing.bfw_rate, #TODO: add this back in
    ]
    if costing_case == "retrofit_NGCC":
        rates.append(fs.costing.steam_rate)

    # resource prices
    prices = {
        # "sorbent": 4 * units.USD_2018 / units.ft**3,  # 4 or 100
        "sorbent": 201 * units.USD_2018 / units.ft**3,
        "aux_power": 0.06 * units.USD_2018 / units.kWh,
        "waste_sorbent": 0.86 * units.USD_2018 / units.ft**3,
        "IP_steam": 0.00733 * units.USD_2018 / units.kg,
        # "boiler_feed_water": 2.45
        # / 1000
        # * units.USD_2018
        # / units.kg,  # Turton et al., 2012
    }

    @fs.costing.Expression()
    def land_cost_exp(b):
        return (
            156000 * (unit.number_beds / 120) ** (0.78)
        ) * 1e-6  # scaled to Millions

    fs.costing.build_process_costs(
        # net_power=None,
        net_power=fs.costing.net_power,
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
        land_cost=fs.costing.land_cost_exp,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel=None,
        # waste=None,
        # chemicals=None,
        tonne_CO2_capture=unit.total_CO2_captured_year,
    )

    fs.costing.costing_initialization()


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


def print_dac_costing(unit):
    fs = unit.flowsheet()

    TPC_list = {}
    for o in fs.component_objects(descend_into=True):
        # look for costing blocks
        if hasattr(o, "costing") and hasattr(o.costing, "total_plant_cost"):
            for k in o.costing.total_plant_cost.keys():
                if k not in ["15.1", "15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k]
                if k in ["15.1"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] / 120 * unit.number_beds
                if k in ["15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] * unit.number_beds / 2

    for i, k in TPC_list.items():
        print(i, value(k))


def _var_dict_costing(unit):

    # get flowsheet
    fs = unit.flowsheet()

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


def dac_costing_summary(unit, export=False):

    fs = unit.flowsheet()

    if not hasattr(fs, "vessels"):
        raise ConfigurationError(f"{unit.name} does not have any costing block.")

    var_dict = _var_dict_costing(unit)

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
        df.to_csv(f"{unit.local_name}_summary_costing.csv")

    print("\n" + "=" * 84)
    print(f"summary costing {unit.local_name}")
    print("-" * 84)
    stdout.write(textwrap.indent(stream_table_dataframe_to_string(df), " " * 4))
    print("\n" + "=" * 84 + "\n")
