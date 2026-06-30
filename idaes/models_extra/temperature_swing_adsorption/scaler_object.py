from zmq import has

from idaes.core.scaling import CustomScalerBase, ConstraintScalingScheme
from pyomo.environ import Constraint


class TSA0DScaler(CustomScalerBase):
    """Scaler for the 0DTSA Model"""

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # call the scaling for each individual port/step
        self.inlet_port_variable_scaling_routine(model, overwrite, submodel_scalers)
        self.outlet_port_variable_scaling_routine(model, overwrite, submodel_scalers)
        self.heating_step_variable_scaling_routine(model, overwrite, submodel_scalers)
        self.cooling_step_variable_scaling_routine(model, overwrite, submodel_scalers)
        self.pressurization_step_variable_scaling_routine(
            model, overwrite, submodel_scalers
        )
        self.adsorption_step_variable_scaling_routine(
            model, overwrite, submodel_scalers
        )
        self.performance_variable_scaling_routine(model, overwrite, submodel_scalers)
        self.design_variable_scaling_routine(model, overwrite, submodel_scalers)

        if hasattr(model, "compressor"):
            self.call_submodel_scaler_method(
                submodel=model.compressor.unit,
                method="variable_scaling_routine",
                submodel_scalers=submodel_scalers,
                overwrite=overwrite,
            )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # call the scaling for each individual port/step
        self.inlet_port_constraint_scaling_routine(model, overwrite, submodel_scalers)
        self.outlet_port_constraint_scaling_routine(model, overwrite, submodel_scalers)
        self.heating_step_constraint_scaling_routine(model, overwrite, submodel_scalers)
        self.cooling_step_constraint_scaling_routine(model, overwrite, submodel_scalers)
        self.pressurization_step_constraint_scaling_routine(
            model, overwrite, submodel_scalers
        )
        self.adsorption_step_constraint_scaling_routine(
            model, overwrite, submodel_scalers
        )
        self.performance_constraint_scaling_routine(model, overwrite, submodel_scalers)
        self.design_constraint_scaling_routine(model, overwrite, submodel_scalers)

        if hasattr(model, "compressor"):
            self.call_submodel_scaler_method(
                submodel=model.compressor.unit,
                method="constraint_scaling_routine",
                submodel_scalers=submodel_scalers,
                overwrite=overwrite,
            )

        for c in model.component_data_objects(Constraint, descend_into=True):
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

    def inlet_port_variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Call scaling methods for variables at the inlet port
        if hasattr(model, "flow_mol_in"):
            for t in model.flowsheet().time:
                for c in model.component_list:
                    if c == "N2":
                        sf = 1e-2
                    elif c == "CO2":
                        sf = 1e-1
                    else:
                        sf = 1
                    self.set_variable_scaling_factor(
                        model.flow_mol_in[t, c], sf, overwrite
                    )
        if hasattr(model, "temperature_in"):
            for t in model.flowsheet().time:
                self.set_variable_scaling_factor(
                    model.temperature_in[t], 1e-2, overwrite
                )
        if hasattr(model, "pressure_in"):
            for t in model.flowsheet().time:
                self.set_variable_scaling_factor(model.pressure_in[t], 1e-5, overwrite)

    def inlet_port_constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale inlet port constraints
        if hasattr(model, "flow_mol_in_total_eq"):
            for c in model.flow_mol_in_total_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model, "mole_frac_in_eq"):
            for c in model.mole_frac_in_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseSum,  # best jac cond. number
                    overwrite=overwrite,
                )

        if hasattr(model, "pressure_in_eq"):
            for c in model.pressure_in_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

    def outlet_port_variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Call scaling methods for variables in the outlet port
        if hasattr(model, "flow_mol_co2_rich_stream"):
            for t in model.flowsheet().time:
                for c in model.isotherm_components:
                    if c == "CO2":
                        sf = 1e-1
                    elif c == "N2":
                        sf = 1
                    else:
                        sf = 1
                    self.set_variable_scaling_factor(
                        model.flow_mol_co2_rich_stream[t, c], sf, overwrite
                    )

        if hasattr(model, "temperature_co2_rich_stream"):
            for t in model.flowsheet().time:
                self.set_variable_scaling_factor(
                    model.temperature_co2_rich_stream[t], 1e-2, overwrite
                )

        if hasattr(model, "pressure_co2_rich_stream"):
            for t in model.flowsheet().time:
                self.set_variable_scaling_factor(
                    model.pressure_co2_rich_stream[t], 1e-5, overwrite
                )

        if hasattr(model, "flow_mol_n2_rich_stream"):  # TODO: default sf
            for t in model.flowsheet().time:
                for c in model.isotherm_components:
                    if c == "CO2":
                        sf = 1e-1
                    if c == "N2":
                        sf = 1e-2
                    else:
                        sf = 1
                    self.set_variable_scaling_factor(
                        model.flow_mol_n2_rich_stream[t, c], sf, overwrite
                    )

        if hasattr(model, "temperature_n2_rich_stream"):
            for t in model.flowsheet().time:
                self.set_variable_scaling_factor(
                    model.temperature_n2_rich_stream[t], 1e-2, overwrite
                )

        if hasattr(model, "pressure_n2_rich_stream"):
            for t in model.flowsheet().time:
                self.set_variable_scaling_factor(
                    model.pressure_n2_rich_stream[t], 1e-5, overwrite
                )

    def outlet_port_constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale outlet port constraints
        if hasattr(model, "flow_mol_co2_rich_stream_eq"):
            for c in model.flow_mol_co2_rich_stream_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseSum,  # best jacobian norm
                    overwrite=overwrite,
                )

        if hasattr(model, "temperature_co2_rich_stream_eq"):
            for c in model.temperature_co2_rich_stream_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model, "pressure_co2_rich_stream_eq"):
            for c in model.pressure_co2_rich_stream_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model, "flow_mol_n2_rich_stream_eq"):
            for c in model.flow_mol_n2_rich_stream_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model, "temperature_n2_rich_stream_eq"):
            for c in model.temperature_n2_rich_stream_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model, "pressure_n2_rich_stream_eq"):
            for c in model.pressure_n2_rich_stream_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        for c in model.flow_mol_h2o_o2_stream_eq.values():
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

        for c in model.temperature_h2o_o2_stream_eq.values():
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

        for c in model.pressure_h2o_o2_stream_eq.values():
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

    def heating_step_variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Call scaling methods for variables in the heating step
        if hasattr(model.heating, "time"):
            self.set_variable_scaling_factor(model.heating.time, 1e-2, overwrite)

        if hasattr(model.heating, "mole_frac"):  # update N2
            for t in model.heating.time_domain:
                for c in model.isotherm_components:
                    if c == "CO2":
                        sf = 1e1
                    if c == "N2":
                        sf = 1e3
                    else:
                        sf = 1
                    self.set_variable_scaling_factor(
                        model.heating.mole_frac[t, c], sf, overwrite
                    )

        if hasattr(model.heating, "temperature"):
            for t in model.heating.time_domain:
                self.set_variable_scaling_factor(
                    model.heating.temperature[t], 1e-2, overwrite
                )

        if hasattr(model.heating, "velocity_out"):
            for t in model.heating.time_domain:
                self.set_variable_scaling_factor(
                    model.heating.velocity_out[t], 1e3, overwrite
                )

        if hasattr(model.heating, "loading"):
            for t in model.heating.time_domain:
                self.set_variable_scaling_factor(
                    model.heating.loading[t, "CO2"], 1e1, overwrite
                )
                self.set_variable_scaling_factor(
                    model.heating.loading[t, "N2"], 1e10, overwrite
                )

        for t in model.heating.time_domain:
            for c in model.isotherm_components:
                if c == "CO2":
                    sf = 1e0
                else:
                    sf = 1e2
                self.set_variable_scaling_factor(
                    model.heating.mole_frac_dt[t, c], sf, overwrite
                )

    def heating_step_constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale heating step constraints
        if hasattr(model.heating, "component_mass_balance_ode"):
            for c in model.heating.component_mass_balance_ode.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.heating, "overall_mass_balance_ode"):
            for c in model.heating.overall_mass_balance_ode.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.heating, "energy_balance_ode"):
            for c in model.heating.energy_balance_ode.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.heating, "equil_loading_eq"):
            for c in model.heating.equil_loading_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.heating, "ic_mole_frac_eq"):
            for c in model.heating.ic_mole_frac_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.heating, "ic_temperature_eq"):
            for c in model.heating.ic_temperature_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.heating, "ic_velocity_out_eq"):
            for c in model.heating.ic_velocity_out_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.heating, "fc_temperature_eq"):
            for c in model.heating.fc_temperature_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

    def cooling_step_variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Call scaling methods for variables in the cooling step
        if hasattr(model.cooling, "time"):
            self.set_variable_scaling_factor(model.cooling.time, 1e-2, overwrite)

        if hasattr(model.cooling, "mole_frac"):
            for t in model.cooling.time_domain:
                self.set_variable_scaling_factor(
                    model.cooling.mole_frac[t, "CO2"], 1e1, overwrite
                )
                self.set_variable_scaling_factor(
                    model.cooling.mole_frac[t, "N2"], 1e6, overwrite
                )

        if hasattr(model.cooling, "temperature"):
            for t in model.cooling.time_domain:
                self.set_variable_scaling_factor(
                    model.cooling.temperature[t], 1e-2, overwrite
                )

        if hasattr(model.cooling, "pressure"):
            for t in model.cooling.time_domain:
                self.set_variable_scaling_factor(
                    model.cooling.pressure[t], 1e-4, overwrite
                )

        if hasattr(model.cooling, "loading"):
            for t in model.cooling.time_domain:
                self.set_variable_scaling_factor(
                    model.cooling.loading[t, "CO2"], 1e1, overwrite
                )
                self.set_variable_scaling_factor(
                    model.cooling.loading[t, "N2"], 1e10, overwrite
                )

        if hasattr(model.cooling, "mole_frac_heating_end"):
            self.set_variable_scaling_factor(
                model.cooling.mole_frac_heating_end, 1e1, overwrite
            )

        for t in model.cooling.time_domain:
            for c in model.isotherm_components:
                if c == "CO2":
                    sf = 1e0
                else:
                    sf = 1e5
                self.set_variable_scaling_factor(
                    model.cooling.mole_frac_dt[t, c], sf, overwrite
                )

        for t in model.cooling.time_domain:
            self.set_variable_scaling_factor(
                model.cooling.pressure_dt[t], 1e-5, overwrite
            )

    def cooling_step_constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale cooling step constraints
        if hasattr(model.cooling, "component_mass_balance_ode"):
            for c in model.cooling.component_mass_balance_ode.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.cooling, "overall_mass_balance_ode"):
            for c in model.cooling.overall_mass_balance_ode.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.cooling, "energy_balance_ode"):
            for c in model.cooling.energy_balance_ode.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.cooling, "equil_loading_eq"):
            for c in model.cooling.equil_loading_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.cooling, "ic_mole_frac_eq"):
            for c in model.cooling.ic_mole_frac_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.cooling, "ic_temperature_eq"):
            for c in model.cooling.ic_temperature_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.cooling, "sum_mole_frac"):
            for c in model.cooling.sum_mole_frac.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.cooling, "ic_pressure_eq"):
            for c in model.cooling.ic_pressure_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )
        if hasattr(model.cooling, "mole_frac_heating_end_eq"):
            for c in model.cooling.mole_frac_heating_end_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

    def pressurization_step_variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Call scaling methods for variables in the pressurization step
        if hasattr(model.pressurization, "time"):
            self.set_variable_scaling_factor(model.pressurization.time, 1, overwrite)

        if hasattr(model.pressurization, "mole_frac"):
            self.set_variable_scaling_factor(
                model.pressurization.mole_frac["CO2"], 1e2, overwrite
            )
            self.set_variable_scaling_factor(
                model.pressurization.mole_frac["N2"], 1e1, overwrite
            )

        if hasattr(model.pressurization, "loading"):
            self.set_variable_scaling_factor(
                model.pressurization.loading["CO2"], 1e1, overwrite
            )
            self.set_variable_scaling_factor(
                model.pressurization.loading["N2"], 1e10, overwrite
            )

        if hasattr(model.pressurization, "mole_frac_cooling_end"):
            self.set_variable_scaling_factor(
                model.pressurization.mole_frac_cooling_end, 1e1, overwrite
            )

        if hasattr(model.pressurization, "pressure_cooling_end"):
            self.set_variable_scaling_factor(
                model.pressurization.pressure_cooling_end, 1e-3, overwrite
            )

        if hasattr(model.pressurization, "loading_cooling_end"):
            self.set_variable_scaling_factor(
                model.pressurization.loading_cooling_end["CO2"], 1e1, overwrite
            )
            self.set_variable_scaling_factor(
                model.pressurization.loading_cooling_end["N2"], 1e10, overwrite
            )

    def pressurization_step_constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale pressurization step constraints

        if hasattr(model.pressurization, "sum_mole_frac_mass_balance_eq"):
            for c in model.pressurization.sum_mole_frac_mass_balance.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.pressurization, "equil_loading_eq"):
            for c in model.pressurization.equil_loading_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.pressurization, "pressurization_time_eq"):
            for c in model.pressurization.pressurization_time_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.pressurization, "mole_frac_cooling_end_eq"):
            for c in model.pressurization.mole_frac_cooling_end_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.pressurization, "pressure_cooling_end_eq"):
            for c in model.pressurization.pressure_cooling_end_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )
        if hasattr(model.pressurization, "loading_cooling_end_eq"):
            for c in model.pressurization.loading_cooling_end_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

    def adsorption_step_variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Call scaling methods for variables in the adsorption step
        if hasattr(model.adsorption, "time"):
            self.set_variable_scaling_factor(model.adsorption.time, 1e-2, overwrite)

        if hasattr(model.adsorption, "mole_frac_pressurization_end"):
            self.set_variable_scaling_factor(
                model.adsorption.mole_frac_pressurization_end, 1e2, overwrite
            )

        if hasattr(model.adsorption, "loading_pressurization_end"):
            self.set_variable_scaling_factor(
                model.adsorption.loading_pressurization_end["CO2"],
                1e1,
                overwrite,
            )
            self.set_variable_scaling_factor(
                model.adsorption.loading_pressurization_end["N2"],
                1e10,
                overwrite,
            )

        if hasattr(model.adsorption, "loading"):
            self.set_variable_scaling_factor(
                model.adsorption.loading["CO2"], 1, overwrite
            )
            self.set_variable_scaling_factor(
                model.adsorption.loading["N2"], 1e10, overwrite
            )

    def adsorption_step_constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Scale adsorption step constraints

        if hasattr(model.adsorption, "equil_loading_eq"):
            for c in model.adsorption.equil_loading_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.adsorption, "adsorption_time_eq"):
            for c in model.adsorption.adsorption_time_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.adsorption, "mole_frac_pressurization_end_eq"):
            for c in model.adsorption.mole_frac_pressurization_end_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        if hasattr(model.adsorption, "loading_pressurization_end_eq"):
            for c in model.adsorption.loading_pressurization_end_eq.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

    def performance_variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        # Call scaling methods for variables in the performance equations
        if hasattr(model, "mole_co2_in"):  # TODO: default sf
            self.set_variable_scaling_factor(model.mole_co2_in, 1, overwrite)

        if hasattr(model, "purity"):
            self.set_variable_scaling_factor(model.purity, 1e3, overwrite)

        if hasattr(model, "recovery"):
            self.set_variable_scaling_factor(model.recovery, 2e1, overwrite)

        if hasattr(model, "productivity"):
            self.set_variable_scaling_factor(model.productivity, 1e-1, overwrite)

        if hasattr(model, "cycle_time"):
            self.set_variable_scaling_factor(model.cycle_time, 1e1, overwrite)

        if hasattr(model, "thermal_energy"):  # TODO: default sf
            self.set_variable_scaling_factor(model.thermal_energy, 1, overwrite)

        if hasattr(model, "specific_energy"):
            self.set_variable_scaling_factor(model.specific_energy, 1, overwrite)

    def performance_constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        for c in model.mole_co2_in_eq.values():
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

        for c in model.purity_eq.values():
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseSum,  # gives better jacobian condition number
                overwrite=overwrite,
            )

        for c in model.recovery_eq.values():
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseRSS,
                overwrite=overwrite,
            )

        for c in model.cycle_time_eq.values():
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

        for c in model.productivity_eq.values():
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

        for c in model.thermal_energy_eq.values():
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

        for c in model.specific_thermal_energy_eq.values():
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

    def design_variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        self.set_variable_scaling_factor(
            model.flow_mol_in_total, 1e-2, overwrite
        )  # TODO:user SF
        self.set_variable_scaling_factor(
            model.mole_frac_in["CO2"], 1e2, overwrite
        )  # TODO:default SF
        self.set_variable_scaling_factor(
            model.mole_frac_in["N2"], 10, overwrite
        )  # TODO:default SF
        self.set_variable_scaling_factor(
            model.pressure_adsorption, 1e-4, overwrite
        )  # TODO:units SF
        self.set_variable_scaling_factor(
            model.temperature_adsorption, 1e-2, overwrite
        )  # TODO:default SF
        self.set_variable_scaling_factor(
            model.temperature_desorption, 1e-2, overwrite
        )  # TODO:default SF
        self.set_variable_scaling_factor(
            model.temperature_heating, 1e-2, overwrite
        )  # TODO:default SF
        self.set_variable_scaling_factor(
            model.temperature_cooling, 1e-2, overwrite
        )  # TODO:default SF
        self.set_variable_scaling_factor(
            model.bed_diameter, 10, overwrite
        )  # TODO:user SF
        self.set_variable_scaling_factor(model.bed_height, 1, overwrite)  # TODO:user SF
        self.set_variable_scaling_factor(
            model.pressure_drop, 1e-4, overwrite
        )  # TODO:default SF
        self.set_variable_scaling_factor(
            model.velocity_in, 10, overwrite
        )  # TODO:default SF
        self.set_variable_scaling_factor(
            model.velocity_mf, 10, overwrite
        )  # TODO:default SF
        for t in model.flowsheet().time:
            self.set_variable_scaling_factor(
                model.temperature_h2o_o2_stream[t], 1e-2, overwrite
            )  # TODO:default SF
            self.set_variable_scaling_factor(
                model.pressure_h2o_o2_stream[t], 1e-5, overwrite
            )  # TODO:default SF
            for k in ["H2O", "O2"]:
                self.set_variable_scaling_factor(
                    model.flow_mol_h2o_o2_stream[t, k], 1, overwrite
                )  # TODO:default SF, maybe?

    def design_constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        for c in model.velocity_mf_eq.values():
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )

        for c in model.pressure_drop_eq.values():
            self.scale_constraint_by_nominal_value(
                c,
                scheme=ConstraintScalingScheme.inverseMaximum,
                overwrite=overwrite,
            )
