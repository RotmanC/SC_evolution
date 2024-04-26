### Rotman A. Criollo Manjarrez
### CO2SINk
### equilibration with k=f(T)
###
#
# PARAMETERS USED
automatic_scaling = true
dbg_actions = false
dbg_residual_norms = true

[Mesh]
  [file]
    type = FileMeshGenerator
    file = '2D_15km_centrado.msh'
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  displacements = 'disp_x disp_y'
  gravity = '0 -9.81 0'
[]

[Variables]
  [pres_water]
    # scaling = 1E-5
  []
  [sgas]
  []
  [temp]
    # scaling = 1E5
  []
[]

[ICs]
  [pres_water]
    type = FunctionIC
    function = pres_wateric
    variable = pres_water
  []
  [temp]
    type = FunctionIC
    function = tempic
    variable = temp
  []
  [disp_y]
    type = FunctionIC
    function = '0.5 * y'
    variable = disp_y
  []
[]

[Functions]
  [pres_wateric]
    type = ParsedFunction
    expression = '101325 - 1005.88*9.81*(y)'
  []
  [tempic]
    type = ParsedFunction
    expression = '288.15 - (0.05*y)' #Geograd 50C/km
  []
  [volcan_func]
    type = ParsedFunction
    symbol_names = 'start_time Tini Tvolcan'
    symbol_values = '86400.0 638.15 823.15'
    expression = if(t<start_time,Tini,Tvolcan)
  []
[]

[BCs]
  # [top_pressure_sink]
  #   type = PorousFlowSink
  #   boundary = 'BCtop'
  #   variable = pres_water
  #   use_mobility = false
  #   use_relperm = false
  #   fluid_phase = 0
  #   flux_function = -101325 # Atmospheric pressure Pa. Negative means a source, rather than a sink
  # []
  # [top_temperature_sink]
  #   type = PorousFlowSink
  #   boundary = 'BCtop'
  #   variable = temp
  #   use_mobility = false
  #   use_relperm = false
  #   fluid_phase = 0
  #   flux_function = -288.15 #Kelvin 15C. Negative means a source, rather than a sink
  # []
  [top_pressure]
    type = DirichletBC
    variable = pres_water
    value = 101325 # Atmospheric pressure Pa
    boundary = 'BCtop'
  []
  [top_temperature]
    type = DirichletBC
    variable = temp
    value = 288.15 #Kelvin 15C
    boundary = 'BCtop'
  []
  [bottom_temperature]
    type = FunctionDirichletBC #DirichletBC #
    variable = temp
    function = volcan_func
    # value = 823.15 #Kelvin, "550C"
    boundary = 'BCvolcan'
  []
[]

[AuxVariables]
  [massfrac_ph0_sp0]
    initial_condition = 1
  []
  [massfrac_ph1_sp0]
    initial_condition = 0
  []
  [massfrac_ph0_sp1]
  []
  [massfrac_ph1_sp1]
  []
  [pgas]
    family = MONOMIAL
    order = FIRST
  []
  [swater]
    family = MONOMIAL
    order = FIRST
  []
  [./density_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./density_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./viscosity_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./viscosity_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./enthalpy_gas]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./enthalpy_water]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [permeability]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [mass_water_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pres_water
  []
  [flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pres_water
  []
  [mass_co2_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = sgas
  []
  [flux_co2]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    variable = sgas
  []
  [energy_dot]
    type = PorousFlowEnergyTimeDerivative
    variable = temp
  []
  [advection]
    type = PorousFlowHeatAdvection
    variable = temp
  []
  [conduction]
    type = PorousFlowHeatConduction
    variable = temp
  []
[]

[AuxKernels]
  [pgas]
    type = PorousFlowPropertyAux
    property = pressure
    phase = 1
    variable = pgas
  []
  [swater]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = swater
  []
  [./enthalpy_gas]
    type = PorousFlowPropertyAux
    property = enthalpy
    phase = 1
    variable = enthalpy_gas
  [../]
  [./enthalpy_water]
    type = PorousFlowPropertyAux
    property = enthalpy
    phase = 0
    variable = enthalpy_water
  [../]
  [./density_water]
    type = PorousFlowPropertyAux
    property = density
    phase = 0
    variable = density_water
  [../]
  [./density_gas]
    type = PorousFlowPropertyAux
    variable = density_gas
    phase = 1
    property = density
  [../]
  [./viscosity_water]
    type = PorousFlowPropertyAux
    property = viscosity
    phase = 0
    variable = viscosity_water
  [../]
  [./viscosity_gas]
    type = PorousFlowPropertyAux
    property = viscosity
    phase = 1
    variable = viscosity_gas
  [../]
  [permeability]
    type = ParsedAux
    # scaling = 1E21
    coupled_variables = 'temp'
    variable = 'permeability' #for the vectorpostprocessor
    expression = 'pow(10,((22-13.5)/(1+exp(0.28*(temp-652.23)))-22))'
    # expression = perm_func
    execute_on = 'INITIAL' #TIMESTEP_END'
  []
[]

[UserObjects]
  [produced_mass_h2o]
    type = PorousFlowSumQuantity
  []
  [produced_mass_co2]
    type = PorousFlowSumQuantity
  []
  [produced_heat_h2o]
    type = PorousFlowSumQuantity
  []
  [produced_heat_co2]
    type = PorousFlowSumQuantity
  []
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'temp pres_water sgas'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
    # type = PorousFlowCapillaryPressureConst
    # pc = 0
   type = PorousFlowCapillaryPressureVG
   alpha = 1e-5
   m = 0.5
  []
[]


[FluidProperties]
  [./water]
    type = Water97FluidProperties
  [../]
  [./tabulated_water]
    type = TabulatedBicubicFluidProperties
    fp = water
    fluid_property_file = water_fluid_properties_Coolprop.csv
    interpolated_properties = 'density enthalpy viscosity'
    temperature_min = 273
    temperature_max = 1400
    pressure_min = 15000000
    pressure_max = 401000000
#    num_T = 1000
#    num_p = 1000
  [../]
  [./co2]
    type = CO2FluidProperties
  [../]
  [./tabulated_co2]
    type = TabulatedBicubicFluidProperties
    fp = co2
    fluid_property_file = co2_fluid_properties_Coolprop.csv
    interpolated_properties = 'density enthalpy viscosity'
    temperature_min = 273
    temperature_max = 1400
    pressure_min = 15000000
    pressure_max = 401000000
#    num_T = 1000
#    num_p = 1000
  [../]
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = temp
  []
  [ppss]
    type = PorousFlow2PhasePS #PorousFlow1PhaseP
    # porepressure = pres_water
    phase0_porepressure = pres_water
    phase1_saturation = sgas
    capillary_pressure = pc
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'massfrac_ph0_sp0 massfrac_ph1_sp0'
  []
  [water]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
  []
  [gas]
    type = PorousFlowSingleComponentFluid
    fp = co2
    phase = 1
  []
  [relperm_liquid]
    type = PorousFlowRelativePermeabilityVG
    phase = 0
    m = 0.5
    s_res = 0.2
    sum_s_res = 0.5
  []
  [relperm_gas]
    type = PorousFlowRelativePermeabilityVG
    phase = 1
    m = 0.5
    s_res = 0.2
    sum_s_res = 0.5
  []
  [porosity]
    type = PorousFlowPorosityConst # only the initial value of this is ever used
    porosity = 0.1
  []
  [biot_modulus]
    type = PorousFlowConstantBiotModulus
    solid_bulk_compliance = 1E-10
    fluid_bulk_modulus = 2E9
  []
  [permeability]
    # type = PorousFlowPermeabilityConst
    # permeability = '1e-15 0 0   0 1e-15 0   0 0 1e-15'
    type = PorousFlowPermeabilityTensorFromVar
    perm = permeability  
  []
  [thermal_expansion]
    type = PorousFlowConstantThermalExpansionCoefficient
    fluid_coefficient = 5E-6
    drained_coefficient = 2E-4
  []
  [thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '1 0 0  0 1. 0  0 0 1'
    wet_thermal_conductivity = '3 0 0  0 3. 0  0 0 3'
  []
  [rock_heat]
    type = PorousFlowMatrixInternalEnergy
    density = 2500.0
    specific_heat_capacity = 1200.0
  []
[]

[VectorPostprocessors]
  [parameters]
    type = LineValueSampler
    start_point = '7500.0 -7000.0 0.0'
    end_point = '7500.0 00.0 0.0'
    num_points = 700
    sort_by = z
    variable = 'pres_water temp permeability'
  []
[]

[Debug]
  show_var_residual_norms = ${dbg_residual_norms}
  show_actions=${dbg_actions}
[]

[Preconditioning]
  active = 'smp'
  [smp]
    type = SMP
    full = true
    # petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'gmres      lu      asm           NONZERO                   2               1E1       1E-3        500'
  []
  [mumps]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -snes_rtol -snes_atol -snes_max_it'
    petsc_options_value = 'gmres      lu       mumps                         NONZERO               1E-3       1E1       500'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK' #NEWTON
  automatic_scaling = ${automatic_scaling}
  end_time = 1104516000000 #35k years # 946728000000 # 30000y
  nl_max_its = 400 # max Non linear iterations before cutback is applied
  l_max_its = 1000
  dtmax = 473364000 # 15y
  dtmin = 1750.0
  nl_abs_tol = 1e-3
#  steady_state_detection = true
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 31557600 # 1y
    growth_factor = 2 # if iterations is less than nl max --> if iteration runs smoothly
    cutback_factor = 0.5 # if iteration exceeds nl max --> runs poorly
  []
[]

[Outputs]
  perf_graph = true
  exodus = true
 [csv]
   type = CSV
   file_base = 2D_equi1_
   execute_on = 'TIMESTEP_END'
 []
[]
