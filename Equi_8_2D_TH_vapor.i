# changing vapor fixed
dbg_actions = false
dbg_residual_norms = true
first_ICs = true

BCintrusion = 'bottom'
BCtop = 'top'

[Mesh]
  # uniform_refine = 1
  type = GeneratedMesh
  dim = 2
  ny = 10
  nx = 15
  ymin = -7000.0
  ymax = 0
  xmax = 15000
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Debug]
  show_var_residual_norms = ${dbg_residual_norms}
  show_actions=${dbg_actions}
[]

[Functions]
  [pres_hydrostatic]
    type = ParsedFunction
    expression = '101325 - 1005.88*9.81*(y)'
  []
[]

[ICs]
  [pliquid]
    type = FunctionIC
    function = pres_hydrostatic
    variable = pliquid
    ignore_uo_dependency = ${first_ICs}
  []
  [h]
    # type = FunctionIC
    # function = 'energy_func'
    # variable = 'h'
    type = PorousFlowFluidPropertyIC
    variable = h
    property = enthalpy
    porepressure = pliquid
    temperature = 300
    temperature_unit = Kelvin
    fp = water97
    ignore_uo_dependency = ${first_ICs}
  []
[]

[BCs]
  [top_pressure]
    type = DirichletBC
    variable = pliquid
    value = 101325 # Atmospheric pressure Pa
    boundary = ${BCtop}
  []
  [top_temperature]
    # type = DirichletBC
    # variable = temperature
    # value = 288.15 #Kelvin 15C
    type = NeumannBC
    variable = h
    value = 109577 #288.15K
    boundary = ${BCtop}
  []
  [bottom_temperature_func]
    type = NeumannBC #
    # type = FunctionDirichletBC
    # function = 'intrusion_func'
    variable = h #temperature #
    value = 3e6 #2.888794e6
    boundary =  ${BCintrusion}
  []
[]

[AuxVariables]
  [pressure_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [pressure_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [enthalpy_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [enthalpy_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [saturation_gas]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 0
  []
  [saturation_water]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1
  []
  [density_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [density_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [viscosity_water]
    order = CONSTANT
    family = MONOMIAL
  []
  [viscosity_gas]
    order = CONSTANT
    family = MONOMIAL
  []
  [temperature]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [enthalpy_water]
    type = PorousFlowPropertyAux
    variable = enthalpy_water
    property = enthalpy
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [enthalpy_gas]
    type = PorousFlowPropertyAux
    variable = enthalpy_gas
    property = enthalpy
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [pressure_water]
    type = PorousFlowPropertyAux
    variable = pressure_water
    property = pressure
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [pressure_gas]
    type = PorousFlowPropertyAux
    variable = pressure_gas
    property = pressure
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [saturation_water]
    type = PorousFlowPropertyAux
    variable = saturation_water
    property = saturation
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [saturation_gas]
    type = PorousFlowPropertyAux
    variable = saturation_gas
    property = saturation
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [density_water]
    type = PorousFlowPropertyAux
    variable = density_water
    property = density
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [density_gas]
    type = PorousFlowPropertyAux
    variable = density_gas
    property = density
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [viscosity_water]
    type = PorousFlowPropertyAux
    variable = viscosity_water
    property = viscosity
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [viscosity_gas]
    type = PorousFlowPropertyAux
    variable = viscosity_gas
    property = viscosity
    phase = 1
    execute_on = 'initial timestep_end'
  []
  [temperature]
    type = PorousFlowPropertyAux
    variable = temperature
    property = temperature
    execute_on = 'initial timestep_end'
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pliquid h'
    number_fluid_phases = 2
    number_fluid_components = 1
  []
  [pc] #from porousflow/test/tests/fluidstate/watervapor.i
    type = PorousFlowCapillaryPressureBC
    pe = 1e5
    lambda = 2
    pc_max = 1e6
  []
  [fs]
    type = PorousFlowWaterVapor
    water_fp = water97
    # moose say 'temperature_from_ph() not implemented for region 5'
    #if I use tabulated values calculated with Coolprop, water_tab (line 228), moose says that 'triplePointPressure()' const not implemented-
    capillary_pressure = pc
  []
[]

[Variables]
  [pliquid]
    # initial_condition = 2.2e7 # 1e6
  []
  [h]
    # initial_condition = 2e6 #440e3
    scaling = 1e-8
  []
[]

[Kernels]
  [mass]
    type = PorousFlowMassTimeDerivative
    variable = pliquid
  []
  [heat]
    type = PorousFlowEnergyTimeDerivative
    variable = h
  []
[]

[FluidProperties]
  [water97]
    type = Water97FluidProperties # IAPWS-IF97
  []
  [./water_tab]
    type = TabulatedBicubicFluidProperties
    fp = water97
    # p_h_variables = true
    fluid_property_file = water_fluid_properties_Coolprop.csv
    # interpolated_properties = 'pressure temperature density enthalpy viscosity c'
    interpolated_properties = 'density enthalpy viscosity c'
    temperature_min = 273
    temperature_max = 1400
    pressure_min = 15000000
    pressure_max = 401000000
#    num_T = 1000
#    num_p = 1000
  [../]
[]

[Materials]
  [watervapor]
    type = PorousFlowFluidStateSingleComponent
    porepressure = pliquid
    enthalpy = h
    capillary_pressure = pc
    fluid_state = fs
    temperature_unit = Kelvin
  []
  [porosity_wallrock]
    type = PorousFlowPorosityConst
    porosity = 0.1
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-13 0 0   0 1E-13 0   0 0 1E-13'
  []
  [relperm_water] # from watervapor.i
    type = PorousFlowRelativePermeabilityCorey
    n = 2
    phase = 0
  []
  [relperm_gas] # from watervapor.i
    type = PorousFlowRelativePermeabilityCorey
    n = 3
    phase = 1
  []
  [internal_energy]
    type = PorousFlowMatrixInternalEnergy
    density = 2500 # kg/m^3
    specific_heat_capacity = 1200 # J/kg/K
  []
[]

# [Executioner]
#   type = Transient
#   solve_type = NEWTON
#   end_time = 1
#   nl_abs_tol = 1e-12
# []

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = -3.15576e+10 #1000 years for stabilization of pliquid
  end_time = 9.467279e+11 # 3.15576e+12 #100k years to recover after switch off the intrusion
  dtmax = 3.15576e+10 #1000y #473364000 # 15y
  dtmin = 1750.0

  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
  nl_max_its = 400 # max Non linear iterations before cutback is applied
  l_tol = 1e-12
  l_abs_tol = 1e-12
  l_max_its = 1000 #1000

  [TimeStepper]
    type = IterationAdaptiveDT
    # optimal_iterations = 5
    dt = 31557600 # 1y
    growth_factor = 2 # if iterations is less than nl max --> if iteration runs smoothly
    cutback_factor = 0.5 # if iteration exceeds nl max --> runs poorly
  []
[]

[Dampers]
  [./limit]
    type = BoundingValueNodalDamper
    variable = pliquid
    max_value = 4e8
    min_value = 1e1
  [../]
[]


[Postprocessors]
  [density_water]
    type = ElementAverageValue
    variable = density_water
    execute_on = 'initial timestep_end'
  []
  [density_gas]
    type = ElementAverageValue
    variable = density_gas
    execute_on = 'initial timestep_end'
  []
  [viscosity_water]
    type = ElementAverageValue
    variable = viscosity_water
    execute_on = 'initial timestep_end'
  []
  [viscosity_gas]
    type = ElementAverageValue
    variable = viscosity_gas
    execute_on = 'initial timestep_end'
  []
  [enthalpy_water]
    type = ElementAverageValue
    variable = enthalpy_water
    execute_on = 'initial timestep_end'
  []
  [enthalpy_gas]
    type = ElementAverageValue
    variable = enthalpy_gas
    execute_on = 'initial timestep_end'
  []
  [sg]
    type = ElementAverageValue
    variable = saturation_gas
    execute_on = 'initial timestep_end'
  []
  [sw]
    type = ElementAverageValue
    variable = saturation_water
    execute_on = 'initial timestep_end'
  []
  [pwater]
    type = ElementAverageValue
    variable = pressure_water
    execute_on = 'initial timestep_end'
  []
  [pgas]
    type = ElementAverageValue
    variable = pressure_gas
    execute_on = 'initial timestep_end'
  []
  [temperature]
    type = ElementAverageValue
    variable = temperature
    execute_on = 'initial timestep_end'
  []
  [enthalpy]
    type = ElementAverageValue
    variable = h
    execute_on = 'initial timestep_end'
  []
  [liquid_mass]
    type = PorousFlowFluidMass
    phase = 0
    execute_on = 'initial timestep_end'
  []
  [vapor_mass]
    type = PorousFlowFluidMass
    phase = 1
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  exodus = true
  # file_base = water_vapor_twophase_tab
  # csv = true
[]

[Controls]
  [./intrusion_on]
      type = TimePeriod
      enable_objects = 'BCs::bottom_temperature_func'
      start_time = '0.0'
      end_time = '9.467279e+11' #test 30ky active
      execute_on ='INITIAL TIMESTEP_END'
  [../]
[]

# ###If you uncomment this it will print out all the kernels and materials that the PorousFlowFullySaturated action generates
# [Problem]
#   type = DumpObjectsProblem
#   dump_path = Equi_8_2D_THM_vapor_PorousFlowMaterials
#  []