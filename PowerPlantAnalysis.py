'''
Power Plant Design Creation and Analysis Tool
Created by Elvis Imamura
ME 123

This program generates information for power cycles using thermophysical properties from the PyroMat python library. The user is able to specify a temperature and pressure range to test, and the program will provide the optimal cycles for the parameter specified.
'''

import pyromat as pm
import numpy as np
import pandas as pd

class State:
	'''
	Creates a State object defined by a given pressure and temperature, enabling calculation of all other thermodynamical values
	'''	
	water = pm.get('mp.H2O')
	
	def __init__(self, p, T):
		self.p = p
		self.T = T
		self.fluid = pm.get('mp.H2O')

	def __str__(self):
		print("Information for state:")
		print("p :", self.get_p())
		print("T :", self.get_T())
		print("h :", self.get_h())
		print("s :", self.get_s())
		print("v :", self.get_v())
		return ""
	
	def __repr__(self):
		water = self.fluid
		print("p :", self.get_p())
		print("T :", self.get_T())
		print("h :", self.get_h())
		print("s :", self.get_s())
		print("v :", self.get_v())
		return ""

	def get_p(self):
		return float(self.p)
	
	def get_T(self):
		return float(self.T)

	def get_h(self):
		water = self.fluid
		return float(water.h(p = self.p, T = self.T))
	
	def get_s(self):
		water = self.fluid
		return float(water.s(p = self.p, T = self.T))

	def get_v(self):
		water = self.fluid
		return float(water.v(p = self.p, T = self.T))

class MixState(State):
	'''
	Creates a State object, with the option to specify a certain quality value (used exclusively for mixtures/saturated states
	'''
	def __init__(self, p, T, x):
		self.p = p
		self.T = T
		self.x = x
		self.fluid = pm.get('mp.H2O')

	def get_h(self):
		water = self.fluid
		return float(water.h(p = self.p, x = self.x))
	
	def get_s(self):
		water = self.fluid
		return float(water.s(p = self.p, x = self.x))

	def get_v(self):
		water = self.fluid
		return float(water.v(p = self.p, x = self.x))

	def get_x(self):
		return self.x

class Cycle:
	'''
	Creates a Cycle object consisting of 4 states (one turbine); calculates all information including power input/output, and thermal efficiency
	'''
	water = pm.get('mp.H2O')
	
	def __init__(self, s1, s2, s3, s4, massflow):
		self.states = [s1, s2, s3, s4]
		self.massflow = massflow
		self.Qin = self.massflow * (self.states[0].get_h() - self.states[3].get_h()) / 1000
		self.Qout = self.massflow * (self.states[1].get_h() - self.states[2].get_h()) / 1000
		self.Wt = self.massflow * (self.states[0].get_h() - self.states[1].get_h()) / 1000
		self.Wp = self.massflow * (self.states[3].get_h() - self.states[2].get_h()) / 1000
		self.thermalEfficiency = (self.Wt - self.Wp) / self.Qin

	def __repr__(self):
		toprint = ''
		toprint += 'State   p      T      v      h      s\n'
		for state in self.states:
			ind = self.states.index(state)
			p = state.get_p()
			T = state.get_T()
			v = state.get_v()
			h = state.get_h()
			s = state.get_s()
			toprint += '{:0d}       {:1.2f}  {:2.2f}  {:3.4f}  {:4.2f}  {:5.4f}\n'.format(ind + 1, p, T, v, h, s)
		return toprint

	def get_p(self, s):
		return self.states[s - 1].get_p()
	
	def get_T(self, s):
		return self.states[s - 1].get_T()

	def get_h(self, s):
		return self.states[s - 1].get_h()

	def get_s(self, s):
		return self.states[s - 1].get_s()

	def get_v(self, s):
		return self.states[s - 1].get_v()

	def get_Qin(self):
		return self.Qin

	def get_Qout(self):
		return self.Qout

	def get_Wt(self):
		return self.Wt

	def get_Wp(self):
		return self.Wp

	def get_eff(self):
		return self.thermalEfficiency

	def get_massflow(self):
		return self.massflow

	def show_info(self):
		col = ['State', 'p (MPa)', 'T (C)', 'v (m3/kg)', 'h (kJ/kg)', 's (kJ/kg*K)']
		data = []
		for i in range(len(self.states)):
			row = [f'{int(i + 1)}']
			rowdata = self.states[i]
			row.append(round(rowdata.get_p(), 3))
			row.append(round(rowdata.get_T(), 2))
			row.append(round(rowdata.get_v(), 5))
			row.append(round(rowdata.get_h(), 3))
			row.append(round(rowdata.get_s(), 5))
			data.append(row)
		cycle_df = pd.DataFrame(np.array(data), columns = col).set_index('State')
		
		print(cycle_df.to_string())

	def export_info(self):
		col = ['State', 'p (MPa)', 'T (C)', 'v (m3/kg)', 'h (kJ/kg)', 's (kJ/kg*K)']
		data = []
		for i in range(len(self.states)):
			row = [f'{int(i + 1)}']
			rowdata = self.states[i]
			row.append(round(rowdata.get_p(), 3))
			row.append(round(rowdata.get_T(), 2))
			row.append(round(rowdata.get_v(), 5))
			row.append(round(rowdata.get_h(), 3))
			row.append(round(rowdata.get_s(), 5))
			data.append(row)
		data.append(['Q in:', self.get_Qin(), '-', '-', '-', '-'])
		data.append(['Q out:', self.get_Qout(), '-', '-', '-', '-'])
		data.append(['Wt:', self.get_Wt(), '-', '-', '-', '-'])
		data.append(['Wp:', self.get_Wp(), '-', '-', '-', '-'])
		data.append(['Cycle efficiency:', self.get_eff(), '-', '-', '-', '-'])
		
		cycle_df = pd.DataFrame(np.array(data), columns = col).set_index('State')
		name = input('Enter a name for this csv file: ')
		cycle_df.to_csv(name, sep = '\t')		

	def show_power(self):
		print("Q in:", self.Qin, "MW")
		print("Q out:", self.Qout, "MW")
		print("W turbine:", self.Wt, "MW")
		print("W pump:", self.Wp, "MW")
		print("Thermal efficiency:", self.get_eff())
		print()
		print("First law check:", self.check_first(), "MW")

	def check_first(self):
		return ((self.get_Qin() - self.get_Qout()) - (self.get_Wp() - self.get_Wt()))/1000

class Cycle2(Cycle):
	'''
	Creates a Cycle object with the same information, but with 6 states (two turbines)
	'''
	def __init__(self, s1, s2, s3, s4, s5, s6, massflow):
		self.states = [s1, s2, s3, s4, s5, s6]
		self.massflow = massflow
		self.Qin1 = self.massflow * (self.states[2].get_h() - self.states[1].get_h()) / 1000
		self.Qin2 = self.massflow * (self.states[0].get_h() - self.states[5].get_h()) / 1000
		self.Qin = self.Qin1 + self.Qin2 
		self.Qout = self.massflow * (self.states[3].get_h() - self.states[4].get_h()) / 1000
		self.Wt1 = self.massflow * (self.states[0].get_h() - self.states[1].get_h()) / 1000
		self.Wt2 = self.massflow * (self.states[2].get_h() - self.states[3].get_h()) / 1000
		self.Wt = self.Wt1 + self.Wt2 
		self.Wp = self.massflow * (self.states[5].get_h() - self.states[4].get_h()) / 1000
		
		self.thermalEfficiency = (self.Wt - self.Wp) / self.Qin

def fuelInfo():
	'''
	Stores information for all 4 fuel types
	'''
	# type = [energy/mass (MJ/kg), cost/mass ($/kg), CO2/mass (kg/kg)]
	anthracite = [24, 0.26, 3.3]
	bituminous = [0.0285, 0.07, 2.2]
	subbituminous = [20.5, 0.02, 1.87]
	lignite = [15, 0.02, 1.43]
	fueltypes = {}
	fueltypes['a'] = anthracite
	fueltypes['b'] = bituminous
	fueltypes['s'] = subbituminous
	fueltypes['l'] = lignite
	return fueltypes

def pmConfig():
	'''
	Configures base units for Pyromat library
	'''
	pm.config['unit_temperature'] = 'C'
	pm.config['unit_pressure'] = 'MPa'

def rangeConfig():
	'''
	Configures range of pressure and temperature values for testing
	'''
	hp_range = [1, 10]
	ip_range = [0.1, 1]
	lp_range = [0.01, 0.07]
	
	p_range = [10, 17.7]
	T_range = [600, 600]

	return [hp_range, ip_range, lp_range, p_range, T_range]

def turbine(state, level):
	'''
	Runs a turbine simulation for a given inlet state and turbine pressure level, returning the outlet states as a list (last entry is nearly isentropic)
	'''
	water = state.fluid
	#level 0 for HP, level 1 for IP, level 2 for LP
	if (level == 0):
		multiplier = 1000
	if (level == 1):
		multiplier = 100
	if (level == 2):
		multiplier = 1000
	p_min = int(multiplier * rangeConfig()[level][0]) #in 0.01 MPa for more increments
	p_max = int(multiplier * state.get_p())
	s_test = state.get_s() #initial state entropy
	turbine_states = []
	for i in range(p_min, p_max):
		p_test = i/multiplier
		if abs(water.s(p = p_test, x = 1)) < s_test : #returns all possible saturated states in p range (last entry is nearly isentropic)
			return turbine_states
		T_test = water.T(p = p_test, x = 1)
		appendset = [p_test, T_test]
		turbine_states.append(appendset)
	return turbine_states

def run1(pressure, temp, massflow, allCycles):
	'''
	Runs a one-turbine cycle simulation for a given initial temperature and pressure
	'''
	s1 = State(pressure, temp)
	turbine_states = turbine(s1, 0)
	if len(turbine_states) == 0:
		return #terminates current pressure if no turbine states area available
	s2 = MixState(p = turbine_states[-1][0], T = turbine_states[-1][1], x = 1) 
	s3 = MixState(p = s2.get_p(), T = s1.fluid.T(p = s2.get_p(), x = 0), x = 0)
	s4 = State(T = s1.fluid.T(p = s1.get_p(), s = s3.get_s()), p = s1.get_p())
	allCycles.append(Cycle(s1, s2, s3, s4, massflow))

def run2(pressure, temp, massflow, allCycles):
	'''
	Runs a two-turbine cycle simulation for a given initial temperature and pressure
	'''
	s1 = State(pressure, temp)
	turbine_states1 = turbine(s1, 0) #use to change between hp and ip for first turbine
	if len(turbine_states1) == 0:
		return
	s2 = MixState(p = turbine_states1[-1][0], T = turbine_states1[0][1], x = 1)
	s3 = State(p = s2.get_p(), T = 630) #assuming use of two-casing mitsubishi turbine with reheat temp of 630 C
	turbine_states2 = turbine(s3, 2)
	if len(turbine_states2) == 0:
		return
	s4 = MixState(p = turbine_states2[-1][0], T = turbine_states2[-1][1], x = 1)
	s5 = MixState(p = s4.get_p(), T = s4.fluid.T(p = s4.get_p(), x = 0), x = 0)
	s6 = State(T = s5.get_T(), p = s1.get_p())
	allCycles.append(Cycle2(s1, s2, s3, s4, s5, s6, massflow))

def runPress(temp, massflow, turbines):
	'''
	Runs and collects all possble isentropic cycles at a given temperature and within the configured pressure range
	'''
	mult = 10
	p_min = mult * rangeConfig()[3][0]
	p_max = mult * rangeConfig()[3][1]
	T_min = rangeConfig()[4][0] // mult
	T_max = rangeConfig()[4][1] // mult
	allCycles = []
	for i in range(p_min, int(1 + p_max)):
		pressure = i / mult
		if (turbines == 1):
			run1(pressure, temp, massflow, allCycles)
		if (turbines == 2):
			run2(pressure, temp, massflow, allCycles)
	return allCycles

def analyzeCycle(allCycles_real, operator):
	'''
	Analyzes a list of cycles to find the optimal cycle for a given parameter (power, efficiency, cost, environment)
	'''
	allCycles = allCycles_real.copy()
	variables = []
	if operator == 'power':
		for cycle in allCycles:
			variables.append(cycle.get_Wt())
		sort(allCycles, variables)
		optimized = allCycles[-1]

	if operator == 'efficiency':
		for cycle in allCycles:
			variables.append(cycle.get_eff())
		sort(allCycles, variables)
		optimized = allCycles[-1]
		
	if operator == 'cost':
		for cycle in allCycles:
			variables.append(cycle.get_Qin())
		sort(allCycles, variables)
		optimized = allCycles[0]

	if operator == 'environment':
		for cycle in allCycles:
			variables.append(cycle.get_Qout())
		sort(allCycles, variables)
		optimized = allCycles[0]
	
	return optimized

def sort(allCycles, variables):
	'''
	Used to sort cycles by order of a given parameter
	'''
	if len(variables) <= 1:
		return
	
	mid = len(variables)//2
	var1 = variables[:mid]
	cycle1 = allCycles[:mid]
	var2 = variables[mid:]
	cycle2 = allCycles[mid:]
	
	sort(cycle1, var1)
	sort(cycle2, var2)

	merge(cycle1, cycle2, var1, var2, allCycles, variables)
	return allCycles

def merge(cycle1, cycle2, var1, var2, allCycles, variables):
	'''
	Used in sort function
	'''
	i = j = k = 0
	while i < len(var1) and j < len(var2):
		if var1[i] < var2[j]:
			variables[k] = var1[i]
			allCycles[k] = cycle1[i]
			i += 1
		else:
			variables[k] = var2[j]
			allCycles[k] = cycle2[j]
			j += 1
		k += 1
	while i < len(var1):
		variables[k] = var1[i]
		allCycles[k] = cycle1[i]
		i += 1
		k += 1
	while j < len(var2):
		variables[k] = var2[j]
		allCycles[k] = cycle2[j]
		j += 1
		k += 1

def costAnalysis(cycle, fuel):
	'''
	Analyzes a cycle's costs and environmental impact given a specified fuel type
	'''
	electricity = 0.15 #$/kW*h
	time = 2592000 # seconds in 30 days
	fuelinfo = fuelInfo()[fuel] 
	energymass = fuelinfo[0] # MJ/kg
	costmass = fuelinfo[1] # $/kg
	carbonmass = fuelinfo[2] # kg/kg
	Qin = cycle.get_Qin() # MJ/kg steam
	Wp = cycle.get_Wp() # MJ/kg steam
	massflow = cycle.get_massflow() # kg steam/s
	
	fuelmass = time * Qin / energymass # mass of fuel per 30 days
	fuelcost = fuelmass * costmass
	pumpcost = time * Wp * electricity * 1000 / 3600
	monthlycost = fuelmass + pumpcost
	monthlyemissions = fuelmass * carbonmass

	return monthlycost, monthlyemissions, fuelmass * 6	

pmConfig()
