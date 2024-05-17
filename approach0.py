import cantera as ct
import math as math

"""Approach 0: Ta3 = 278 K, m_acid = 0.0005 kg/s, m_N2O4 = 0.0085 kg/s"""

def compress(fluid, p_out, eta):
    """Adiabatically compress an ideal gas to pressure p_out, using
    a compressor with isentropic efficiency eta."""
    h_in = fluid.h
    s_in = fluid.s

    fluid.SP =s_in, p_out

    h_out_s = fluid.h # isentropic enthalpy OUT
    isentropic_work = h_out_s - h_in
    actual_work = isentropic_work * eta
    h_out = h_in + actual_work
    fluid.HP = h_out, p_out

    return actual_work


def expand(fluid, p_out, eta):
    """Adiabatically expand an ideal gas to pressure p_out, using
    a turbine with isentropic efficiency eta."""
    h_in = fluid.h
    s_in = fluid.s

    fluid.SP =s_in, p_out

    h_out_s = fluid.h # isentropic enthalpy OUT
    isentropic_work = h_in - h_out_s
    actual_work = isentropic_work * eta
    h_out = h_in - actual_work
    fluid.HP = h_out, p_out
    
    return actual_work

def printState(n, fluid):
    print('\n***************** State {0} ******************'.format(n))
    print(fluid.report())

"""REACTION MIXTURE (rm) VOLUME TANK 1"""
h_rm = 0.21 # m height of reaction mixture
d_rm = 0.3 # m diameter of tank 1
r_rm = d_rm/2 # m radius of tank 1
V_rm = math.pi*r_rm**2*h_rm # m^3 volume of reaction mixture

"""CONSTANTS OUTER LOOP"""
cp_h2o = 4.18*1000 # J/kg*K specific heat capacity of water
D = 0.0068 # m inner diameter of all pipes 
r = D/2 # radius of pipes
v = 8.917e-7 # m^2/s kinematic viscocity of water
A_cross = math.pi*r**2 # Cross area of pipes
density_h2o = 997.05 # kg/m^3 density of water

m_acid = 0.0005 # kg/s mass flowrate of acid solution

PUS = 1000 # [W] Power of ultrasound probe 

Ta1 = 25 + 273 # Temp of acid before tank 1
Ta2 = 50 + 273 # Temp of acid after tank 1 
Ta3 = 5 + 273 # Temp of acid after tank 2 

"""Reynolds number acid solution"""
q = m_acid/density_h2o # m^3/s VOLUMETRIC FLOWRATE
u = q/A_cross # m/s VELOCITY
Re_acid = u*D/v

"""Residence time"""
t = V_rm/q/60

"""Energy balances outer loop"""
Q1_acid = m_acid*cp_h2o*(Ta2-Ta1) # Heat added to acid in tank 1 
Q2_acid = m_acid*cp_h2o*(Ta3-Ta2) # Heat removed from acid in tank 2 
Q3_acid = m_acid*cp_h2o*(Ta1 - Ta3) # Heat added to acid in tank 3 

"""Energy balances inner loop"""
Q1_N2O4 = PUS - Q1_acid # Heat added to N2O4 in tank 1 
Q2_N2O4 = -Q2_acid # Heat added to N2O4 in tank 2 
Q_ratio = Q1_N2O4/Q2_N2O4 

"""Inner loop values"""
m_N2O4 = 0.0085 # kg/s mass flowrate of working fluid 

eta_turb = 0.7
eta_comp = 0.7

P1 = 1e5 
T1 = -2 + 273 

P2 = P1
P5 = P1

P3 = 3.3e5 
P4 = P3

"""Simulate N2O4 fluid"""
wf = ct.Solution('./kth/src/TDA/N2O4_eq.yaml')

"""State 1 - initial state"""
wf.TP = T1, P1
wf.equilibrate('TP')
H1 = wf.h
printState(1,wf)

"""State 2: heating fluid in tank 2"""
H2 = H1 + Q2_N2O4/m_N2O4
wf.HP = H2, P2
wf.equilibrate('HP')
T2 = wf.T
printState(2, wf)

"""State 3: compressing fluid"""
Compressor = compress(wf, P3, eta_comp) # Defined as positive value
wf.equilibrate('HP')
H3 = wf.h
T3 = wf.T
printState(3, wf)

"""State 4: heating fluid in tank 2"""
H4 = H3 + Q1_N2O4/m_N2O4
wf.HP = H4, P4
wf.equilibrate('HP')
T4 = wf.T
printState(4, wf)

"""State 5: expanding fluid"""
Turbine = expand(wf, P5, eta_turb) # Defined as positive value
wf.equilibrate('HP')
H5 = wf.h
T5 = wf.T
printState(5, wf)
print('\n',"Inner Loop complete")

Cooler = H5-H1 # Defined as positive value

"""Energy calculations"""
Net_power_spent_OL = PUS + Q3_acid # [W], OL = OUTER LOOP
Net_power_spent_IL = (Compressor + Cooler)*m_N2O4 # J/kg * kg/s = J/s = [W], IL = INNER LOOP
Net_power_generated_IL = Turbine*m_N2O4 # J/kg * kg/s = J/s = [W], IL = INNER LOOP

Net_power = Net_power_generated_IL - Net_power_spent_IL - Net_power_spent_OL
print("Net power = ", Net_power)
print("Turbine = ", Turbine*m_N2O4)
print("Compressor = ", Compressor*m_N2O4)
print("Q1 = ", Q1_N2O4)
print("Q2 = ", Q2_N2O4)
print("Cooler = ", Cooler*m_N2O4)
print("Q3 = ", Q3_acid)
print("Reynolds = ", Re_acid)
print("Residence time = ", t, " [min]")
print("Q_ratio = ", Q_ratio)
print("Total heat added = ", Q1_N2O4+Q2_N2O4, " [W]")
print("Power production = ", (Turbine-Compressor)*m_N2O4, " [W]")