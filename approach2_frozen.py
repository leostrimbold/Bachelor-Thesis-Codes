import cantera as ct
import numpy as np
import math as math
import matplotlib.pyplot as plt

"""Approach 2: frozen conditions"""

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

PUS = 1000 # [W] Power of ultrasound probe 

Ta1 = 25 + 273 # Temp of acid before tank 1 
Ta2 = 50 + 273 # Temp of acid after tank 1 
Ta3 = 5 + 273 # Temp of acid after tank 2


"""CONSTANT INNER LOOP"""
eta_turb = 0.7
eta_comp = 0.7 

P1 = 1e5 
T1 = -2 + 273 

P2 = P1 
P5 = P1 

P3 = 3.3e5 
P4 = P3 

"""VARIABLES"""
t_vec = np.linspace(1440,1,100) # min, residency time
m_N2O4_vec = np.linspace(0.02,0.1,100)

"""SAVING VALUES"""
Net_power_vec = [[np.NaN for _ in range(len(m_N2O4_vec))] for _ in range(len(t_vec))]

"""CREATE FLUID: N2O4"""
wf = ct.Solution('./kth/src/TDA/N2O4_eq.yaml') # wf = working fluid

i = 0
for t in t_vec:
    j = 0
    q = V_rm/(t*60) # m^3/s
    m_acid = q*density_h2o # kg/s
    
    """ENERGY BALANCES OUTER LOOP: ACID SOLUTION"""
    Q1_acid = m_acid*cp_h2o*(Ta2-Ta1) # Heat added to acid in tank 1 
    Q2_acid = m_acid*cp_h2o*(Ta3-Ta2) # Heat removed from acid in tank 2 
    Q3_acid = m_acid*cp_h2o*(Ta1 - Ta3) # Heat added to acid in tank 3 

    """ENERGY BALANCES INNER LOOP: N2O4"""
    Q1_N2O4 = PUS - Q1_acid # Heat added to N2O4 in tank 1 
    Q2_N2O4 = -Q2_acid # Heat added to N2O4 in tank 2 

    for m_N2O4 in m_N2O4_vec:
        """State 1 - initial state"""
        wf.TP = T1, P1
        wf.equilibrate('TP')
        H1 = wf.h

        """State 2: heating of fluid in tank 2"""
        H2 = H1 + Q2_N2O4/m_N2O4
        wf.HP = H2, P2
        T2 = wf.T

        if T2 > Ta2:
            j += 1
            continue

        """State 3: compressing fluid"""
        Compressor = compress(wf, P3, eta_comp) # Defined as positive value
        H3 = wf.h
        T3 = wf.T

        if T3 > Ta2:
            j += 1
            continue

        """State 4: heating of fluid in tank 1"""
        H4 = H3 + Q1_N2O4/m_N2O4
        wf.HP = H4, P4
        T4 = wf.T

        if T4 > Ta2:
            j += 1
            continue

        """State 5: expanding fluid"""
        Turbine = expand(wf, P5, eta_turb) # Defined as positive value
        H5 = wf.h
        T5 = wf.T

        Cooler = H5-H1 # Defined as positive value

        """Energy calculations"""
        Net_power_spent_OL = PUS + Q3_acid # [W], OL = OUTER LOOP
        Net_power_spent_IL = (Compressor + Cooler)*m_N2O4 # J/kg * kg/s = J/s = [W], IL = INNER LOOP
        Net_power_generated_IL = Turbine*m_N2O4 # J/kg * kg/s = J/s = [W], IL = INNER LOOP

        Net_power = Net_power_generated_IL - Net_power_spent_IL - Net_power_spent_OL

        Net_power_vec[i][j] = Net_power
        
        j += 1
    i += 1

fig1 = plt.figure()
plt.title("Net power (frozen working fluid)")
plt.imshow(Net_power_vec, extent=[m_N2O4_vec[0], m_N2O4_vec[-1] ,t_vec[-1],t_vec[0]], aspect = "auto", cmap = "rainbow")
plt.colorbar(label="[W]")
plt.xlabel("Mass flowrate of working fluid [kg/s]")
plt.ylabel("Residence time [min]")

plt.show()