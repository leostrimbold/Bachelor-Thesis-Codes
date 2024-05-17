import cantera as ct
import numpy as np
import math as math
import matplotlib.pyplot as plt

"""Approach 1: Ta3 as variables. m_acid = 0.0005 kg/s, m_N2O4 = 0.0085 kg/s"""

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

"""CONSTANTS OUTER LOOP"""
m_acid = 0.0005 # kg/s mass flowrate of acid solution 
cp_h2o = 4.18*1000 # J/kg*K specific heat capacity of water 
density_h2o = 997.05 # kg/m^3 density of water

PUS = 1000 # [W] Power of ultrasound probe 

Ta1 = 25 + 273 # Temp of acid before tank 1 
Ta2 = 50 + 273 # Temp of acid after tank 1 

"""REACTION MIXTURE (rm) VOLUME TANK 1"""
h_rm = 0.21 # m height of reaction mixture
d_rm = 0.3 # m diameter of tank 1
r_rm = d_rm/2 # m radius of tank 1
V_rm = math.pi*r_rm**2*h_rm # m^3 volume of reaction mixture

"""CONSTANT INNER LOOP"""
m_N2O4 = 0.0085 # kg/s mass flowrate of working fluid 
eta_turb = 0.7
eta_comp = 0.7

P1 = 1e5 
T1 = -2 + 273 

P2 = P1
P5 = P1

P3 = 3.3e5 
P4 = P3

"""VARIABLE"""
Ta3_vec = range(-2+273,26+273,1) #Temp of acid after tank 2

"""SAVING VALUES"""
Net_power_vec = [np.NaN for _ in range(len(Ta3_vec))]
Q_ratio_vec = [np.NaN for _ in range(len(Ta3_vec))]
W_turbine = [np.NaN for _ in range(len(Ta3_vec))]
W_compressor = [np.NaN for _ in range(len(Ta3_vec))]
Q_cooler = [np.NaN for _ in range(len(Ta3_vec))]
Qa3_vec = [np.NaN for _ in range(len(Ta3_vec))]
totheat = [np.NaN for _ in range(len(Ta3_vec))]
Pprod = [np.NaN for _ in range(len(Ta3_vec))]

"""Creating N2O4 fluid"""
wf = ct.Solution('./kth/src/TDA/N2O4_eq.yaml')

i = 0
for Ta3 in Ta3_vec:
    """ENERGY BALANCES OUTER LOOP"""
    Q1_acid = m_acid*cp_h2o*(Ta2-Ta1) # Heat added to acid in tank 1 
    Q2_acid = m_acid*cp_h2o*(Ta3-Ta2) # Heat removed from acid in tank 2 
    Q3_acid = m_acid*cp_h2o*(Ta1 - Ta3) # Heat added to acid in tank 3 

    """ENERGY BALANCES INNER LOOP"""
    Q1_N2O4 = PUS - Q1_acid # Heat added to N2O4 in tank 1  
    Q2_N2O4 = -Q2_acid # Heat added to N2O4 in tank 2 
    Q_ratio = Q1_N2O4/Q2_N2O4 
    totheat[i] = Q1_N2O4 + Q2_N2O4

    """State 1 - initial state"""
    wf.TP = T1, P1
    wf.equilibrate('TP')
    H1 = wf.h

    """State 2 - heating fluid in tank 2"""
    H2 = H1 + Q2_N2O4/m_N2O4
    wf.HP = H2, P2
    wf.equilibrate('HP')
    T2 = wf.T

    if T2 > Ta2:
        i += 1
        continue

    """State 3 - compressing fluid"""
    Compressor = compress(wf, P3, eta_comp) # Defined as positive value
    wf.equilibrate('HP')
    H3 = wf.h
    T3 = wf.T

    if T3 > Ta2:
        i += 1
        continue

    """State 4 - heating fluid in tank 1"""
    H4 = H3 + Q1_N2O4/m_N2O4
    wf.HP = H4, P4
    wf.equilibrate('HP')
    T4 = wf.T

    if T4 > Ta2:
        i += 1
        continue

    """State 5 - expanding fluid"""
    Turbine = expand(wf, P5, eta_turb) # Defined as positive value
    wf.equilibrate('HP')
    H5 = wf.h
    T5 = wf.T

    Cooler = H5-H1 # Defined as positive value

    """Energy calculations"""
    Net_power_spent_OL = PUS + Q3_acid # [W], OL = OUTER LOOP
    Net_power_spent_IL = (Compressor + Cooler)*m_N2O4 # J/kg * kg/s = J/s = [W], IL = INNER LOOP
    Net_power_generated_IL = Turbine*m_N2O4 # J/kg * kg/s = J/s = [W], IL = INNER LOOP

    Net_power = Net_power_generated_IL - Net_power_spent_IL - Net_power_spent_OL
    
    Net_power_vec[i] = Net_power
    Q_ratio_vec[i] = Q_ratio
    W_turbine[i] = Turbine*m_N2O4
    W_compressor[i] = Compressor*m_N2O4
    Q_cooler[i] = Cooler*m_N2O4
    Qa3_vec[i] = Q3_acid
    Pprod[i] = (Turbine - Compressor)*m_N2O4

    if i == 0:
        print("Ta3 = ", Ta3, " yields: ", "Net power = ", Net_power)
        print("                     Q_cooler =", Q_cooler[i])
        print("                     Turbine = ", W_turbine[i])
        print("                     Compressor = ", W_compressor[i])
        print("                     Qa3 = ", Qa3_vec[i])
        print("                     Q_ratio = ", Q_ratio)
        
    
    if i == len(Ta3_vec)-1:
        print("Ta3 = ", Ta3, " yields: ", "Net power = ", Net_power)
        print("                     Q_cooler =", Q_cooler[i])
        print("                     Turbine = ", W_turbine[i])
        print("                     Compressor = ", W_compressor[i])
        print("                     Qa3 = ", Qa3_vec[i])
        print("                     Q_ratio = ", Q_ratio)
    i += 1

fig1 = plt.figure()
plt.title("Net power")
plt.plot(Ta3_vec,Net_power_vec)
plt.xlabel("Ta3 [K]")
plt.ylabel("[W]")

fig2 = plt.figure()
plt.title("Q_ratio = Q1/Q2")
plt.plot(Ta3_vec,Q_ratio_vec)
plt.xlabel("Ta3 [K]")
plt.ylabel("[-]")

fig3, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
ax1.plot(Ta3_vec,W_turbine)
ax1.set_title('a) Turbine')
ax1.set(xlabel = 'Ta3 [K]', ylabel = '[W]')

ax2.plot(Ta3_vec,W_compressor)
ax2.set_title('b) Compressor')
ax2.set(xlabel = 'Ta3 [K]', ylabel = '[W]')

ax3.plot(Ta3_vec, Q_cooler)
ax3.set_title('c) Cooler')
ax3.set(xlabel = 'Ta3 [K]', ylabel = '[W]')

ax4.plot(Ta3_vec, Qa3_vec)
ax4.set_title('d) Qa3')
ax4.set(xlabel = 'Ta3 [K]', ylabel = '[W]')

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.2, 
                    hspace=0.4)

fig4 = plt.figure()
plt.title("Qtot = Q1 + Q2")
plt.plot(Ta3_vec,totheat)
plt.xlabel("Ta3 [K]")
plt.ylabel("[W]")

fig5 = plt.figure()
plt.title("Power production")
plt.plot(Ta3_vec,Pprod)
plt.xlabel("Ta3 [K]")
plt.ylabel("[W]")

plt.show()