import cantera as ct
import numpy as np
import math as math
import matplotlib.pyplot as plt

"""Approach 3: Varying Ta3, t and m_N2O4"""

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
m_N2O4_vec = np.linspace(0.006,0.02,100)
Ta3_vec = [5+273, 10+273, 15+273, 25+273]


"""CREATE FLUID: N2O4"""
wf = ct.Solution('./kth/src/TDA/N2O4_eq.yaml') # wf = working fluid

"""Create subplots"""
fig = ['fig1', 'fig2', 'fig3', 'fig4']
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
axs = [ax1, ax2, ax3, ax4]
fig.suptitle("Net power")

fig1 = ['fig1', 'fig2', 'fig3', 'fig4']
fig1, ((ax11,ax22),(ax33,ax44)) = plt.subplots(2,2)
axs1 = [ax11, ax22, ax33, ax44]
fig1.suptitle("Power production")

fig2 = ['fig1', 'fig2', 'fig3', 'fig4']
fig2, ((ax111,ax222),(ax333,ax444)) = plt.subplots(2,2)
axs2 = [ax111, ax222, ax333, ax444]
fig2.suptitle("Cooler")


k = 0
for Ta3 in Ta3_vec:
    """SAVING VALUES"""
    Net_power_vec = [[np.NaN for _ in range(len(m_N2O4_vec))] for _ in range(len(t_vec))]
    Pprod = [[np.NaN for _ in range(len(m_N2O4_vec))] for _ in range(len(t_vec))]
    Qcooler = [[np.NaN for _ in range(len(m_N2O4_vec))] for _ in range(len(t_vec))]

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
            wf.equilibrate('HP')
            T2 = wf.T

            if T2 > Ta2:
                j += 1
                continue

            """State 3: compressing fluid"""
            Compressor = compress(wf, P3, eta_comp) # Defined as positive value
            wf.equilibrate('HP')
            H3 = wf.h
            T3 = wf.T

            if T3 > Ta2:
                j += 1
                continue

            """State 4: heating of fluid in tank 1"""
            H4 = H3 + Q1_N2O4/m_N2O4
            wf.HP = H4, P4
            wf.equilibrate('HP')
            T4 = wf.T

            if T4 > Ta2:
                j += 1
                continue

            """State 5: expanding fluid"""
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

            Net_power_vec[i][j] = Net_power
            Pprod[i][j] = (Turbine-Compressor)*m_N2O4
            Qcooler[i][j] = Cooler*m_N2O4

            j += 1
        i += 1
    axs[k].set_title("Ta3 = " + str(Ta3) + " [K]")
    axs[k].set(xlabel = 'Mass flowrate of working fluid [kg/s]', ylabel = 'Residence time [min]')
    im = axs[k].imshow(Net_power_vec, extent=[m_N2O4_vec[0], m_N2O4_vec[-1] ,t_vec[-1], t_vec[0]], aspect = "auto", cmap="rainbow")
    cbar = plt.colorbar(im, label="[W]")

    axs1[k].set_title("Ta3 = " + str(Ta3) + " [K]")
    axs1[k].set(xlabel = 'Mass flowrate of working fluid [kg/s]', ylabel = 'Residence time [min]')
    im1 = axs1[k].imshow(Pprod, extent=[m_N2O4_vec[0], m_N2O4_vec[-1] ,t_vec[-1], t_vec[0]], aspect = "auto", cmap="rainbow")
    cbar1 = plt.colorbar(im1, label="[W]")

    axs2[k].set_title("Ta3 = " + str(Ta3) + " [K]")
    axs2[k].set(xlabel = 'Mass flowrate of working fluid [kg/s]', ylabel = 'Residence time [min]')
    im2 = axs2[k].imshow(Qcooler, extent=[m_N2O4_vec[0], m_N2O4_vec[-1] ,t_vec[-1], t_vec[0]], aspect = "auto", cmap="rainbow")
    cbar2 = plt.colorbar(im2, label="[W]")

    k += 1

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.2, 
                    hspace=0.4)

plt.show() 