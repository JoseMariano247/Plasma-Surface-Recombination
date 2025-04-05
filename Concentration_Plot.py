# Import necessary libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the data from output.txt
data_8 = pd.read_csv("Real_Test_MC.txt", delimiter='\t')
data_9 = pd.read_csv("Real_Test_RK.txt", delimiter='\t')
data_10 = pd.read_csv("recomb_prob.txt", delimiter='\t')

time          = data_10.iloc[:, 0]
gamma_ER      = data_10.iloc[:, 1]
gamma_LHS     = data_10.iloc[:, 2]
gamma_LHF     = data_10.iloc[:, 3]       
gamma_total   = data_10.iloc[:, 4]

time_8 = data_8.iloc[:, 0]
A_r_2    = data_8.iloc[:, 1]
Fv_r_2   = data_8.iloc[:, 2]
Af_r_2   = data_8.iloc[:, 3]       
Sv_r_2   = data_8.iloc[:, 4]   
As_r_2   = data_8.iloc[:, 5]      
A2_r_2   = data_8.iloc[:, 6]

R1 = data_8.iloc[:, 7]
R2 = data_8.iloc[:, 8]
R3 = data_8.iloc[:, 9]
R4 = data_8.iloc[:, 10]
R5 = data_8.iloc[:, 11]
R6 = data_8.iloc[:, 12]
CR7 = data_8.iloc[:, 13]

time_9 = data_9.iloc[:, 0]
A_r_3    = data_9.iloc[:, 1]
Fv_r_3   = data_9.iloc[:, 2]
Af_r_3   = data_9.iloc[:, 3]       
Sv_r_3   = data_9.iloc[:, 4]   
As_r_3   = data_9.iloc[:, 5]      
A2_r_3   = data_9.iloc[:, 6]

R12 = data_8.iloc[:, 7]
R22 = data_8.iloc[:, 8]
R32 = data_8.iloc[:, 9]
R42 = data_8.iloc[:, 10]
R52 = data_8.iloc[:, 11]
R62 = data_8.iloc[:, 12]
CR72 = data_8.iloc[:, 13]

plt.figure(figsize=(9, 16))

plt.plot(time, gamma_ER/np.max(gamma_total), linestyle='-', color='black', label='$\gamma_{ER}$')
plt.plot(time, gamma_LHS/np.max(gamma_total), linestyle='-', color='cyan', label='$\gamma_{LH}^S$')
plt.plot(time, gamma_LHF/np.max(gamma_total), linestyle='-', color="palevioletred", label='$\gamma_{LH}^F$')
plt.plot(time, gamma_total/np.max(gamma_total), linestyle='-', color="orange", label='$\gamma_{total}$')

plt.xlabel('Tw [K]', fontsize=18)
plt.ylabel('$\gamma$', fontsize=18)
plt.yscale('log')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(200, 400)
#plt.title('Concentration vs Time')
plt.legend(fontsize=14)
plt.grid(True)
plt.show()

plt.figure(figsize=(9, 16))

plt.plot(time_8, A_r_2/np.max(A_r_2), linestyle='-', color='black', label='[$A$]')
#plt.plot(time_9, A_r_3/np.max(A_r_3), linestyle='--', color='cyan', label='[$A$]')
plt.plot(time_8, Fv_r_2/np.max(Fv_r_2), linestyle='-', color="palevioletred", label='[$F_v$]')
#plt.plot(time_9, Fv_r_3/np.max(Fv_r_3), linestyle='--', color="black", label='[$F_v$]')
plt.plot(time_8, Af_r_2/np.max(Af_r_2), linestyle='-', color='blue', label='[$A_f$]')
#plt.plot(time_9, Af_r_3/np.max(Af_r_3), linestyle='--', color='palevioletred', label='[$A_f$]')
plt.plot(time_8, Sv_r_2/np.max(Sv_r_2), linestyle='-', color='lime', label='[$S_v$]')
#plt.plot(time_9, Sv_r_3/np.max(Sv_r_3), linestyle='--', color='blue', label='[$S_v$]')
plt.plot(time_8, As_r_2/np.max(As_r_2), linestyle='-', color='orange', label='[$A_s$]')
#plt.plot(time_9, As_r_3/np.max(As_r_3), linestyle='--', color='lime', label='[$A_s$]')
plt.plot(time_8, A2_r_2/np.max(A2_r_2), linestyle='-', color='cyan', label='[$A_2$]')
#plt.plot(time_9, A2_r_3/np.max(A2_r_3), linestyle='--', color='orange', label='[$A_2$]')

plt.xlabel('Time [s]', fontsize=18)
plt.ylabel('Concentration', fontsize=18)
plt.xscale('log')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(10e-16, 10e-11)
#plt.title('Concentration vs Time')
plt.legend(fontsize=14)
plt.grid(True)
plt.show()



plt.figure(figsize=(9, 16))

plt.plot(time_8, R1/np.max(R1), linestyle='-', color='black', label='R1')
#plt.plot(time_8, R12/np.max(R12), linestyle='--', color='cyan', label='R1')
plt.plot(time_8, R2/np.max(R2), linestyle='-', color="palevioletred", label='R2')
#plt.plot(time_8, R22/np.max(R22), linestyle='--', color="black", label='R2')
plt.plot(time_8, R3/np.max(R3), linestyle='-', color='blue', label='R3')
#plt.plot(time_8, R32/np.max(R32), linestyle='--', color='palevioletred', label='R3')
plt.plot(time_8, R4/np.max(R4), linestyle='-', color='yellow', label='R4')
#plt.plot(time_8, R42/np.max(R42), linestyle='--', color='blue', label='R4')
plt.plot(time_8, R5/np.max(R5), linestyle='-', color='orange', label='R5')
#plt.plot(time_8, R52/np.max(R52), linestyle='--', color='yellow', label='R5')
plt.plot(time_8, R6/np.max(R6), linestyle='-', color='violet', label='R6')
#plt.plot(time_8, R62/np.max(R62), linestyle='--', color='orange', label='R6')
plt.plot(time_8, CR7/np.max(CR7), linestyle='-', color='cyan', label='R7')
#plt.plot(time_8, CR72/np.max(CR72), linestyle='--', color='violet', label='R7')

plt.xlabel('Time [s]', fontsize=18)
plt.ylabel('Reaction Rate', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xscale('log')
plt.xlim(10e-16, 10e-11)
#plt.title('Concentration vs Time')
plt.legend(fontsize=14)
plt.grid(True)
plt.show()
