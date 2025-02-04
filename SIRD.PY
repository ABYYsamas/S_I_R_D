import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt 

def methode_euler(Beta, Gamma, MU, T, dt, S0, I0, R0, D0,):
    itération = int(T / dt)
    S=[S0]
    I=[I0]
    R=[R0]
    D=[D0]
    
    for T in range (1, itération ):
        New_s=S[-1]-Beta*S[-1]*I[-1] *dt
        New_i=I[-1]+Beta*S[-1]*I[-1]-MU*I[-1]-Gamma*I[-1] *dt
        New_r=R[-1]+Gamma*I[-1] *dt
        New_d=D[-1]+MU*I[-1] *dt
        
        S.append(New_s)  
        I.append(New_i)  
        R.append(New_r)  
        D.append(New_d)  
     
    return S, I, R, D

Beta=0.5
Gamma=0.15
MU=0.15
S0=0.99
I0=0.1 
R0=0
D0=0
T = 50 
dt = 0.1

S, I, R, D = methode_euler(Beta, Gamma, MU, T, dt, S0, I0, R0, D0)
temps = np.linspace(0, T, len(S))


plt.plot(temps, S, label="Susceptibles", color='blue', linestyle='-')
plt.plot(temps, I, label="Infectés", color='red', linestyle='--')
plt.plot(temps, R, label="Rétablis", color='green', linestyle='-.')
plt.plot(temps, D, label="Décès", color='black', linestyle=':')
plt.xlabel('Time (days)')
plt.ylabel('Percentage of Population')
plt.show()
