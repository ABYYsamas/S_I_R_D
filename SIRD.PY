import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import product

S0=0.99
I0=0.1 
R0=0
D0=0
beta_for_grid=np.linspace(0.25, 0.5, 8)
Gamma_for_grid=np.linspace(0.08, 0.15, 8)
Mu_for_grid=np.linspace(0.005, 0.015, 8)
Great_params = None
Great_mse = float('inf')
T=50
T = 50  
dataset_size = 90  
dt = T / dataset_size  
data = pd.read_csv("sird_dataset.csv")


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


#J'ai choisi la MSE car elle permet de réduire au maximum l'écart entre les valeurs prédites 
# et les valeurs réelles, avec une pénalité plus élevée pour les erreurs les plus importantes,
# ce qui est particulièrement judicieux dans un cadre épidémiologique. À l'inverse, 
# la CCE est généralement utilisée pour des probabilités de classes discrètes, 
# ce qui ne correspond donc pas aux besoins de notre projet.


def cost(Beta, Gamma, MU, T, dt, S0, I0, R0, D0):
    data = pd.read_csv("sird_dataset.csv") 
    S, I, R, D = methode_euler(Beta, Gamma, MU, T, dt, S0, I0, R0, D0)  
    New_mse_s = np.mean((S - data["Susceptibles"]) ** 2)
    New_mse_i = np.mean((I - data["Infectés"]) ** 2)
    New_mse_r = np.mean((R - data["Rétablis"]) ** 2)
    New_mse_d = np.mean((D - data["Décès"]) ** 2)
    M_S_E = (New_mse_s, New_mse_i, New_mse_r, New_mse_d)

    return New_mse_s, New_mse_i, New_mse_r, New_mse_d

for Beta, Gamma, MU in product(beta_for_grid, Gamma_for_grid, Mu_for_grid):
    current_mse = cost(Beta, Gamma, MU, T, dt, S0, I0, R0, D0)
    total_mse = np.mean(current_mse)
    if total_mse < Great_mse:
        Great_mse = total_mse
        Great_params = (Beta, Gamma, MU)


print("Meilleurs paramètres:", Great_params)
print("MSE minimum:", Great_mse)

S, I, R, D = methode_euler(*Great_params, T, dt, S0, I0, R0, D0)
M_S_E = cost(*Great_params, T, dt, S0, I0, R0, D0)


print("MSE pour les Susceptibles:", M_S_E[0])
print("MSE pour les Infectés:", M_S_E[1])
print("MSE pour les Rétablis:", M_S_E[2])
print("MSE pour les Décès:", M_S_E[3])


temps = np.linspace(0, T, len(S))


plt.figure(figsize=(12, 6))
plt.plot(temps, S, label="susceptible")
plt.plot(temps, I, label="infect")
plt.plot(temps, R, label="Recovered")
plt.plot(temps, D, label="Deceased")
plt.xlabel('Time (days)')
plt.ylabel('Percentage of Population')
plt.legend()
plt.show()
