# In[42]:


from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt



# describe the model
def deriv(y, t, N, beta, gamma, delta, theta, my):
    S, E, I, R ,B = y
    dBdt = my * I / N
    dSdt = -beta * (1-B) * S * I / N + theta * R# S(t) – susceptible (de som är mottagliga för infektion).
    #dEdt = beta *(1-B) * S * I / N - gamma * E #Exposed (har blivit utsatta för sjukdomen)
    dIdt = delta * E - gamma * I # I(t) – infected (de som har pågående infektion)
    dRdt = gamma * I - theta * R #Removed (tillfrisknar)
    return dSdt, dIdt, dRdt, dBdt




# describe the parameters
N =  100000 #Totala befolkningen N=s(t)+I(t)+R(t)
D = 7 #Hur länge infektionen varar
gamma = 1.0 / D #Removal (hur många som tillfrisknar per dag)
delta = 1.0 / 5.0 #incubation period of five days
theta = 0.01 / D #tidigarer immuna som insjuknar igen
my = 0.3
R_0 = 2.5 #Reproduktionstalet  
beta = R_0 * gamma #r_0=beta/gamma. antal som smittas per infekterad och per tid (beror på virusets egenskaper samt hur vi beter oss).  
S0,I0, R0, B0 = N-1, 1, 0, 0, 0  # initial conditions: one infected, rest susceptible
 
#Rt = R0 * S(t)/Ntot* (1 – b). b = effekt av policy och beteendeförändringar



t = np.linspace(0, 999, 1000) # Grid of time points (in days)
y0 = S0,I0, R0, B0 # Initial conditions vector
 
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta,theta,my))
S, I, R, B = ret.T




def plotsir(t, S, E, I, R):
  f, ax = plt.subplots(1,1,figsize=(10,4))
  ax.plot(t, S, 'b', alpha=0.7, linewidth=2, label='Susceptible')
  ax.plot(t, E, 'y', alpha=0.7, linewidth=2, label='Exposed')
  ax.plot(t, I, 'r', alpha=0.7, linewidth=2, label='Infected')
  ax.plot(t, R, 'g', alpha=0.7, linewidth=2, label='Recovered')

  ax.set_xlabel('Time (days)')

  ax.yaxis.set_tick_params(length=0)
  ax.xaxis.set_tick_params(length=0)
  ax.grid(b=True, which='major', c='w', lw=2, ls='-')
  legend = ax.legend()
  legend.get_frame().set_alpha(0.5)
  for spine in ('top', 'right', 'bottom', 'left'):
      ax.spines[spine].set_visible(False)
  plt.savefig('Plot.png')
  plt.show();


# plot the graph



plotsir(t, S, E, I, R)


# In[ ]:




