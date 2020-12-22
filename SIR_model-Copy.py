# In[42]:


from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# In[43]:


# describe the model
def deriv(y, t, N, beta, gamma, delta):
    S, E, I, R = y
    dSdt = -beta * S * I / N # S(t) – susceptible (de som är mottagliga för infektion).
    dEdt = beta * S * I / N - gamma * E
    dIdt = delta * E - gamma * I # I(t) – infected (de som har pågående infektion)
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt


# In[44]:


# describe the parameters
N =  2283 #Totala befolkningen N=s(t)+I(t)+R(t)
D = 4.0 #infections last four days
gamma = 1.0 / D #Reoval rate (Hur många som tillfrisknar)
delta = 1.0 / 5.0 #incubation period of five days
R_0 = 2.5 #Reproduktionstalet  
beta = R_0 * gamma #r_0=beta/gamma. antal som smittas per infekterad och per tid (beror på virusets egenskaper samt hur vi beter oss).  
S0, E0, I0, R0 = N-1, 1, 0, 0  # initial conditions: one infected, rest susceptible
 
#Rt = R0 * S(t)/Ntot* (1 – b). b = effekt av policy och beteendeförändringar

# In[45]:


t = np.linspace(0, 99, 100) # Grid of time points (in days)
y0 = S0, E0, I0, R0 # Initial conditions vector
 
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta))
S, E, I, R = ret.T


# In[46]:


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

# In[47]:



plotsir(t, S, E, I, R)


# In[ ]:




