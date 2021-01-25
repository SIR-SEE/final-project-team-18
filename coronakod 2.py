#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 14:24:39 2021

@author: annabeljerre
"""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

#Hej på dig

# describe the model
def deriv(y, t, N, beta, gamma, delta, theta, my, dead):
    S, E, I, R ,B, D = y
    dBdt = my * I / N  #B står för beteendeförändringar, vi gör så att de beror på andelen infekterande med en faktor my
    dSdt = -beta * (1-B) * S * I / N + theta * R # "Susceptible", dvs antal människor som kan bli smittade
    dEdt = beta *(1-B) * S * I / N - gamma * E #"Exposed", dvs antal människor som har blivit utsatta för sjukdomen
    dIdt = delta * E - gamma * I - dead * I #"Infected", dvs antal människor som har sjukdomen
    dRdt = gamma * I - theta * R #"Recovered", antal människor som har tillfriskat
    dDdt = dead * I #"Dead", antal människor som har dött
    return dSdt, dEdt, dIdt, dRdt, dBdt, dDdt


# Parametrar och konstanter
N =  10000000 #Totala befolkningen i Sverige
delta = 1.0 / 5.0 #Inkubationstiden är 5 dagar
D = delta + 5 #Infektionen varar inkubationsdgar + antal dagar man är sjuk (4 dagar)
gamma =  1 / D #Konstant för "Recovered" per dag, antal som tillfrisknar per dag
dead = 0.01 #Risken att dö av COVID
theta = 0.01 / D #tidigarer immuna som insjuknar igen
my = 10**-12 #konstant för beteendeförändringar
R_0 = 2.5 #Reproduktionstalet  
beta = R_0 * gamma #r_0=beta/gamma. antal som smittas per infekterad och per tid (beror på virusets egenskaper samt hur vi beter oss).  
S0, E0, I0, R0, B0, D0 = N-1, 1, 0, 0, 0, 0 # Initialvärden: en smittad, resten noll



t = np.linspace(0, 1999, 2000) # Antal dagar
y0 = S0, E0, I0, R0, B0, D0 # vektor för initalvärden
 
# Integrerar våra differentialfunktioner
ret = odeint(deriv, y0, t, args=(N, beta, gamma, delta,theta,my,dead))
S, E, I, R, B ,D = ret.T


#plottar graferna
def plotsir(t, S, E, I, R):
  f, ax = plt.subplots(1,1,figsize=(10,4))
  ax.plot(t, S, 'b', alpha=0.7, linewidth=2, label='Susceptible')
  ax.plot(t, E, 'y', alpha=0.7, linewidth=2, label='Exposed')
  ax.plot(t, I, 'r', alpha=0.7, linewidth=2, label='Infected')
  ax.plot(t, R, 'g', alpha=0.7, linewidth=2, label='Recovered')
  ax.plot(t, D, 'k', alpha=0.7, linewidth=2, label='Dead')

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



plotsir(t, S, E, I, R)


# In[ ]:






