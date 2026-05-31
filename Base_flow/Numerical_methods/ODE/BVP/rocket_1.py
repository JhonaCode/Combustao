import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

plt.style.use('seaborn-poster')
#matplotlib inline

F = lambda t, s: \
  np.dot(np.array([[0,1],[0,-9.8/s[1]]]),s)

t_span = np.linspace(0, 10, 100)
y0 = 0
v0 = 25

#initial time
t0=0
#final time 
tf=6
t_eval = np.linspace(t0, tf, 20)


sol1 = solve_ivp(F, [t0, tf], \
                [y0, v0], t_eval = t_eval)
v0 = 35
sol2 = solve_ivp(F, [t0, tf], \
                [y0, v0], t_eval = t_eval)
v0 = 34.5
sol3 = solve_ivp(F, [t0, tf], \
                [y0, v0], t_eval = t_eval)


plt.figure(figsize = (10, 8))
plt.plot(sol1.t, sol1.y[0],color='red', label='first')
plt.plot(sol2.t, sol2.y[0],color='blue',label='second')
plt.plot(sol3.t, sol3.y[0],color='green',label='tirst')
plt.legend()
plt.plot(5, 50, 'ro')
plt.xlabel('time (s)')
plt.ylabel('altitude (m)')
plt.title(f'guessing v={v0} m/s')
plt.show()
