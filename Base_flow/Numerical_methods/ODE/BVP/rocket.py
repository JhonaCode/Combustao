import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

plt.style.use('seaborn-poster')
#matplotlib inline

def objective(v0,y0,yf,t0,tf,t_eval):

    #print(v0[0],y0)
    #print(t0,tf)
    #print(t_eval)

    sol = solve_ivp(F, [t0, tf], \
            [y0, v0[0]], t_eval = t_eval)
    y = sol.y[0]


    diff=-1
    k=0
    for ti in sol.t :

        if tf==5:

            diff=y[k] - yf

        k+=1

    print(v0,diff)

    return diff 


F = lambda t, s: \
  np.dot(np.array([[0,1],[0,-9.8/s[1]]]),s)

y0 = 0
yf = 50

v0 = 1.0

#initial time
t0=0
#final time 
tf=5

t_span = np.linspace(0, 10, 100)
t_eval = np.linspace(t0, tf, 20)

v0, = fsolve(objective,v0,args=(y0,yf,t0,tf,t_eval))
print(v0)


F = lambda t, s: \
  np.dot(np.array([[0,1],[0,-9.8/s[1]]]),s)

                #f    , x0
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
"""
