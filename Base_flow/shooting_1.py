import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

plt.style.use('seaborn-poster')
#matplotlib inline

def objective(guess,guess2,u,l,l_eval):

    #S0[0]=\psi, at l=0 does not cross the center line 
    #guessF'(0), and F''(0)
    
    up = fsolve(objective,guess2,args=(fpp0,u[0],[l0,lf],l_eval))

    S0=[0,guess[0],guess2]

    sol = solve_ivp(F, [l[0],l[1]], \
            S0, t_eval = None)
            #S0, t_eval = l_eval)

    #the solution at $lf=infty
    y = sol.y[1]

    diff=y[-1] - u

    print(u,diff)

    #"""
    fig = plt.figure(figsize=(10, 6))
    plt.legend()
    #plt.plot(5, 50, 'ro')
    plt.xlabel('U')
    plt.ylabel('L')
    #plt.title(f'guessing v={v0} m/s')
    plt.plot(sol.y[1],sol.t, label=diff)
    plt.show()
    #"""

    return diff


F = lambda t, s: \
  np.dot(np.array([[0,1,0],[0,0,1],[-s[2],0,0]]),s)


#upper bc
up = 1.0 
#lower bc
ul = 0.5 #u2/u1 

#guess values
#firt derivative in L=0
fp0=1 
#second derivative in L=0
fpp0=0.1


#initial time
l0=0
#final time 
lf=40

#Solving upper part until infinity

l_span = np.linspace(l0, lf, 100)
l_eval = np.linspace(l0, lf, 20)

#guess=fp0
#up = fsolve(objective,guess,args=(fpp0,up,[l0,lf],l_eval))

guess=fpp0
ul = fsolve(objective,guess,args=(fp0,[up,ul],[l0,lf],l_eval))

exit()

#"""
