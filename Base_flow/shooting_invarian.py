import numpy as np
import matplotlib.pyplot as plt
import fuctions as fun
from   scipy.integrate import solve_ivp
from   scipy.optimize import fsolve
from   scipy.interpolate import interp1d
import scipy.integrate as inte
import plotparameters as pp 


# --- Parâmetros do Problema ---
# Adimensionalizados pela corrente de cima
###########################################
# --- Perfil de Velocidade   --- 
# --- upper
Uu = 1.0
# --- lower
Ul = 0.5
###########################################
# --- Perfil de species e chama Z 
# --- upper
Y_F= 1.0
Y_O= 1.0 
# --- n stechiometry mass ratio  
n  = 4.0
#________________________________________
Zu = Y_F 
# --- lower
Zl = -1*Y_O/n # Y_o/n
###########################################
# --- Ração de valor específico para o Combustível (F)
gamma_F= 1.4
Cp_F   = 1.0
# --- Ração de valor específico para o Oxidante (O)
gamma_O= 1.4
Cp_O   = 1.0
# --- Ração de valor específico para o Produtos (P)
gamma_P= 1.4
Cp_P   = 1.0
#Todos são constantes é iguais 
#-------------------------------
# --- Constante do gas para o Combustível
R_F    = 1.0
# --- Constante do gas para o Oxidante
R_O    = 1.0
# --- Constante do gas para os Produtos 
R_P    = 1.0
#-------------------------------
# --- Entalpia de formação Combustível 
Ef_F   = 1.0
# --- Entalpia de formação Oxidante 
Ef_O   = 0.1
# --- Entalpia de formação Produtos 
Ef_P   = 3.0
#-------------------------------
# --- Temperatura corrente de cima 
Tu     = 1.0
# --- Temperatura corrente de baixo 
Tl     = 1.0
#_______________________________
# --- Perfil Entapia Total 
# --- Upper
#Hu =  gamma_F*R_F/(gamma_F-1)*Tu+Y_F*Ef_F+Uu**2/2.0
Hu =  Cp_F*Tu+Y_F*Ef_F+Uu**2/2.0
# --- Lower
#Hl =  gamma_O*R_O/(gamma_O-1)*Tl+Y_O*Ef_O+Ul**2/2.0
Hl =  Cp_O*Tl+Y_O*Ef_O+Ul**2/2.0
print(f'Zu:{Zu},Zl:{Zl}')
print(f'Hu:{Hu},Hl:{Hl}')

#exit()

###########################################
#Reynolds Number 
Re=10000.0
###########################################
#Heat Release 
theta=((n+1)*Ef_P-n*Ef_O-Ef_F)/((1.0/Y_F+n/Y_O)*Cp_F*Tu)

print(r'\theta:%s'%theta)
#exit()

#Malha
# --- -infinity in L
Ll=  -40
# --- infinity in L
Lu=   40

# Chute inicial para o fsolve (apenas para o F''(0))
chute_fpp0 = [0.5]

# Executa o loop de iteração externa do fsolve
# ate encontrar o chuting para fpp0, partindo de un 
# chute initial para fp0 que é solucãop pela invarinça, só tem que 
# ser rescalado.

fpp0_convergido = fsolve(fun.objective_shooting_F, chute_fpp0, args=(Uu, Ul,Lu,Ll))[0]

print(f"-> Convergência alcançada!")
print(f"-> F''(0) original de partida encontrado pelo fsolve: {fpp0_convergido:.6f}")

# ==========================================
# RECONSTRUÇÃO DA SOLUÇÃO FINAL 
# ==========================================
# Executa o passo final para pegar os valores de escala corretos gerados pelo fpp0_convergido
#Reconstruindo a solução, pq o fsolve na devolve isto!
fp0_arbitrary = 1.0
sol_dir_final = solve_ivp(fun.F_system, [0, Lu], [0.0, fp0_arbitrary, fpp0_convergido])
c_final = np.sqrt(Uu / sol_dir_final.y[1][-1])

fp0_scaled_final = fp0_arbitrary * (c_final**2)
fpp0_scaled_final = fpp0_convergido * (c_final**3)

# Integração final para gráficos com malha fina de pontos (t_eval)
l_eval_dir = np.linspace(0, Lu, 200)
l_eval_esq = np.linspace(0, Ll, 200)

sol_dir = solve_ivp(fun.F_system, [0, Lu], [0.0, fp0_scaled_final, fpp0_scaled_final], t_eval=l_eval_dir)
sol_esq = solve_ivp(fun.F_system, [0, Ll], [0.0, fp0_scaled_final, fpp0_scaled_final], t_eval=l_eval_esq)

fpp0_convergido = fsolve(fun.objective_shooting_F, chute_fpp0, args=(Uu, Ul,Lu,Ll))[0]

# Unifica os dados da esquerda para a direita de forma ordenada
#l_total = np.concatenate((sol_esq.t[::-1]   , sol_dir.t))
#u_total = np.concatenate((sol_esq.y[1][::-1], sol_dir.y[1]))
#print(f"Verificação F'(-40): {u_total[0]:.4f} (Alvo: {ul})")
#print(f"Verificação F'(40): {u_total[-1]:.4f} (Alvo: {uu})")

# Unifica os dados da esquerda para a direita de forma ordenada.
# Assím, a solução final para F é 
# O segundo laço começa em 1 para nao repetir o zero
L  = np.concatenate((sol_esq.t[::-1]   , sol_dir.t[1::]))

F  = np.concatenate((sol_esq.y[0][::-1], sol_dir.y[0][1::]))
Fp = np.concatenate((sol_esq.y[1][::-1], sol_dir.y[1][1::]))
Fpp= np.concatenate((sol_esq.y[2][::-1], sol_dir.y[2][1::]))

# 1. Cria uma função que interpola seus dados
F_interpolado   = interp1d(L, F  , kind='cubic', fill_value="extrapolate")
Fp_interpolado  = interp1d(L, Fp , kind='cubic', fill_value="extrapolate")
Fpp_interpolado = interp1d(L, Fpp, kind='cubic', fill_value="extrapolate")


Z,LZ = fun.Z_H_invariant(F_interpolado,Ll,Lu,L,Zl,Zu) 


H,LH = fun.Z_H_invariant(F_interpolado,Ll,Lu,L,Hl,Hu) 

#Then, Mass fraction and temperature can be calculated
#This means that Y_i are function of Z, the to calculate Cp=Cp(Z)
YF = np.maximum(Z, 0)
YO = -n*np.minimum(Z, 0)
YP = 1-YF-YO

#Calculating the constant heat coeficient. 
#cp_F= gamma_F*R_F/(gamma_F-1) 
#cp_O= gamma_O*R_O/(gamma_O-1) 
#cp_P= gamma_P*R_P/(gamma_P-1) 

#Calculando U na mesma malha do Z
U   = Fp_interpolado(LZ)

#print(YF*Cp_F+YO*Cp_O+YP*Cp_P)
#exit()

#Calculando a Temperatura usando a Energia total 
T  = (H-YF*Ef_F-YO*Ef_O-U**2.0/2.0)/(YF*Cp_F+YO*Cp_O+YP*Cp_P)


# Reverse the Howarth's transformation

deta=LZ[1]-LZ[0]

# 1. Run a standard cumulative integration from the very beginning (-40)
#raw_integral = inte.cumulative_trapezoid(T,dx=deta, initial=0)
raw_integral = inte.cumulative_trapezoid(T,x=LZ, initial=0)

# 2. Find the index where eta is exactly 0
zero_index = np.where(LZ == 0)[0][0]

# 3. Subtract the value at eta=0. 
# This forces y*(0) to be exactly 0, makes values before it negative, and values after it positive!
#At the index where $\eta = 0$, raw_integral[zero_index] - raw_integral[zero_index] equals exactly 0. Thus, $y(0) = 0$.
ystar = raw_integral - raw_integral[zero_index]


#Calculing the vorticity thickness
temporal_func = Fpp/T

max_val= max(abs(temporal_func))

x = (Re/2.0)*(max_val/(1-Ul))**2.0


#Calculated x, now y is calculated
y = np.sqrt(2*x/Re)*ystar

#print(x)
#exit()

#plt.style.use('seaborn-poster')

'''

plt.figure(figsize=(10, 6))
plt.plot(Z,y, label=r"Z($y$)", color='red'  , lw=2.5)
plt.xlabel(r"Z($y$)")
plt.ylabel(r'$y$')
plt.title('Base flow')
plt.xlim( 0, 1.5)
plt.ylim(-5, 5)
plt.legend()
plt.grid(True)

plt.figure(figsize=(10, 6))
plt.plot(Z,y, label=r"Z($y$)", color='red'  , lw=2.5)
plt.xlabel(r"Z($y$)")
plt.ylabel(r'$y$')
plt.title('Base flow')
plt.xlim( 0, 1.5)
plt.ylim(-5, 5)
plt.grid(True)


plt.figure(figsize=(10, 6))
plt.plot(H,y, label=r"H($y$)", color='red'  , lw=2.5)
plt.xlabel(r"H($y$)")
plt.ylabel(r'$y$')
plt.title('Base flow')
plt.xlim( 0, 1.5)
plt.ylim(-5, 5)
plt.legend()
plt.grid(True)

'''

size_wg = 0.5
size_hf = 1.0
cmas    = 0.0

plotdef = 'diurnal2' 
    
tama    = pp.plotsize(size_wg,size_hf, cmas,plotdef)

plt.figure()
plt.plot(T,y, label=r"$\mathrm{T(y)}$", color='red'  ,dashes=[1,0], lw=1.5)
plt.plot(U,y, label=r"$\mathrm{u(y)}$", color='blue' ,dashes=[1,0], lw=1.5)
plt.xlabel(r"$\mathrm{T, u}$")
plt.ylabel(r'$\mathrm{y}$')
#plt.title('Base flow')
plt.xlim( 0.4, 1.5)
plt.ylim(-2, 2)
plt.grid(color = 'gray',linewidth=0.5,alpha=0.5,dashes=[1,1,0,0] )
plt.legend(frameon=False)
plt.savefig('%s/%s.pdf'%(pp.out_fig,'base_flow_u_T'),bbox_inches='tight',dpi=200, format='pdf')


#plt.figure(figsize=(10, 6))
#plt.plot(U,y, label=r"U($y$)", color='red'  , lw=2.5)
#plt.xlabel(r"U($y$)")
#plt.ylabel(r'$y$')
#plt.title('Base flow')
#plt.xlim( 0, 1.5)
#plt.ylim(-5, 5)
#plt.legend()
#plt.grid(True)

tama    = pp.plotsize(size_wg,size_hf, cmas,plotdef)
plt.figure()
plt.plot(YF,y, label=r"$\mathrm{Y_F(\mathrm{y})}$", color='red'    ,dashes=[1,0], lw=1.5)
plt.plot(YO,y, label=r"$\mathrm{Y_O(\mathrm{y})}$", color='blue'   ,dashes=[1,0], lw=1.5)
plt.plot(YP,y, label=r"$\mathrm{Y_P(\mathrm{y})}$", color='magenta',dashes=[1,1], lw=1.5)
plt.xlabel(r"$\mathrm{Y_i(y)}$")
plt.ylabel(r'$\mathrm{y}$')
#plt.title('Base flow')
plt.xlim( -0.1, 1.1)
plt.ylim(-2, 2)
#plt.grid(True)
plt.grid(color = 'gray',linewidth=0.5,alpha=0.5,dashes=[1,1,0,0] )
plt.legend(frameon=False)
plt.savefig('%s/%s.pdf'%(pp.out_fig,'base_flow_Y_i'),bbox_inches='tight',dpi=200, format='pdf')
#'''

plt.show()
exit()

plt.figure()
plt.plot(YF,LZ, label=r"Y_F($\gamma$)", color='red'  , lw=2.5)
plt.plot(YO,LZ, label=r"Y_O($\gamma$)", color='blue' , lw=2.5)
plt.plot(YP,LZ, label=r"Y_P($\gamma$)", color='green', lw=2.5)
plt.xlabel(r"Y_i'($\gamma$)")
plt.ylabel(r'$\gamma$')
#plt.title('Base flow')
plt.xlim( -0.1, 1.1)
plt.ylim(-10, 10)
plt.grid(color = 'gray',linewidth=0.5,alpha=0.5,dashes=[1,1,0,0] )
plt.legend()
#plt.grid(True)

plt.figure(figsize=(10, 6))
plt.plot(T,LZ, label=r"T($\gamma$)", color='red'  , lw=2.5)
plt.xlabel(r"T'($\gamma$)")
plt.ylabel(r'$\gamma$')
plt.title('Base flow')
plt.xlim( 0, 10)
plt.ylim(-1, 1)
plt.legend()
plt.grid(True)
plt.show()
exit()


#"""
plt.figure(figsize=(10, 6))
plt.plot(Fp,L, label=r"F'($\gamma$)", color='purple', lw=2.5)
plt.xlabel(r"F'(y)=U")
plt.ylabel(r'$\gamma$')
plt.title('Base flow')
plt.legend()
plt.grid(True)
#plt.show()

plt.figure(figsize=(10, 6))
plt.plot(Z,LZ, label="Z($\gamma$)", color='purple', lw=2.5)
plt.xlabel(r"Z")
plt.ylabel(r'$\gamma$')
plt.title('Base flow')
plt.legend()
plt.grid(True)

plt.figure(figsize=(10, 6))
plt.plot(H,LH, label="H($\gamma$)", color='purple', lw=2.5)
plt.xlabel(r"H")
plt.ylabel(r'$\gamma$')
plt.title('Base flow')
plt.legend()
plt.grid(True)


plt.show()
#"""
