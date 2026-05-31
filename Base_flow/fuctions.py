import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# Como Z e H tem a mesma forma, o método de solução é o mesmo 
def Z_H_invariant(F_inter,ll,lu,L,zl,zu):


    #Resolvido F, podemos resolve Z.
    # ==========================================================
    # PASSO 1: Integração Teste Arbitrária (De -40 até 40)
    # ==========================================================
    # Chutamos Z(-40) = 0 e Z'(-40) = 1 (valores limpos para facilitar a álgebra)
    S0_teste = [0.0, 1.0] 
    
    sol_teste = solve_ivp(
                          Z_system,
                          [ll,lu],
                          S0_teste,
                          args=(F_inter,),
                          t_eval=L
                         )
    
    Z_teste_final = sol_teste.y[0, -1] # Valor obtido em l = 40
    
    # ==========================================================
    # PASSO 2: Aplicação Analítica da Invariância Linear
    # ==========================================================
    # Calculamos as constantes a e b para encaixar nos alvos exatos
    b = zl
    # pq esta errado, segundo minhas contas é zu+zl, mas assim esta ok? poq?
    a = (zu - zl) / Z_teste_final
    
    # Transformamos o perfil inteiro de uma vez só!
    Z_final = a * sol_teste.y[0] + b
    Zp_final = a * sol_teste.y[1] # Caso precise da derivada Z' também

    return Z_final,sol_teste.t 


# Definição do sistema de EDOs (Equação de Blasius)
def F_system(t, s):
    # s = [f, fp, fpp] -> f''' + f * f'' = 0 (ou variante)
    return [s[1], s[2], -s[0] * s[2]]

# Definição do sistema de EDOs (Equação de Blasius)
def Z_system(t,s,F_inter):

    Z, Zp = s
    F_atual = F_inter(t)
    # s = [Z, Zp] -> Z'' + F * Z'' = 0 (ou variante)
    return [s[1], -F_atual*s[1]]

def objective_shooting_F(fpp0_guess, up, ul,lu,ll):
    # 1. Chute inicial arbitrário para fp0 e o fpp0 controlado pelo fsolve
    fp0_arbitrary = 1.0  

    S0_dir = [0.0, fp0_arbitrary, fpp0_guess[0]]
    
    # Integra para a direita até l=40
    sol_dir = solve_ivp(F_system, [0, lu], S0_dir, t_eval=None)
    
    # Valor obtido no "infinito" à direita
    fp_inf_dir = sol_dir.y[1][-1]

    # Valor obtido no "meio" para fpp
    fpp_0_dir = sol_dir.y[2][0]
    
    # 2. Aplica a propriedade de invariância de escala (Transformação de Töpfer)
    # Para a equação de Blasius f''' + f * f'' = 0, a escala é: fp_real = c^2 * fp_arbitrary
    # Queremos que fp_real(lu) == up. Portanto: c = sqrt(up / fp_inf_dir)
    c = np.sqrt(up / fp_inf_dir)
    
    # Valores reais reescalados na origem (l=0) que garantem F'(lu) == up
    fp0_scaled = fp0_arbitrary * (c**2)
    fpp0_scaled = fpp_0_dir * (c**3)
    
    # 3. Integra para a esquerda até ll usando os valores corrigidos e conectados
    S0_esq = [0.0, fp0_scaled, fpp0_scaled]
    sol_esq = solve_ivp(F_system, [0, ll], S0_esq, t_eval=None)
    
    # Valor obtido na ponta esquerda pra fp
    fp_inf_esq = sol_esq.y[1][-1]
    
    # 4. O fsolve compara o resultado com a condição de contorno da esquerda (ul = 0.5)
    residual = fp_inf_esq - ul
    return [residual]

"""
def objective_shooting_F(fpp0_guess, up, ul,lu,ll):
    # 1. Chute inicial arbitrário para fp0 e o fpp0 controlado pelo fsolve
    fp0_arbitrary = 1.0  
    S0_dir = [0.0, fp0_arbitrary, fpp0_guess[0]]
    
    # Integra para a direita até l=40
    sol_dir = solve_ivp(F_system, [0, lu], S0_dir, t_eval=None)
    
    # Valor obtido no "infinito" à direita
    fp_inf_dir = sol_dir.y[1][-1]

    # Valor obtido no "meio" para fpp
    fpp_0_dir = sol_dir.y[2][0]
    
    # 2. Aplica a propriedade de invariância de escala (Transformação de Töpfer)
    # Para a equação de Blasius f''' + f * f'' = 0, a escala é: fp_real = c^2 * fp_arbitrary
    # Queremos que fp_real(lu) == up. Portanto: c = sqrt(up / fp_inf_dir)
    c = np.sqrt(up / fp_inf_dir)
    
    # Valores reais reescalados na origem (l=0) que garantem F'(lu) == up
    fp0_scaled = fp0_arbitrary * (c**2)
    fpp0_scaled = fpp_0_dir * (c**3)
    
    # 3. Integra para a esquerda até ll usando os valores corrigidos e conectados
    S0_esq = [0.0, fp0_scaled, fpp0_scaled]
    sol_esq = solve_ivp(F_system, [0, ll], S0_esq, t_eval=None)
    
    # Valor obtido na ponta esquerda pra fp
    fp_inf_esq = sol_esq.y[1][-1]
    
    # 4. O fsolve compara o resultado com a condição de contorno da esquerda (ul = 0.5)
    residual = fp_inf_esq - ul
    return [residual]
"""

