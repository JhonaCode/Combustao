import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# Definição do sistema de EDOs
def F_system(t, s):
    # s = [f, fp, fpp]
    return [s[1], s[2], -s[0] * s[2]]

def objective_shooting(guesses, up, ul):
    # O fsolve controla fp0 e fpp0 simultaneamente na origem (l=0)
    fp0_guess, fpp0_guess = guesses
    
    # --- PASSO 1: Chuta e integra para a DIREITA (0 até 40) ---
    S0_dir = [0.0, fp0_guess, fpp0_guess] # f(0)=0, f'(0)=fp0, f''(0)=fpp0
    sol_dir = solve_ivp(F_system, [0, 40], S0_dir, t_eval=None)
    fp_final_dir = sol_dir.y[1][-1] # F'(40)
    
    # --- PASSO 2: Usa as MESMAS condições iniciais e integra para a ESQUERDA (0 até -40) ---
    # Isso garante que a curva seja contínua e perfeitamente conectada no centro (l=0)
    S0_esq = [0.0, fp0_guess, fpp0_guess]
    sol_esq = solve_ivp(F_system, [0, -40], S0_esq, t_eval=None)
    fp_final_esq = sol_esq.y[1][-1] # F'(-40)
    
    # --- PASSO 3: Calcula os dois resíduos (erros) ---
    # Queremos que F'(40) chegue em up E que F'(-40) chegue em ul
    residual_dir = fp_final_dir - up
    residual_esq = fp_final_esq - ul
    
    # Retorna os dois erros. O fsolve vai alterar fp0 e fpp0 até ambos zerarem.
    return [residual_dir, residual_esq]

# --- Parâmetros do Problema ---
up = 1.0
ul = 0.5

# Chutes iniciais para as duas variáveis na origem: [F'(0), F''(0)]
chutes_iniciais = [0.6, 0.2]

# O fsolve ajusta fp0 e fpp0 juntos para satisfazer as duas pontas
valores_convergidos = fsolve(objective_shooting, chutes_iniciais, args=(up, ul))

fp0_final, fpp0_final = valores_convergidos
print("-> Convergência alcançada!")
print(f"-> F'(0)  encontrado: {fp0_final:.6f}")
print(f"-> F''(0) encontrado: {fpp0_final:.6f}")

# ==========================================
# RECONSTRUÇÃO DA SOLUÇÃO FINAL PARA PLOT
# ==========================================
l_eval_dir = np.linspace(0, 40, 200)
l_eval_esq = np.linspace(0, -40, 200)

S0_final = [0.0, fp0_final, fpp0_final]
sol_dir_plot = solve_ivp(F_system, [0, 40], S0_final, t_eval=l_eval_dir)
sol_esq_plot = solve_ivp(F_system, [0, -40], S0_final, t_eval=l_eval_esq)

# Unifica os dados da esquerda para a direita de forma ordenada
l_total = np.concatenate((sol_esq_plot.t[::-1], sol_dir_plot.t))
u_total = np.concatenate((sol_esq_plot.y[1][::-1], sol_dir_plot.y[1]))

# Gráfico explicativo
plt.figure(figsize=(10, 6))
plt.plot( u_total, l_total,label="F'(l) - Perfil de Velocidade", color='blue', lw=2.5)
#plt.axhline(y=up, color='r', linestyle='--', label=f'Alvo up = {up}')
#plt.axhline(y=ul, color='g', linestyle='--', label=f'Alvo ul = {ul}')
#plt.axvline(x=0, color='gray', linestyle=':', label='Origem l=0')
plt.ylabel('l')
plt.ylabel("F'(l)")
#plt.title("Shooting Method Tradicional (Sem Invariância) - Integração Partida")
plt.legend()
plt.grid(True)
plt.show()

