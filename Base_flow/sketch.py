import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import plotparameters as pp 
size_wg = 0.5
size_hf = 1.0
cmas    = 0.0

plotdef = 'diurnal2' 
    
tama    = pp.plotsize(size_wg,size_hf, cmas,plotdef)

# ==========================================================
# Parameters
# ==========================================================
xf = 2.5                 # flame starts downstream of x=0
xmax = 12.0

# ==========================================================
# Figure
# ==========================================================
fig, ax = plt.subplots(figsize=(10,4))

# ==========================================================
# Mixing-layer thickness
# ==========================================================
x = np.linspace(0, xmax, 800)

delta = 0.08 + 0.30*np.sqrt(x)

# ==========================================================
# Separator plate
# ==========================================================
plate = np.array([
    [-4.0,  0.12],
    [-0.25, 0.12],
    [ 0.0,  0.00],   # trailing edge at x=0
    [-0.25,-0.12],
    [-4.0,-0.12]
])

ax.add_patch(
    Polygon(
        plate,
        closed=True,
        facecolor='dimgray',
        edgecolor='black',
        linewidth=2,
        zorder=30
    )
)

# ==========================================================
# Upstream streams (x < 0)
# ==========================================================
ax.fill_between(
    [-4, 0],
    [3, 3],
    [0.12, 0.12],
    color='goldenrod',
    alpha=0.45,
    zorder=1
)

ax.fill_between(
    [-4, 0],
    [-0.12, -0.12],
    [-3, -3],
    color='lightskyblue',
    alpha=0.40,
    zorder=1
)

# ==========================================================
# Downstream streams (x > 0)
# ==========================================================
ax.fill_between(
    x,
    delta,
    3,
    color='goldenrod',
    alpha=0.45,
    zorder=1
)

ax.fill_between(
    x,
    -3,
    -delta,
    color='lightskyblue',
    alpha=0.40,
    zorder=1
)

# ==========================================================
# Mixing region before flame
# ==========================================================
mask_mix = x < xf

ax.fill_between(
    x[mask_mix],
    -delta[mask_mix],
    delta[mask_mix],
    color='khaki',
    alpha=0.55,
    zorder=2
)

# ==========================================================
# Reacting region
# ==========================================================
mask_react = x >= xf

ax.fill_between(
    x[mask_react],
    -delta[mask_react],
    delta[mask_react],
    color='orange',
    alpha=0.18,
    zorder=2
)

# ==========================================================
# Mixing-layer boundaries
# ==========================================================
ax.plot(
    x,
    delta,
    '--',
    color='black',
    linewidth=2,
    dashes=(6,4)
)

ax.plot(
    x,
    -delta,
    '--',
    color='black',
    linewidth=2,
    dashes=(6,4)
)

# ==========================================================
# Diffusion flame
# ==========================================================
xflame = np.linspace(xf, xmax, 400)

flame_thickness = (
    0.03
    + 0.025*np.sqrt(xflame-xf)
)

# outer glow
ax.fill_between(
    xflame,
    -4*flame_thickness,
    4*flame_thickness,
    color='gold',
    alpha=0.30,
    zorder=10
)

# intermediate glow
ax.fill_between(
    xflame,
    -2*flame_thickness,
    2*flame_thickness,
    color='orange',
    alpha=0.55,
    zorder=11
)

# flame core
ax.fill_between(
    xflame,
    -flame_thickness,
    flame_thickness,
    color='red',
    alpha=0.95,
    zorder=12
)

# ==========================================================
# Flow arrows (upper stream)
# ==========================================================
for xpos in [-4.0, -2.0, -0.6]:
    ax.arrow(
        xpos, 1.6,
        0.8, 0,
        head_width=0.08,
        head_length=0.18,
        fc='darkmagenta',
        ec='darkmagenta',
        linewidth=1.5
    )

# ==========================================================
# Flow arrows (lower stream)
# ==========================================================
for xpos in [-4.0, -2.0, -0.6]:
    ax.arrow(
        xpos, -1.6,
        0.8, 0,
        head_width=0.08,
        head_length=0.18,
        fc='darkblue',
        ec='darkblue',
        linewidth=1.5
    )


# ==========================================================
# Labels
# ==========================================================
ax.text(
    -3.6,
    2.25,
    "Fuel stream"+"\n"+r'$\mathrm{u_1^{*},Y_{F1}^{*},T_1^{*},\rho_1^{*}, M_{w1}^{*}}$',
    fontsize=12,
    color='darkmagenta'
)

ax.text(
    -3.45,
    -2.65,
    "Oxidizer stream"+"\n"+r'$\mathrm{u_2^{*},Y_{F2}^{*},T_2^{*},\rho_2^{*}, M_{w2}^{*}}$',
    fontsize=12,
    color='darkblue'
)

ax.annotate(
    'Diffusive flame',
    xy=(8,0),
    xytext=(8.5,0.45),
    color='red',
    fontsize=12,
    arrowprops=dict(
        arrowstyle='->',
        color='red',
        lw=2
    )
)

#ax.annotate(
#    'Mixing layer',
#    xy=(2.5,delta[-1]),
#    xytext=(2.3,2.2),
#    fontsize=12
#)

ax.annotate(
    'Separating plate',
    xy=(-0.05,0),
    xytext=(-3.1,-0.75),
    fontsize=12,
    arrowprops=dict(
        arrowstyle='->',
        lw=1.8
    )
)

# x=0 location
ax.axvline(
    0,
    color='k',
    linestyle='--',
    linewidth=1.5,
    dashes=(6,6)
)

ax.text(
     0.20,
    -0.12,
    r'$(\mathrm{x^{*},y^{*}=(0,0)}$',
    fontsize=10
)
# =====================================================
# Velocity-profile sketch at x = 4
# =====================================================

r = 0.3

y_prof = np.linspace(-5.2, 5.2, 200)

U_prof = (
    0.5*(1+r)
    + 0.5*(1-r)*np.tanh(y_prof/0.35)
)

# Position of inset
x0 = 1.5      # location in main figure
y0 = 0.0

# Scaling
xscale = 1.0
yscale = 1.2


# Velocity profile
ax.plot(
    x0 + xscale*(U_prof-r)/(1-r),
    y0 + y_prof,
    color='black',
    lw=2.5
)


# ---------------------------------------
# delta_99 thickness
# ---------------------------------------
y99 = 0.35*np.arctanh(0.98)

ax.annotate(
    '',
    xy=(x0+1.25, y99),
    xytext=(x0+1.25, -y99),
    arrowprops=dict(
        arrowstyle='<->',
        lw=2,
        color='black'
    )
)

ax.text(
    x0+1.30,
    0.28,
    r'$\delta^{*}$',
    fontsize=14,
    va='center'
)

# Optional guide lines
ax.plot(
    [x0+1.15, x0+1.35],
    [ y99, y99],
    'k',
    lw=1
)

ax.plot(
    [x0+1.15, x0+1.35],
    [-y99,-y99],
    'k',
    lw=1
)


# ==========================================================
# Axes formatting
# ==========================================================
ax.set_xlim(-4, xmax)
ax.set_ylim(-3, 3)

ax.set_xlabel(r'$\mathrm{x}^*$', fontsize=18)
ax.set_ylabel(r'$\mathrm{y}^*$', fontsize=18)

ax.set_xticks([])
ax.set_yticks([])

ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)

plt.title("Compressible Reacting Mixing Layer with a Diffusive Flame", fontsize=14)

plt.tight_layout()

plt.savefig('%s/%s.pdf'%(pp.out_fig,'skecht'),bbox_inches='tight',dpi=200, format='pdf')

plt.show()
