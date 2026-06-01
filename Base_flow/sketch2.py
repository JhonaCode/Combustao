import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# =====================================================
# Velocity-profile sketch at x = 4
# =====================================================

r = 0.3

y_prof = np.linspace(-1.2, 1.2, 200)

U_prof = (
    0.5*(1+r)
    + 0.5*(1-r)*np.tanh(y_prof/0.35)
)

# Position of inset
x0 = 4.0      # location in main figure
y0 = 0.0

# Scaling
xscale = 1.0
yscale = 1.2

# Draw local axes
ax.plot(
    [x0, x0],
    [y0-yscale, y0+yscale],
    color='black',
    lw=1.5
)

ax.plot(
    [x0, x0+xscale],
    [y0-yscale, y0-yscale],
    color='black',
    lw=1.5
)

# Velocity profile
ax.plot(
    x0 + xscale*(U_prof-r)/(1-r),
    y0 + y_prof,
    color='black',
    lw=2.5
)

# Labels
ax.text(
    x0+0.3,
    y0+1.4,
    r'$U(y)$',
    fontsize=14
)

ax.text(
    x0-0.15,
    y0+1.25,
    r'$y$',
    fontsize=12
)

ax.text(
    x0+1.05,
    y0-1.35,
    r'$U$',
    fontsize=12
)

ax.text(
    x0+1.2,
    y0,
    r'$x=4$',
    fontsize=12
)
