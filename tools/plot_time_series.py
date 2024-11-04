"""
This script produces plots of the Nusselt number and
Reynolds number over time, using the simulation output
stored in `means.h5`.
"""

import matplotlib.pyplot as plt
import afidtools as afid

simdir = '/scratch/seismo/dave/test_yang/2D'

# Load the time values, grid, and input parameters
t = afid.mean_time(simdir)
grid = afid.Grid(simdir)
inputs = afid.InputParams(simdir)

## Nusselt number

Tbar = afid.read_mean(simdir, 'Tbar')
print(Tbar)
# Calculate Nusselt numbers at lower and upper plates
Nu_pl = (0.5 + Tbar[0,:])/grid.xm[0]
print(Nu_pl)
Nu_pu = (0.5 - Tbar[-1,:])/(inputs.alx3 - grid.xm[-1])
print(Nu_pu)

# Define the Peclet number
Pec = (inputs.RayT*inputs.PraT)**0.5

# Calculate Nusselt number from scalar dissipation rate
chiT = afid.read_mean(simdir, 'chiT')
Nu_chi = afid.xmean(chiT, grid.xc)*Pec

# Calcuate Nusselt number from KE dissipation rate
epsilon = afid.read_mean(simdir, 'epsilon')
Nu_eps = 1.0 + afid.xmean(epsilon, grid.xc)*Pec

# Calculate Nusselt number from global turbulent heat flux
wT = afid.read_mean(simdir, 'vxT')
Nu_vol = 1.0 - afid.xmean(wT, grid.xc)*Pec

fig, ax = plt.subplots(figsize=(6.0,2.0), layout='constrained')
#ax.plot(t, Nu_pl, label="$Nu_\mathrm{pl}$")
#ax.plot(t, Nu_pu, label="$Nu_\mathrm{pu}$")
ax.plot(t, Nu_chi, label="$Nu_\chi$")
ax.plot(t, Nu_eps, label="$Nu_\\varepsilon$")
#ax.plot(t, Nu_vol, label="$Nu_\mathrm{vol}$")
ax.grid()
ax.legend(ncols=2)
ax.set(
    xlim=[0,t[-1]],
    ylim=[0,100],
    xlabel="$t/(H/U_f)$",
    ylabel='$Nu$'
)
fig.savefig('Nusselt.png')

"""
## Reynolds number

# Define the dimensionless inverse of viscosity
inu = (inputs.RayS/inputs.PraS)**0.5

# Compute the Reynolds number based on vertical KE
wrms = afid.read_mean(simdir, 'vxrms')
Re_v = afid.xmean(wrms**2, grid.xc)**0.5*inu

# Compute the Reynolds number based on horizontal KE
vrms = afid.read_mean(simdir, 'vzrms')
urms = afid.read_mean(simdir, 'vyrms')
Re_h = afid.xmean(urms**2 + vrms**2, grid.xc)**0.5*inu

# Make the time series plot
fig, ax = plt.subplots(figsize=(6.0,2.0), layout='constrained')
ax.plot(t, Re_v, label='vertical')
ax.plot(t, Re_h, label='horizontal')
Re_t = (Re_h**2 + Re_v**2)**0.5
ax.plot(t, Re_t, label='total')
ax.grid()
ax.legend()
ax.set(
    xlim=[0,t[-1]],
    ylim=[0,800],
    xlabel='$t/(H/U_f)$',
    ylabel=r'$Re = \sqrt{2\mathcal{K}} H / \nu$'
)
fig.savefig('Reynolds.svg')
"""