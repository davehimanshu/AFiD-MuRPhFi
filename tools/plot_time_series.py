"""
This script produces plots of the Nusselt number and
Reynolds number over time, using the simulation output
stored in `means.h5`.
"""

import matplotlib.pyplot as plt
import afidtools as afid
import numpy as np

#----------------------------------------------------------------------------------#
#------------------------- Set the simulation directory ---------------------------#
simdir = '/scratch/seismo/dave/rotate_run'
#----------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------#
#-------------- Load the time values, grid, and input parameters ------------------#
#----------------------------------------------------------------------------------#
t = afid.mean_time(simdir)
grid = afid.Grid(simdir)
inputs = afid.InputParams(simdir)
# Define the start and end times for averaging
t_start = 1500
t_end = 2000
# Find the indices corresponding to t_start and t_end
start_idx = np.searchsorted(t, t_start)
end_idx = np.searchsorted(t, t_end)
#----------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------#
# calculate mean values of Temperature and velocity and obtain temperature fluctuations #
#---------------------------------------------------------------------------------------#
Tbar = afid.read_mean(simdir, 'Tbar')
Trms = afid.read_mean(simdir, 'Trms')
Trms_fluc = np.sqrt(Trms**2 - Tbar**2)
vyrms = afid.read_mean(simdir, 'vyrms')
vzrms = afid.read_mean(simdir, 'vzrms')
vxrms = afid.read_mean(simdir, 'vxrms')
#---------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------#
#--------------------------------- Nusselt number ---------------------------------#
#----------------------------------------------------------------------------------#
# Calculate Nusselt numbers at lower and upper plates
Nu_pl = -(Tbar[0,:] - 0.5)/grid.xm[0]
Nu_pu = -(-0.5 - Tbar[-1,:])/(1.0 - grid.xm[-1])

# Define the Peclet number
Pec = (inputs.RayT*inputs.PraT)**0.5

# Calculate Nusselt number from scalar dissipation rate
chiT = afid.read_mean(simdir, 'chiT')
Nu_chi = afid.xmean(chiT, grid.xc)*Pec

# Calcuate Nusselt number from KE dissipation rate
#epsilon = afid.read_mean(simdir, 'epsilon')
#Nu_eps = 1.0 + afid.xmean(epsilon, grid.xc)*Pec

# Calculate Nusselt number from global turbulent heat flux
xT = afid.read_mean(simdir, 'vxT')
Nu_vol = 1.0 + afid.xmean(xT, grid.xc)*Pec

# Make the time series plot for Nusselt numbers
fig, ax = plt.subplots(figsize=(6.0,2.0), layout='constrained')
ax.plot(t, Nu_pl, label="$Nu_\mathrm{pl}$")
ax.plot(t, Nu_pu, label="$Nu_\mathrm{pu}$")
ax.plot(t, Nu_chi, label="$Nu_\chi$")
#ax.plot(t, Nu_eps, label="$Nu_\\varepsilon$")
ax.plot(t, Nu_vol, label="$Nu_\mathrm{vol}$")
ax.grid()
ax.legend(ncols=2)
ax.set(
    xlim=[0,t[-1]],
    ylim=[0,20],
    xlabel="$t/(H/U_f)$",
    ylabel='$Nu$'
)
fig.savefig('Nusselt.png')
#----------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------#
#--------------------------------- Averaged Nusselt number ------------------------#
#----------------------------------------------------------------------------------#
# Calculate the averaged Nusselt numbers
Nu_pl_avg = np.mean(Nu_pl[start_idx:end_idx])
Nu_pu_avg = np.mean(Nu_pu[start_idx:end_idx])
Nu_chi_avg = np.mean(Nu_chi[start_idx:end_idx])
Nu_vol_avg = np.mean(Nu_vol[start_idx:end_idx])

print(f"Averaged Nusselt numbers between t={t_start} and t={t_end}:")
print(f"Nu_pl_avg: {Nu_pl_avg}")
print(f"Nu_pu_avg: {Nu_pu_avg}")
print(f"Nu_chi_avg: {Nu_chi_avg}")
print(f"Nu_vol_avg: {Nu_vol_avg}")
#----------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------#
#--------------------------------- Reynolds number --------------------------------#
#----------------------------------------------------------------------------------#
# Define the dimensionless inverse of viscosity
inu = (inputs.RayT/inputs.PraT)**0.5

# Compute the Reynolds number based on vertical KE
Re_v = afid.xmean(vxrms**2, grid.xc)**0.5*inu

# Compute the Reynolds number based on horizontal KE
Re_h = afid.xmean(vyrms**2 + vzrms**2, grid.xc)**0.5*inu

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
    ylim=[0,150],
    xlabel='$t/(H/U_f)$',
    ylabel=r'$Re = \sqrt{2\mathcal{K}} H / \nu$'
)
fig.savefig('Reynolds.png')
#----------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------#
#--------------------------------- Averaged Reynolds number -----------------------#
#----------------------------------------------------------------------------------#
# Calculate the averaged Reynolds numbers
Re_v_avg = np.mean(Re_v[start_idx:end_idx])
Re_h_avg = np.mean(Re_h[start_idx:end_idx])
Re_t_avg = np.mean(Re_t[start_idx:end_idx])

print(f"Averaged Reynolds numbers between t={t_start} and t={t_end}:")
print(f"Re_v_avg: {Re_v_avg}")
print(f"Re_h_avg: {Re_h_avg}")
print(f"Re_t_avg: {Re_t_avg}")
#----------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------#
#------------------------Calculate average profiles over time ---------------------#
#----------------------------------------------------------------------------------#
# Calculate the time-averaged profiles for vyrms and vzrms
vyrms_avg = np.mean(vyrms[:, start_idx:end_idx], axis=1)
vzrms_avg = np.mean(vzrms[:, start_idx:end_idx], axis=1)
Trms_fluc_avg = np.mean(Trms_fluc[:, start_idx:end_idx], axis=1)

# Obtain the average horizontal profile from the above calculate profiles
vhrms_avg = (vyrms_avg + vzrms_avg)*0.5

# Plot the averaged profiles
fig3, ax3 = plt.subplots(figsize=(6.0,6.0), layout='constrained')
ax3.plot(grid.xm, vhrms_avg, label='vhrms_avg')
ax3.plot(grid.xm, Trms_fluc_avg, label='Trms_fluc_avg')
ax3.set(
    xlim=[0, 1],
    xlabel='xm',
    ylabel='Averaged horizontal velocity and temperature RMS fluctuations',
    title='Time-averaged velocity and temperature RMS fluctuation Profiles'
)
ax3.legend()
ax3.grid()
fig3.savefig('Velocity_Temp_RMS_Avg.png')

# Assume symmetry and calculate the average of the top and bottom halves
half_len = len(vyrms_avg) // 2
vyrms_avg_sym = (vyrms_avg[:half_len] + vyrms_avg[-1:-half_len-1:-1]) * 0.5
vzrms_avg_sym = (vzrms_avg[:half_len] + vzrms_avg[-1:-half_len-1:-1]) * 0.5
Temp_fluc_avg_sym = (Trms_fluc_avg[:half_len] + Trms_fluc_avg[-1:-half_len-1:-1]) * 0.5
# Update the horizontal profile accordingly
vhrms_avg_sym = (vyrms_avg_sym + vzrms_avg_sym) * 0.5

# Plot the averaged symmetrical profiles
fig4, ax4 = plt.subplots(figsize=(6.0,6.0), layout='constrained')
ax4.plot(grid.xm[:half_len], vhrms_avg_sym, label='vh_rms_avg_sym')
ax4.plot(grid.xm[:half_len], Temp_fluc_avg_sym, label='Temp_fluc_avg_sym')
ax3.set(
    xlim=[0, grid.xm[half_len-1]],
    xlabel='xm',
    ylabel='Velocity RMS',
    title='Time-averaged Symmetrical Velocity and temp RMS Profiles'
)
ax4.legend()
ax4.grid()
fig4.savefig('Velocity_Temp_RMS_Avg_Sym.png')
#----------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------#
#------------------------Calculate maxima using gradient ascent -------------------#
#----------------------------------------------------------------------------------#
def gradient_ascent(y, x, learning_rate=0.005, tolerance=1e-7, max_iter=2000):
   
    """
    Perform gradient ascent to find the maximum of a function given x and y arrays.

    Args:
        x (np.ndarray): Array of independent variable values.
        y (np.ndarray): Array of dependent variable values (function values).
        learning_rate (float): Step size for gradient ascent.
        tolerance (float): Convergence tolerance.
        max_iter (int): Maximum number of iterations.

    Returns:
        np.ndarray: Array containing the x and y values of the maximum [x_max, y_max].
    """
 
    # Ensure inputs are numpy arrays
    x = np.array(x)
    y = np.array(y)

    # Initial guess for x
    x_current = x[np.argmax(y)]  # Start near the highest observed point

    for _ in range(max_iter):
        # Estimate the gradient using finite differences
        idx = np.searchsorted(x, x_current)
        if idx <= 0 or idx >= len(x) - 1:
            break  # Exit if gradient can't be computed

        gradient = (y[idx + 1] - y[idx - 1]) / (x[idx + 1] - x[idx - 1])

        # Update x using the gradient
        x_next = x_current + learning_rate * gradient

        # Check for convergence
        if abs(x_next - x_current) < tolerance:
            break

        x_current = x_next

    # Find the corresponding y value
    y_max = np.interp(x_current, x, y)

    return np.array([x_current, y_max])

# Find maxima using gradient ascent
vhrms_max_idx = gradient_ascent(vhrms_avg_sym, grid.xm[:half_len])
print(f"The viscous boundary layer is: {vhrms_max_idx[0]}")
Temp_fluc_max_idx = gradient_ascent(Temp_fluc_avg_sym, grid.xm[:half_len])
print(f"The thermal boundary layer value is: {Temp_fluc_max_idx[0]}")

# Plot the maxima
fig5, ax5 = plt.subplots(figsize=(6.0,6.0), layout='constrained')
ax5.plot(grid.xm[:half_len], vhrms_avg_sym, label='vh_rms_avg_sym')
ax5.plot(grid.xm[:half_len], Temp_fluc_avg_sym, label='Temp_fluc_avg_sym')
ax5.scatter(vhrms_max_idx[0], vhrms_max_idx[1], color='red', label='vh_rms_max')
ax5.scatter(Temp_fluc_max_idx[0], Temp_fluc_max_idx[1], color='blue', label='Temp_fluc_max')
ax5.set(
    xlim=[0, grid.xm[half_len-1]],
    xlabel='xm',
    ylabel='Value',
    title='Maxima of Time-averaged Symmetrical Profiles'
)
ax5.legend()
ax5.grid()
fig5.savefig('Maxima_Profiles.png')
#----------------------------------------------------------------------------------#
