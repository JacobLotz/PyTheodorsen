from theodorsen import *
import matplotlib.pyplot as plt

"""
Example
"""
# Input
k = 2.0
extrusion = np.pi/k
A = 0.1;               # Heave amplitude
n = 50				   # Number of time steps

# extrusion = 10*np.pi   # Period T

# Calc var
t = np.linspace(0,extrusion, n)
k = np.pi/extrusion

# Plot
fig, axs = plt.subplots(1)
axs.set_ylabel(r'$C_L$ / heave motion' )
axs.set_xlabel('t')


plt.plot(t, CLAddedMassSimple(extrusion, A, n), color='red', linestyle=':', linewidth = 1.0, marker='', markersize='3', label = r"$C_l$ AM")
plt.plot(t, CLCircSimple(extrusion, A, n), color='green', linestyle=':', linewidth = 1.0, marker='', markersize='3',  label = r"$C_l$ circ")
plt.plot(t, CLTotalSimple(extrusion, A, n), label = r"$C_l$ total")
plt.plot(t, HeaveMotionSimple(extrusion, A, n),  linewidth = 1, color = 'black',  label = r"$y$")
plt.plot(t, HeaveVelocitySimple(extrusion, A, n),  linewidth = 1, color = 'black',  linestyle ="--", label = r"$\dot{y}$")
plt.plot(t, HeaveAccelerationSimple(extrusion, A, n),  linewidth = 1, color = 'black',  linestyle =":", label = r"$\ddot{y}$")


# Shrink current axis by 20%
box = axs.get_position()
axs.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
axs.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.title("k = " + str(k) + ' -- theo phase = ' + "{:.3f}".format(TheoPhase(k, extrusion)) + ' [s] -- theo fact = ' + "{:.3f}".format(TheoAmpl(k)), fontsize=10)
plt.savefig('k'+str(k)+'.png', bbox_inches='tight')
plt.show();