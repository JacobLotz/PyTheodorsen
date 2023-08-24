import matplotlib.pyplot as plt
import numpy as np
import scipy.special as scsp


"""
Theodorsen function and resulting phase and ampl
"""

def theo_fun(k):
    r"""Returns the value of Theodorsen's function at a reduced frequency :math:`k`.

    .. math:: \mathcal{C}(jk) = \frac{H_1^{(2)}(k)}{H_1^{(2)}(k) + jH_0^{(2)}(k)}

    where :math:`H_0^{(2)}(k)` and :math:`H_1^{(2)}(k)` are Hankel functions of the second kind.

    Args:
        k (np.array): Reduced frequency/frequencies at which to evaluate the function.

    Returns:
        np.array: Value of Theodorsen's function evaluated at the desired reduced frequencies.

    from: https://www.programcreek.com/python/?code=ImperialCollegeLondon%2Fsharpy%2Fsharpy-master%2Fsharpy%2Futils%2Fanalytical.py
    """

    j = 1.0j
    H1 = scsp.hankel2(1, k)
    H0 = scsp.hankel2(0, k)

    C = H1 / (H1 + j * H0)

    return C

def TheoPhase(k, T):
    theo = theo_fun(k)
    theor = np.real(theo)
    theoi = np.imag(theo)
    return np.arctan(theoi/theor)/(2*np.pi)*T

def TheoAmpl(k):
    theo = theo_fun(k)
    theor = np.real(theo)
    theoi = np.imag(theo)
    return np.sqrt(theor**2 + theoi**2)


"""
Motion functions of the foil
"""
def HeaveMotion(omega, A, t):
	return A*np.sin(omega*t)

def HeaveVelocity(omega, A, t):
	return A*omega*np.cos(omega*t)

def HeaveAcceleration(omega, A, t):
	return -A*(omega**2)*np.sin(omega*t)


"""
Motion functions of the foil: easier to use
"""

def HeaveMotionSimple(T, A, n):
    t = np.linspace(0,T,n)
    omega = 2*np.pi/T
    return HeaveMotion(omega, A, t)

def HeaveVelocitySimple(T, A, n):
    t = np.linspace(0,T,n)
    omega = 2*np.pi/T
    return HeaveVelocity(omega, A, t)

def HeaveAccelerationSimple(T, A, n):
    t = np.linspace(0,T,n)
    omega = 2*np.pi/T
    return HeaveAcceleration(omega, A, t)



"""
Theodoorsen lift functions --> Easier: use the simple functions
"""
def CLAddedMass(omega, A, t):
	return np.pi *1.0/2.0 * HeaveAcceleration(omega, A, t)  

def CLCirc(omega, A, t, k):
    # The minus sign here is debatable. In other derivations it is not present. But I think it should be
    # out of phase with the heave velocity
	return -2.0 * np.pi * HeaveVelocity(omega, A, t)*TheoAmpl(k)

def CLTotal(omega, A, t, t_cir, k):
	return CLAddedMass(omega, A, t_cir) + CLCirc(omega, A, t, k)

"""
Simple functions that requires an easier input for for instance external scripts
"""
def CLAddedMassSimple(T, A, n):
    t = np.linspace(0,T,n)
    k = np.pi/T
    omega = 2*np.pi/T
    return CLAddedMass(omega, A, t)  

def CLCircSimple(T, A, n):
    t = np.linspace(0,T,n)
    k = np.pi/T
    t_cir = t + TheoPhase(k, T) # phase --> frame of reference thing
    omega = 2*np.pi/T	
    return CLCirc(omega, A, t_cir, k)



def CLTotalSimple(T, A, n):
    """
	Requires: 
	Heave period T
	Heave amplitude A
	Number of evaluation points n

	Returns total force
	"""
    t = np.linspace(0,T,n)
    k = np.pi/T
    t_cir = t + TheoPhase(k, T) # phase --> frame of reference thing
    omega = 2*np.pi/T	
    return CLAddedMass(omega, A, t) + CLCirc(omega, A, t_cir, k)


"""
QS Lift functions
"""
def AlphaSimple(T, A, n):
    v = HeaveVelocitySimple(T, A, n)
    return np.tan(-v)

def CLQTSimple(T, A, n):
    return 2*np.pi*AlphaSimple(T, A, n)

"""
Example
"""

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
plt.savefig('potential-theo/k'+str(k)+'.png', bbox_inches='tight')
plt.show();

"""
