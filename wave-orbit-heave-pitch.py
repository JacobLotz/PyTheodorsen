import numpy as np
import matplotlib.pyplot as plt

"""
Comparing the wave orbital to a combined heave and pitch motion.
"""



##############################################################################
# Wave orbital velocities
##############################################################################
# Requires scalars:
# - Wave frequency omega
# - Wave amplitude a
# - Wave number k
# - Water depth d
# - Depth of evaluation location z
# Time t or location x can be numpy arrays 
# but only one at a time to return a numpy array.

class Orbit:
   def __init__(self, omega_, a_, k_, d_, z_, t_, x_):
      self.omega = omega_ 
      self.a     = a_
      self.k     = k_
      self.d     = d_
      self.z     = z_
      self.t     = t_ 
      self.x     = x_ 

   def UAmplitude(self):
      # Shallow water
      # return self.omega * self.a * np.cosh(self.k*(self.d+self.z)) / np.sinh(self.k*self.d) 
      # assuming deep water
      return self.omega * self.a * np.exp(k*z)

   def VAmplitude(self):
      # Shallow water
      # return self.omega * self.a * np.sinh(self.k *(self.d+self.z)) / np.sinh(self.k*self.d) 
      # assuming deep water
      return self.omega * self.a * np.exp(k*z)

   def U(self):
      return self.UAmplitude() * np.sin(self.omega*self.t - self.k*self.x)

   def V(self):
      return self.VAmplitude()* np.cos(self.omega*self.t - self.k*self.x)


##############################################################################
# Hydrofoil heave motion induced velocities
##############################################################################
# Requires scalars:
# - Heave frequency omega
# - Heave amplitude a
# - Heave phase ph
# Time t or location x can be numpy arrays 
# but only one at a time to return a numpy array.

class Heave:
   def __init__(self, omega_, a_, ph_, t_, x_):
      self.omega = omega_ 
      self.a     = a_
      self.ph    = ph_ 
      self.t     = t_ 
      self.x     = x_ 

   def U(self):
      return 0;

   def V(self):
      return self.a * self.omega * np.cos(self.omega * self.t + self.ph)


##############################################################################
# Hydrofoil pitch motion induced velocities
##############################################################################
# Requires scalars:
# - Pitch frequency omega
# - Pitch amplitude a
# - Rotation fraction of the chord of the rotation point
# - Chord c
# - Pitch phase ph
# Time t or location x can be numpy arrays 
# but only one at a time to return a numpy array.

class Pitch:
   def __init__(self, omega_, a_, rf_, c_, ph_, t_, x_):
      self.omega = omega_ 
      self.a     = a_ 
      self.rf    = rf_ 
      self.c     = c_ 
      self.ph    = ph_
      self.t     = t_ 
      self.x     = x_ 

   def Angle(self):
      return self.a * np.sin(self.omega * self.t + self.ph);

   def AngularVelocity(self):
      return self.a * self.omega* np.cos(self.omega * self.t + self.ph);

   def UV(self):
      # We are neglecting the thickness of the foil section
      xrot = rf * c;
      return self.AngularVelocity() * (self.x - xrot)

   def U(self):
      return self.UV() * np.cos(self.Angle())
    
   def V(self):
      return self.UV() * np.sin(self.Angle())


##############################################################################
# Plotting of velocities
##############################################################################

# ----------------------------
# Input
# ----------------------------
# Time duration
T = 20 

# Precision
n = 300;

# Spatial details
x = 0.0;


# Hydrofoil details
c  = 1.014                   # Chord
rf = 1/2		                 # Rotation point
z  = -0.5                    # Submersion

# Wave
l      = c                   # Wave lenght
k      = 2*np.pi/l           # Wave number
womega = np.sqrt(k/9.81)     # Wave frequency [rad/s]: disp. rel. deep water
wa     = 0.3                 # Wave amplitude
d      = 10                  # Water depth: deep water

# Heave motion 1
homega1 = womega;            # Heave frequency
ha1     = 0.014 		        # Heave amplitude
hph1    = 0.0                # Heave phase

# Heave motion 2
homega2 = womega;            # Heave frequency
ha2     = 0.00               # Heave amplitude
hph2    = np.pi              # Heave phase

# Pitch motion
tomega = womega;             # Pitch frequency
ta     = 1.5			        # Pitch amplitude [degrees]
tph    = np.pi/2             # Pitch phase

# ----------------------------
# Initialize functions
# ----------------------------

# Compensate pitch phase
tph += np.pi * x


# Convert to radians
ta = ta/180*np.pi

# Time vector
t = np.linspace(0, T, n);

# Spatial vector
# Uncommend for spatial plot:
# t = 0
# x = np.linspace(0,c, n)

# Velocities
orbit  = Orbit(womega, wa, k, d, z, t, x)
heave1 = Heave(homega1, ha1, hph1, t, x)
heave2 = Heave(homega2, ha2, hph2, t, x)
pitch  = Pitch(tomega, ta, rf, c, tph, t, x)


# ----------------------------
# Initialize functions
# ----------------------------
wu = orbit.U()
wv = orbit.V()

foilu = heave1.U() + heave2.U() + pitch.U()
foilv = heave1.V() + heave2.V() + pitch.V()


# ----------------------------
# Plot functions
# ----------------------------
# Uncommend for spatial plot:
# t = x

plt.figure()
plt.subplot(211)
plt.plot(t, wu)
plt.plot(t, foilu)
plt.ylabel("u [m/s]")


plt.subplot(212)
plt.plot(t, wv, label='wave')
plt.plot(t, foilv, label='motion')
plt.xlabel("t [s]")
plt.ylabel("v [m/s]")
plt.legend()
plt.show()