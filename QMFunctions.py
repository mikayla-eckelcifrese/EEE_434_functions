from numpy import pi

#constants
h = 4.136 * 10**(-15) # in eV*s
h_bar = h/(2*pi) # in eV/Hz
m_e = 9.109 * 10**(-31) # in kg
c = 2.99 * 10**8 # in m/s
q = 1.6022*10**(-19) #in coloumbs
k_b = 8.617333262*10**(-5) #in eV/Kelvin

def momentum(E, mass = m_e):
  return (2*mass*E)**0.5

def wavelength(E, mass = m_e):
  return h/momentum(E, mass)

def energy(f):
  return h*f

def freq_light(wavelength_in_nanometers):
  return c/(wavelength_in_nanometers*10**(-9))