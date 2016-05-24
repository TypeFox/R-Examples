# Constants used in functions
# These are not exported
# Modified 7 April 2016 SDH

# Standard chemical formulas for mcomp
std.forms <- c(
	       # Standard macromolecules
	       vfa = 'C2H4O2', 
	       protein = 'C5H7O2N', 
	       carbohydrate = 'C6H10O5', 
	       lipid = 'C57H104O6', 
	       lignin = 'C10H13O3', 
	       # Cell biomass 
	       biomass = 'C5H7O2N', 
	       # Organic acids
	       acetic = 'CH3COOH', 
	       lactic = 'C3H6O3', 
	       # Alcohols
	       ethanol = 'CH3CH2OH', 
	       cellulose = 'C6H10O5', 
	       # Carbohydrates
	       glucose = 'C6H12O6')

# Standard atomic weights
# These are from CIAAW (http://www.ciaaw.org/atomic-weights.htm), significant digits only
atom.weights <- c(
  Al = 26.9815385,
  Sb = 121.760,
  Ar = 39.948,
  As = 74.921595,
  Ba = 137.328,
  Be = 9.0121831,
  Bi = 208.98040,
  B  = 10.8,
  Br = 79.9,
  Cd = 112.414,
  Cs = 132.90545197,
  Ca = 40.078,
  C  = 12.01,
  Ce = 140.116,
  Cl = 35.4,
  Cr = 51.9962,
  Co = 58.933194,
  Cu = 63.546,
  Dy = 162.500,
  Er = 167.259,
  Eu = 151.964,
  F  = 18.998403164,
  Gd = 157.25,
  Ga = 69.723,
  Ge = 72.631,
  Au = 196.966569,
  Hf = 178.49,
  He = 4.002602,
  Ho = 164.93033,
  H  = 1.008,
  In = 114.818,
  I  = 126.90447,
  Ir = 192.217,
  Fe = 55.845,
  Kr = 83.798,
  La = 138.90548,
  Pb = 207.2,
  Li = 7,
  Lu = 174.9668,
  Mg = 24.30,
  Mn = 54.938044,
  Hg = 200.592,
  Mo = 95.95,
  Nd = 144.242,
  Ne = 20.1798,
  Ni = 58.6934,
  Nb = 92.90637,
  N  = 14.007,
  Os = 190.23,
  O  = 16.00,
  Pd = 106.42,
  P  = 30.973761998,
  Pt = 195.085,
  K  = 39.0983,
  Pr = 140.90766,
  Pa = 231.03588,
  Re = 186.207,
  Rh = 102.90550,
  Rb = 85.4678,
  Ru = 101.07,
  Sm = 150.36,
  Sc = 44.955908,
  Se = 78.972,
  Si = 28.08,
  Ag = 107.8682,
  Na = 22.98976928,
  Sr = 87.62,
  S  = 32.1,
  Ta = 180.94788,
  Te = 127.60,
  Tb = 158.92535,
  Tl = 204.38,
  Th = 232.0377,
  Tm = 168.93422,
  Sn = 118.711,
  Ti = 47.867,
  W  = 183.84,
  U  = 238.02891,
  V  = 50.9415,
  Xe = 131.294,
  Yb = 173.054,
  Y  = 88.90584,
  Zn = 65.38,
  Zr = 91.224
  )

# Molar gas volumes (mL/mol, at 273.15 K and 101.325 kPa (1 atm))
# From NIST (search webbook, then select fluid properties, then isobaric, etc.)
vol.mol <- c(
  CH4 = 22360.588, 
  CO2 = 22263.009, 
  N2 = 22403.863, 
  H2 = 22427.978
  )
