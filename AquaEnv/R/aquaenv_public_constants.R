# a list of physical-chemical constants used in AquaEnv
PhysChemConst <- list(
                  R         = 83.14472,        # (bar*cm3)/(mol*K) the gas constant (corrected after Lewis1998, in Millero1995: R = 83.131; digits extended after DOE1994, Dickson2007)
                  F         = 96485.3399,      # C/mol the Faraday constant (charge per mol of electrons) (N_A*e-) (Dickson2007)
                  uMolToMol = 1e-6,            # conversion factor from umol to mol
                  absZero   = -273.15,         # absolute zero in degrees centigrade
                  e         = 79,              # relative dielectric constanf of seawater (Zeebe2001)
                  
                  K_HNO2    = 1.584893e-3,     # dissociation constant of HNO2:  mol/l, NBS pH scale, hybrid constant (Riordan2005)
                  K_HNO3    = 23.44,           # dissociation constant of HNO3:  assumed on mol/kg-soln and free pH scale, stoichiometric constant (Soetaert pers. comm.)
                  K_H2SO4   = 100,             # dissociation constant of H2SO4: assumed on mol/kg-soln and free pH scale, stoichiometric constant (Atkins1996)
                  K_HS      = 1.1e-12          # dissociation constant of HHS:   assumed on mol/kg-soln and free pH scale, stoichiometric constant (Atkins1996)
                  )

# a list of technical constants used in AquaEnv
Technicals <- list(
                   Haccur   = 1e-12,               # accuracy for iterative (Follows2006) pH calculations (max. deviation in [H+])
                   Hstart   = 1e-8,                # start [H+] for an iterative pH calculation
                   maxiter  = 100,                 # maximum number of iterations for iterative (Follows2006) pH calculation method as well as for the application of the standard R function uniroot
                                               
                   unirootinterval = c(1e-18, 1),  # the interval (in terms of [H+]) for pH calculation using the standard R function uniroot
                   uniroottol      = 1e-20,        # the interval (in terms of [H+]) for pH calculation using the standard R function uniroot

                   epsilon_fraction = 0.1,         # fraction of disturbance for the numerical calculation of derivatives of TA with respect to changes in the dissociation constants
                   revelle_fraction = 1e-5         # fraction of disturbance for the numerical calculation of the revelle factor
                   )
               

# coefficients for the pressure correction of dissociation constants and solubility products (Millero1995 WITH CORRECTIONS BY Lewis1998 (CO2Sys)!!!!!!)
DeltaPcoeffs <- data.frame(  
                           K_HSO4        = c(-18.03, 0.0466, 0.3160e-3,- 4.53, 0.0900,0),
                           K_HF          = c( -9.78,-0.0090,-0.9420e-3,- 3.91, 0.0540,0),
                           K_CO2         = c(-25.50, 0.1271, 0.0000e-3,- 3.08, 0.0877,0),
                           K_HCO3        = c(-15.82,-0.0219, 0.0000e-3,  1.13,-0.1475,0),
                           K_W           = c(-25.60, 0.2324,-3.6246e-3,- 5.13, 0.0794,0),
                           K_BOH3        = c(-29.48, 0.1622, 2.6080e-3,- 2.84, 0.0000,0),
                           K_NH4         = c(-26.43, 0.0889,-0.9050e-3,- 5.03, 0.0814,0),
                           K_H2S         = c(-14.80, 0.0020,-0.4000e-3,  2.89, 0.0540,0),
                           K_H3PO4       = c(-14.51, 0.1211,-0.3210e-3,- 2.67, 0.0427,0),
                           K_H2PO4       = c(-23.12, 0.1758,-2.6470e-3,- 5.15, 0.0900,0),
                           K_HPO4        = c(-26.57, 0.2020,-3.0420e-3,- 4.08, 0.0714,0),
                           K_SiOH4       = c(-29.48, 0.1622, 2.6080e-3,- 2.84, 0.0000,0), # same as K_BOH3
                           K_SiOOH3      = c(-29.48, 0.1622, 2.6080e-3,- 2.84, 0.0000,0), # same as K_BOH3
                           Ksp_calcite   = c(-48.76, 0.5304, 0.0000e-3,-11.76, 0.3692,0), 
                           Ksp_aragonite = c(-45.96, 0.5304, 0.0000e-3,-11.76, 0.3692,0)
                           )


# mean molecular mass of key chemical species in seawater in g/mol (DOE1994, Dickson2007)
MeanMolecularMass <- data.frame(
                                  Cl  = 35.453,
                                  SO4 = (32.065+4*(15.999)),
                                  Br  = 79.904,
                                  F   = 18.998,
                                  Na  = 22.990,
                                  Mg  = 24.3050,
                                  Ca  = 40.078,
                                  K   = 39.098,
                                  Sr  = 87.62,
                                  B   = 10.811
                                  )


# concentrations of key chemical species in seawater, relative with respect to chlorinity (DOE1994)
ConcRelCl <- data.frame(
                        Cl  = 0.99889,
                        SO4 = 0.1400,
                        Br  = 0.003473,
                        F   = 0.000067,
                        Na  = 0.55661,
                        Mg  = 0.06626,
                        Ca  = 0.02127,
                        K   = 0.0206,
                        Sr  = 0.00041,
                        B   = 0.000232
                        )
 

