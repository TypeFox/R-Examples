#' Fiducial cosmology parameter list
#' 
#' @name parameter.fidcosmo
#' @docType data
#' @description the parameter list is a crucial input for other calculations
parameter.fidcosmo <-
  list(omegaM0=0.3,
       omegaL0=0.7,
       omegaK=0.0,
       h=0.7)

#' Physical constants in SI units, source: [PDG]
#' @name parameter.constant
#' @docType data
#' @description Useful physical constants in SI units.
parameter.constant <- 
  list( c = 299792458.,            # speed of light 
        me = 9.10938291e-31,       # electron mass
        h = 6.62606957e-34,        # Planck constant
        mp = 1.672621777e-27,      # proton mass
        e = 1.602176565e-19,       # electron charge
        G = 6.67384e-11,           # gravitational constant
        k = 1.3806488e-23,         # Blotzmann constant
        sigma.T = 6.652458734e-29  # Thomson scatter cross-section 
  )

#' useful conversion of astrophysical units
#' @docType data
parameter.unit =
  list( Mpc = 3.0856776e22, # Megaparsecs [Wiki]
        Msun = 1.98892e30,  # Mass of Sun [Wiki]
        eV = 1.602176565e-19 # electron Volt [PDG]
  )

#' Hubble Distance
#' \eqn{D_H = c/H_0} [Mpc/h]
#' @docType data
parameter.DH = 2997.92458