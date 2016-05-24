CheckAtm.lin = function(ATM = list()){
  # checks linear ATM variable to make sure it has all necessary info
  
  # base elevation, at which c0, wx0, wy0, and rho0 are valid
  if(is.null(ATM$z0)){
    ATM$z0 = 0
  }
  # base sound speed--330 m/s corresponds to 0 C
  if(is.null(ATM$c0)){
    ATM$c0 = 330
  }

  # sound speed gradient--use very small nonzero value
  if(is.null(ATM$gc)){
    ATM$gc = -10^-9
  }else if(ATM$gc == 0){
    ATM$gc = -10^-9
  }

  # zonal wind: assume zero, with zero gradient.
  if(is.null(ATM$wx0)){
    ATM$wx0 = 0
  }
  if(is.null(ATM$gwx)){
    ATM$gwx = 0
  }

  # meridional wind: assume zero, with zero gradient
  if(is.null(ATM$wy0)){
    ATM$wy0 = 0
  }
  if(is.null(ATM$gwy)){
    ATM$gwy = 0
  }

  # density: assume rho = 1 kg/m^3, zero gradient
  # equation: rho(z) = 1.013 * 10^5/(287.04 * 273) e^(-z/6800): Ps/(RT) e^(-z/H)
  # R: gas constant
  # T: temperature (K) (safe to assume constant--very small errors)
  # Ps: Surface pressure
  # H: scale length for atmosphere (m)
  # take a linear approximation to this equation
  
  if(is.null(ATM$rho0)){
    ATM$rho0 = 1.2929 * exp(-ATM$z0/6800)
  }
  if(is.null(ATM$grho)){
    ATM$grho = -1.013*10^5/(287.04*273*6800) * exp(-ATM$z0/6800)
  }

  return(ATM)
}
