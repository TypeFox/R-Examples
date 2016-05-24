# Modified: 14 July 2015 SDH

watVap <- function(temp.k) {

  # Saturated vapor pressure of water
  #                                                                            A      B         C
  # This is based on NIST parameter values from Gubkov, Fermor, et al. 1964: 6.20963 2354.731 7.559
  # Oirignal (NIST) was in bar and K, so A in Pa (only change) is therefore: \Sexpr{6.20963 + log10(1E5)} = 6.20963 + 5 = 11.20963.
  
  if(any(temp.k < 273.15 | temp.k > 373.15)) 
    warning('in low level function WatVap(), temp.k is ', temp.k, ' K. Is this really correct?')
  
  return(10^(11.20963 - 2354.731/(temp.k + 7.559)))

}

