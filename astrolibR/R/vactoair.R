vactoair = function(wave_vac) {
  wave_air = wave_vac
  g = wave_vac>=2000     #Only modify above 2000 A
  
    sigma2 = (1e4/wave_vac[g])^2.   #Convert to wavenumber squared
    fact = 1 +  5.792105e-2/(238.0185 - sigma2) + 
      1.67917e-3/( 57.362 - sigma2)
    
    wave_air[g] = wave_vac[g]/fact

  return(wave_air)
}                        
