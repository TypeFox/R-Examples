airtovac = function(wave_air) {
  
  wave_vac = wave_air
  g =which(wave_vac>=2000);     #Only modify above 2000 A
  
  for(iter in 0:1) {
    sigma2 = (1e4/wave_vac[g])^2     #Convert to wavenumber squared
    fact = 1 +  5.792105e-2/(238.0185 - sigma2) + 
      1.67917e-3/( 57.362 - sigma2)
    
    wave_vac[g] = wave_air[g]*fact              #Convert Wavelength
  }

  
  return(wave_vac)            
}
