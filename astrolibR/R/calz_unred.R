calz_unred = function( wave, flux, ebv, R_V=4.05) {
 
 w1 =which((wave>=6300) & (wave<=22000)); c1=length( w1 )
 w2 =which((wave>= 912) & (wave< 6300)); c2=length( w2 )
 x  = 10000.0/wave                      #Wavelength in inverse microns
 if( (c1 + c2)!=length(wave))
   warning('some elements of wavelength vector outside valid domain')

 klam = 0.0*flux
 klam[w1] = 2.659*(-1.857 + 1.040*x[w1]) + R_V
   
 klam[w2] = 2.659*(polyidl(x[w2], c(-2.156, 1.509e0, -0.198e0, 0.011e0))) + R_V
 
 funred = flux*10.0^(0.4*klam*ebv)
 return(funred)
}
