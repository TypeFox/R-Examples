readline("Hit Return\n")

##############   set up a velocity model

### layer 1:
alpha1 = 4.98
  beta1 =  2.9
   rho1 = 2.667

##########  layer 2
  alpha2 = 8.0
   beta2 = 4.6
   rho2 = 3.38

##########  visualize the scattering coefficients

  App =  pzoeppritz( "Amplitude" , alpha1, alpha2, beta1, beta2, rho1 ,rho2, "P", "ALL");

readline("Hit Return\n")


###########  change incoming wave to S-wave:
App =  pzoeppritz( "Amplitude" , alpha1, alpha2, beta1, beta2, rho1 ,rho2, "S", "ALL");
#####################################################
#############  Incident wave in high velocity layer
  alpha1 = 8.0
  beta1 =  4.6
   rho1 = 3.38

  alpha2 = 4.98
   beta2 = 2.9
   rho2 = 2.667


readline("Hit Return\n")

 App =  pzoeppritz( "Amplitude" , alpha1, alpha2, beta1, beta2, rho1 ,rho2, "P", "ALL");

readline("Hit Return\n")

###########  change incoming wave to S-wave:
App =  pzoeppritz( "Amplitude" , alpha1, alpha2, beta1, beta2, rho1 ,rho2, "S", "ALL");

