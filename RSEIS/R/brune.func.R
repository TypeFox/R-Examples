`brune.func` <-
function( freq,  omega0,  tstar0, fc, alpha, gamma)
{

  SMALL.NUMBER = 1e-300
  gam2 = 2.0*gamma;
   a1 = freq^(-alpha);
   a2 = (freq/fc)^(gam2);

   tstar = tstar0* a1 ;
   e2ft = exp((-pi*freq*tstar));


   tmod = (omega0*e2ft) / sqrt(1+a2 );

  ## plot(freq, tmod, type='l', log='xy')
  
   tmod[tmod<=0.0] = SMALL.NUMBER ;
  

   return(tmod);

}

