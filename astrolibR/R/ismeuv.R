ismeuv = function(wave,hcol,heicol=0.1*hcol,heiicol=0*hcol,fano=F) {

  ratio = wave/911.75
  tauh = wave*0.
  good =which(ratio<1); ngood=length( good )
  minexp = .Machine$double.min.exp #min exponent to avoid underflow

  r = ratio[good]
  z = sqrt( r/(1.0-r) )
  denom = rep(1.0, ngood)
  y = -2.*pi*z
  good1 =y>minexp;
  denom[good1] = (1.0 - exp(y[good1]))        
  tauh[good] = hcol * 3.44e-16 * (r^4)*exp(-4.0*z*atan(1/z)) /  denom
  
  tauheii = wave*0.
  ratio = 4. * wave/911.75
  good =ratio<1;

  r = ratio[good]
  z = sqrt( r/(1.0-r) )
  denom = rep(1.0, ngood)
  y = -2*pi*z
  good1 =which(y>minexp);
  denom[good1] = 1.0 - exp(y[good1])*4.
  tauheii[good] = heiicol * 3.44e-16 * (r^4)*exp(-4.0*z*atan(1/z)) / denom
  
  
  c1 = c(-2.953607e+01, 7.083061e+00, 8.678646e-01,-1.221932e+00,  
    4.052997e-02, 1.317109e-01, -3.265795e-02, 2.500933e-03 )
  c2 = c( -2.465188e+01, 4.354679e+00, -3.553024e+00, 5.573040e+00, 
    -5.872938e+00, 3.720797e+00, -1.226919e+00, 1.576657e-01 )
  q  = c(2.81, 2.51, 2.45, 2.44 )
  nu = c(1.610, 2.795, 3.817, 4.824 )
  fano_gamma = c(2.64061e-03, 6.20116e-04, 2.56061e-04, 1.320159e-04 )
  esubi = 3.0 - 1.0/nu^2 + 1.807317
  tauhei = wave*0.
  good =which( wave<503.97)

  x = log10(wave[good])
  y = x*0.
  good1 = (wave<46.0)
  y[good1] = polyidl( x[good1], c2)      
  good2 =wave>=46.0
  y[good2] = polyidl( x[good2], c1)
  if(fano){
    epsilon = 911.2671/wave
    for(i in 1:4) {       #loop over first four hei resonances
      x = 2.0 * ((epsilon-esubi[i] )/ fano_gamma[i] ) 
      y = y + log10( (x - q[i])^2/ (1 + x*x ) )
    }
  }
  
  tauhei[good] = heicol * 10^y
  
  return(tauh + tauheii + tauhei)
}
