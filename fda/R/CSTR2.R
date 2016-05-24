CSTR2 <- function(Time, y, parms){
#function [Dy, DDy] = CSTR2(Time, y, parms)
#%  CSTR2 is called by the ODE solving function and evaluates the
#%  right hand side of the equation.
#
# 'lsoda' requires this function to return a list
# whose first component is a vector giving the derivatives of 'y'
#      with respect to 'Time' and 
# whose second component (if any) "is a vector ... of global values
#      that are required at each point in 'Time'".
# CSTR2 provides no second component.    

#%  R port 20 April 2007 by Spencer Graves
#%  Matlab version last modified 5 May 2005

  tau = 1;

#[F,CA0,T0,Tcin,Fc] = CSTR2in(Time, condition, tau);
  CSTR.in <- CSTR2in(Time, parms$condition, tau) 

  fitstruct <- parms$fitstruct
  kref   = fitstruct$kref;
  EoverR = fitstruct$EoverR;
  a      = fitstruct$a;
  b      = fitstruct$b;

  V      = fitstruct$V;
  FoverV = CSTR.in[, "F."]/V;

  aFc2bp = a*CSTR.in[, "Fc"]^b;
  const1 = V*fitstruct$rho*fitstruct$Cp;
  const2 = 2*fitstruct$rhoc*fitstruct$Cpc;

  betaTTcin = (CSTR.in[, "Fc"]*aFc2bp/(    
    const1*(CSTR.in[, "Fc"]+aFc2bp/const2)) )
  betaTT0   = FoverV;

#%  set up current values of variables
  {
    if( length(y) == 2 ){
      Ci = y[1];
      Ti = y[2];
    }
    else{
      Ci = y[,1];
      Ti = y[,2];
    }
  }
# end if(length(y)==2) 

  Tdif = 1/as.vector(Ti) - 1/fitstruct$Tref;
  betaCC = kref*exp(-1e4*EoverR*Tdif);
  betaTC = (-fitstruct$delH/const1)*betaCC;
  betaTT = FoverV + betaTTcin;

  DC = (-(FoverV + betaCC)*as.vector(Ci)  +
    (CSTR.in[, "F."]/V)*CSTR.in[, "CA0"]); 
  DT = (-betaTT*Ti + betaTC*Ci + (CSTR.in[, "F."]/V) * 
        CSTR.in[, "T0"] + betaTTcin*CSTR.in[, "Tcin"]) ;
#
  x <- cbind(Conc=DC, Temp=DT)
# 
  list(x)
}
