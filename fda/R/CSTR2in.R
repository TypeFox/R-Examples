CSTR2in <- function(Time, condition = 
   c('all.cool.step', 'all.hot.step', 'all.hot.ramp', 'all.cool.ramp',
     'Tc.hot.exponential', 'Tc.cool.exponential', 'Tc.hot.ramp',
     'Tc.cool.ramp', 'Tc.hot.step', 'Tc.cool.step'),
   tau=1){
##
## Simulated Input Vectors for
## Continuously Stirred Temperature Reactor:
##
## Returns data.frame(Fvec, CA0vec, T0vec, Tcinvec, Fcvec)

##
## condition = c( 'all_cool_step', 'all_hot_step',
##                'all_hot_ramp', 'all_cool_ramp',
##     'Tc_hot_exponential', 'Tc_cool_exponential',
##                 'Tc_hot_ramp', 'Tc_cool_ramp'
##                 'Tc_hot_step', 'Tc_cool_step'

#  Last modified 20 April 2007 by Spencer Graves
#%  Matlab version last modified 2 May 2005
  rtnMat <- function(Fvec=Fvec, CA0vec=CA0vec,
        T0vec=T0vec, Tcinvec=Tcinvec, Fcvec=Fcvec){
    x <- cbind(F.=Fvec, CA0=CA0vec, T0=T0vec,
               Tcin=Tcinvec, Fc=Fcvec)
#    Fnames <- dimnames(Fvec)[[1]]
#    dimnames(x) <- list(Fnames,
#        c("Fvec", "CA0vec", "T0vec", "Tcinvec", "Fcvec"))
    x
  }
#
  n       = length(Time);
# defaults used for condition == Tc* 
  CA0vec  =   rep(2, n) 
  T0vec   = rep(323, n) 
  Fcvec   =  rep(15, n);
  
#switch condition
  if(condition[1] == 'all.cool.step'){
#        %  compute F
        
    Fvec = rep(1.0, n);
    Fvec[(4 <= Time) & (Time <  8)] = 1.5;
    Fvec[(8 <= Time) & (Time < 12)] = 0.5;
        
#        %  compute C_{A0}
        
    CA0vec = rep(2.0, n);
    CA0vec[(16 <= Time) & (Time < 20)] = 2.2;
    CA0vec[(20 <= Time) & (Time < 24)] = 1.8;
        
#        %  compute T0
        
    T0vec = rep(323,n);
    T0vec[(28 <= Time) & (Time < 32)] = 343;
    T0vec[(32 <= Time) & (Time < 36)] = 303;
        
#        %  compute Tcin
        
    Tcinbase = 335;
    Tcinvec  = rep(Tcinbase,n);
    Tcinvec[(40 <= Time) & (Time < 44)] = Tcinbase+5;
    Tcinvec[(44 <= Time) & (Time < 48)] = Tcinbase-5;
        
#        %  compute Fc
        
    Fcvec = rep(15,n);
    Fcvec[(52 <= Time) & (Time < 56)] = 20;
    Fcvec[(56 <= Time) & (Time < 60)] = 10;

    return(rtnMat(Fvec, CA0vec, T0vec, Tcinvec, Fcvec))
  }
#
  if(condition[1] == 'all.hot.step'){
        
#        %  compute F
        
    Fvec = rep(1.0,n);
    Fvec[(4 <= Time) & (Time <  8)] = 1.5;
    Fvec[(8 <= Time) & (Time < 12)] = 0.5;
        
#        %  compute C_{A0}
        
    CA0vec = rep(2.0, n);
    CA0vec[(16 <= Time) & (Time < 20)] = 2.2;
    CA0vec[(20 <= Time) & (Time < 24)] = 1.8;
        
#        %  compute T0
        
    T0vec = rep(323, n);
    T0vec[(28 <= Time) & (Time < 32)] = 343;
    T0vec[(32 <= Time) & (Time < 36)] = 303;
        
#        %  compute Tcin
        
    Tcinbase = 365;
    Tcinvec  = rep(Tcinbase, n);
    Tcinvec[(40 <= Time) & (Time < 44)] = Tcinbase+5;
    Tcinvec[(44 <= Time) & (Time < 48)] = Tcinbase-5;
        
#        %  compute Fc
        
    Fcvec = rep(15, n);
    Fcvec[(52 <= Time) & (Time < 56)] = 20;
    Fcvec[(56 <= Time) & (Time < 60)] = 10;

    return(rtnMat(Fvec, CA0vec, T0vec, Tcinvec, Fcvec))

  }

  if(condition[1] == 'all.hot.ramp'){

    Fvec    = rep(1.0, n);
    Tcinvec = rep(365, n);
        
    index = ((2 <= Time) & (Time < 10));
    Fvec[index] = 1.5;
    index = ((Time >= 10) & (Time < 14));
    Fvec[index] = -0.2*(Time[index]-10) + 1.5;
    index = ((Time >= 14) & (Time < 18));
    Fvec[index] = -0.2*4 + 1.5;
    
    index = ((Time >= 26) & (Time < 34));
    CA0vec[index] = 2.5;
    index = ((Time >= 34) & (Time < 38));
    CA0vec[index] = -0.2*(Time[index]-34) + 2.5;
    index = ((Time >= 38) & (Time < 42));
    CA0vec[index] = -0.2*4 + 2.5;
    
    index = ((Time >= 50) & (Time < 58));
    T0vec[index] = 353;
    index = ((Time >= 58) & (Time < 62));
    T0vec[index] = -20*(Time[index]-58) + 353;
    index = ((Time >= 62) & (Time < 66));
    T0vec[index] = -20*4 + 353;
        
    index = ((Time >= 74) & (Time < 82));
    Tctop = 390;
    Tcinvec[index] = Tctop;
    index = ((Time >= 82) & (Time < 86));
    Tcinvec[index] = -10*(Time[index]-82) + Tctop;
    index = ((Time >= 86) & (Time < 90));
    Tcinvec[index] = -10*4 + Tctop;
        
    index = ((Time >=  98) & (Time < 106));
    Fcvec[index] = 25;
    index = ((Time >= 106) & (Time < 110));
    Fcvec[index] = -5*(Time[index]-106) + 25;
    index = ((Time >= 110) & (Time < 114));
    Fcvec[index] = -5*4 + 25;

    return(rtnMat(Fvec, CA0vec, T0vec, Tcinvec, Fcvec))

  }

  if(condition[1] == 'all.cool.ramp'){

    Fvec    = rep(0.05, n);
    Tcinvec = rep(330, n);
        
    index = (( 2 <= Time) & (Time < 10));
    Fvec[index] = 1.5;
    index = ((Time >= 10) & (Time < 14));
    Fvec[index] = -0.2*(Time[index]-10) + 1.5;
    index = ((Time >= 14) & (Time < 18));
    Fvec[index] = -0.2*4 + 1.5;
    
    index = ((Time >= 26) & (Time < 34));
    CA0vec[index] = 2.5;
    index = ((Time >= 34) & (Time < 38));
    CA0vec[index] = -0.2*(Time[index]-34) + 2.5;
    index = ((Time >= 38) & (Time < 42));
    CA0vec[index] = -0.2*4 + 2.5;
    
    index = ((Time >= 50) & (Time < 58));
    T0vec[index] = 353;
    index = ((Time >= 58) & (Time < 62));
    T0vec[index] = -20*(Time[index]-58) + 353;
    index = ((Time >= 62) & (Time < 66));
    T0vec[index] = -20*4 + 353;
    
    index = ((Time >= 74) & (Time < 82));
    Tctop = 355;
    Tcinvec[index] = Tctop;
    index = ((Time >= 82) & (Time < 86));
    Tcinvec[index] = -10*(Time[index]-82) + Tctop;
    index = ((Time >= 86) & (Time < 90));
    Tcinvec[index] = -10*4 + Tctop;
    
    index = ((Time >=  98) & (Time < 106));
    Fcvec[index] = 25;
    index = ((Time >= 106) & (Time < 110));
    Fcvec[index] = -5*(Time[index]-106) + 25;
    index = ((Time >= 110) & (Time < 114));
    Fcvec[index] = -5*4 + 25;

    return(rtnMat(Fvec, CA0vec, T0vec, Tcinvec, Fcvec))

  }

  if(condition[1] == 'Tc.hot.exponential'){
        
    Fvec    = rep(1.0, n);
    Tcinvec = rep(365, n);
    
    index = (Time < 10); 
    Tcinvec[index] = 400 - (400 - 365)*exp(-(Time[index])/tau);
    index = ((10 <= Time) & (Time < 20)); 
    Tcinvec[index] = 344 - (344 - 400)*exp(-(Time[index] - 10)/tau);
    index = (20 <= Time); 
    Tcinvec[index] = 365 - (365 - 344)*exp(-(Time[index] - 20)/tau);

    return(rtnMat(Fvec, CA0vec, T0vec, Tcinvec, Fcvec))

  }

  if(condition[1] == 'Tc.cool.exponential'){
        
    Fvec    = rep(0.05, n);
    Tcinvec =  rep(330, n);
    
    index = (Time/5 < 10); 
    Tcinvec[index] = 400 - (400 - 365)*exp(-(Time[index] )/tau);
    index = ((10 <= Time/5) & (Time/5 < 20)); 
    Tcinvec[index] = 344 - (344 - 400)*exp(-(Time[index] - 10)/tau);
    index = (20 <= Time/5); 
    Tcinvec[index] = 365 - (365 - 344)*exp(-(Time[index] - 20)/tau);

    return(rtnMat(Fvec, CA0vec, T0vec, Tcinvec, Fcvec))
    
  }

  if(condition[1] == 'Tc.hot.ramp'){

    Fvec    = rep(1.0, n);
    Tcinvec = rep(365, n);
        
    index = ((2 <= Time) & (Time < 10)); 
    Tcinvec[index] = 400;
    index = ((10 <= Time) & (Time < 14)); 
    Tcinvec[index] = -14*(Time[index]-10) + 400;
    index = ((14 <= Time) & (Time < 18)); 
    Tcinvec[index] = -14*4 + 400;

    return(rtnMat(Fvec, CA0vec, T0vec, Tcinvec, Fcvec))
    
  }

  if(condition[1] == 'Tc.cool.ramp'){

    Fvec    = rep(0.05, n);
    Tcinvec =  rep(330, n);
        
    index = ((2 <= Time/5) & (Time/5 < 10)); 
    Tcinvec[index] = 340;
    index = ((10 <= Time/5) & (Time/5 < 14)); 
    Tcinvec[index] = -0.8*(Time[index]-5*10) + 340;
    index = ((14 <= Time/5) & (Time/5 < 22)); 
    Tcinvec[index] = -5*0.8*4 + 340;

    return(rtnMat(Fvec, CA0vec, T0vec, Tcinvec, Fcvec))
    
  }

  if(condition[1] == 'Tc.hot.step'){

    Fvec    = rep(1.0, n);
    Tcinvec = rep(365, n);
    
    index = ((2 <= Time & Time < 12)); 
    Tcinvec[index] = 350;

    return(rtnMat(Fvec, CA0vec, T0vec, Tcinvec, Fcvec))
    
  }

  if(condition[1] == 'Tc.cool.step'){

    Fvec    = rep(0.05, n);
    Tcinvec =  rep(335, n);
    
    index = ((2 <= Time) & (Time < 12)); 
    Tcinvec[index] = 320;
        
    return(rtnMat(Fvec, CA0vec, T0vec, Tcinvec, Fcvec))
    
  }
}

