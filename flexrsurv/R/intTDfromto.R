intTDft_NC <- function(func=function(x) return(x), fromT, toT,
                     step, Nstep, intweightsfunc = intweights_CAV_SIM,
                      fail=NULL,
                      debug=FALSE,
                     ...){
  # compute numerical integral of func in [fromT , toT] following Newton_Cote method
  # func : (vector of) function to integrate, func(t, ...)
  # fromT : lower bound (vector)
  # toT   : upper bound (vector)
  # step : vector of the steps (one row per toT)
  # Nstep : vector of the number of steps ((toT - fromT) = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # fail : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-vector("numeric", length(toT))
  for(i in 1:length(toT)){
    # vector of evaluated t
    theT <- fromT[i] + (0:Nstep[i])*step[i]
    # matrix of the evaluated functions (nt row, nfunc col)
    FF <- func(theT, i, fromT=fromT, toT=toT, ...)
    # weights 1 * nt matrix
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i] <- crossprod(w , FF)
  }
  res * step 
}

intTDft_NC_debug<- function(func=function(x) return(x),  fromT, toT,
                           step, Nstep, intweightsfunc = intweights_CAV_SIM,
                           fail=NULL,
                           debug=FALSE, ...){
  # compute numerical integral of func in [fromT , toT] following Newton_Cote method
  # func : function to integrate, func(t, ...)
  # fromT : lower bound (vector)
  # toT   : upper bound (vector)
  # step : vector of the steps (one row per T)
  # Nstep : vector of the number of steps ((toT - fromT) = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # fail : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-vector("numeric", length(toT))
  for(i in 1:length(toT)){
    # vector of evaluated t
    theT <- fromT[i] + (0:Nstep[i])*step[i]
    # vector of the evaluated functions
    FF <- func(theT, i, fromT=fromT, toT=toT, ...)
    # weights
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i] <- crossprod(w , FF)
  }
  res * step 
}

intTDft_SIM3_8 <- function(func=function(x) x, fromT, toT,
                          step, Nstep, ...){
  # compute numerical integral of func in [fromT , toT] following cavalieri Simpson method
  # func : function to integrate
  # fromT : lower bound (vector)
  # toT   : upper bound (vector)
  # step : vector of the steps (one lig per T)
  # Nstep : vector of the number of steps ((toT - fromT) = Nstep * step), Nstep = 3 * k
  # weights are (1 3 3 2 3 3 2 3 3  ... 3 3 2 3 3 2 3 3 1 )* step * 3 / 8
  res<-rep(0, length(toT))
           for(i in 1:length(toT)){
             res[i] <- 3 * sum(func( fromT + (1:(Nstep[i]-1))*step, i, fromT=fromT, toT=toT, ...)) -
                           sum(func(fromT + (3*(1:(Nstep[i]/3-1)))*step, i, fromT=fromT, toT=toT, ...))
                               
           }
           (res + func(fromT, i, fromT=fromT, toT=toT, ...) + func(toT, i, fromT=fromT, toT=toT, ...) )* step * 3 / 8
}


intTDft_BOOLE <- function(func=function(x) x, fromT, toT,
                          step, Nstep, ...){
  # compute numerical integral of func in [fromT , toT] following cavalieri Simpson method
  # func : function to integrate
  # fromT : lower bound (vector)
  # toT   : upper bound (vector)
  # step : vector of the steps (one ligne per T)
  # Nstep : vector of the number of steps ((toT - fromT) = Nstep * step), Nstep = 4 * k
  # weights are (7 32 12 32 14 32 12 32 14 ... 14 32 12 32 14 32 12 32 7 ) * 4 step / 90 
  res<-rep(0, length(toT))
           for(i in 1:length(toT)){
             res[i] <- 32 * sum(func(fromT + (1:(Nstep[i]/2-1))*step, fromT=fromT, toT=toT, ...)) + 
                       12 * sum(func(fromT + (4*(1:(Nstep[i]/4))-2)*step, fromT=fromT, toT=toT, ...)) +
                       14 * sum(func(fromT + (4*(1:(Nstep[i]/4-1)))*step, fromT=fromT, toT=toT, ...))
           }
           (res + 7 * (func(fromT, fromT=fromT, toT=toT, ...) + func(toT, fromT=fromT, toT=toT, ...)) )* step / 90
}

intTDft_GL <- function(func=function(x) x, fromT, toT,
                      step, Nstep, ...){
  # compute numerical integral of func in [0 , T] following Gauss Legendre quadrature
}




intTDft_GLM <- function(func=function(x) return(x), fromT, toT,
                       step,
                       Nstep, 
                       intweightsfunc=NULL,
                      fail,
                      debug=FALSE, #Zalphabeta,
                     ...){
  # compute numerical integral of func in [0 , T] for equivalence with the poisson GLM trick
  # func : function to integrate, func(t, ...)
  # fromT : lower bound (vector)
  # toT   : upper bound (vector)
  # step : object of class GLMStepParam
  # Nstep : index of the first and last complete band ( fromT[i] < step@cuts[Nstep[i,1]] <= step@cuts[Nstep[i,2]+1] < toT)
  #                                                    ( fromT[i] < step@points[Nstep[i,1]] <= step@points[Nstep[i,2]] < toT)
  # Nstep  : index of the first complete band 
  # toT is in the Nstep'th band
  # fromT is in the first band
  # intweightfunc function for computing weights : 
  #     weights are (b_i - b_(i-1))
  #                  with b0=0 and bn = T[j]
  # intweightsfunc=NULL, not used, for compatibility with ind_TD_base_NC
  # fail : statuts at toT
  # ... : parameters of func()
  res<-vector("numeric", length(toT))
  for(i in 1:length(toT)){
    # vector of evaluated t
    if(Nstep[i,2]>= Nstep[i,1]){
      # at least one complete step
      if( fail[i] ){
        theT <- c((step@cuts[Nstep[i,1]]+fromT[i])/2,
                  step@points[Nstep[i,1]:Nstep[i,2]] ,
                  toT[i])
      }
      else {
        theT <- c((step@cuts[Nstep[i,1]]+fromT[i])/2,
                  step@points[Nstep[i,1]:Nstep[i,2]] ,
                  (step@cuts[1+Nstep[i,2]]+toT[i])/2)
      }
                                        # vector of the evaluated functions
      FF <- func(theT, i, fromT=fromT, toT=toT, ...)
                                        # weights
      w<- c(step@cuts[Nstep[i,1]] - fromT[i],
            step@steps[Nstep[i,1]:Nstep[i,2]],
            toT[i]-step@cuts[1+Nstep[i,2]])

     # numerical integration
      res[i] <- crossprod(w , FF)
    }
    else if((Nstep[i,2] - Nstep[i,1]) == -1L){
# fromT and toT are in 2 successive bands
      # Nstep[i,2] + 1 =  Nstep[i,1]
      if( fail[i] ){
        theT <- c((step@cuts[Nstep[i,1]]+fromT[i])/2,
                  toT[i])
      }
      else {
        theT <- c((step@cuts[Nstep[i,1]]+fromT[i])/2,
                  (step@cuts[Nstep[i,1]]+toT[i])/2)
      }
                                        # vector of the evaluated functions
      FF <- func(theT, i, fromT=fromT, toT=toT, ...)
                                        # weights
      w<- c(step@cuts[Nstep[i,1]] - fromT[i],
            toT[i]-step@cuts[Nstep[i,1]])

     # numerical integration
      res[i] <- crossprod(w , FF)
    }
    else {   #if((Nstep[i,2] - Nstep[i,1]) == -2L){
# fromT and toT are in the same band
      if( fail[i] ){
        res[i] <- (toT[i] - fromT[i]) * func(toT[i], i, fromT=fromT, toT=toT,  ...)
      }
      else {
        res[i] <- (toT[i] - fromT[i]) * func((toT[i] + fromT[i])/2, i, fromT=fromT, toT=toT,  ...)
      }
    }
  }

  
res  
}




