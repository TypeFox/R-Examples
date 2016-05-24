intTD_NC <- function(func=function(x) return(x), T,
                     step, Nstep, intweightsfunc = intweights_CAV_SIM,
                      fail=NULL,
                      debug=FALSE,
                     ...){
  # compute numerical integral of func in [0 , T] following Newton_Cote method
  # func : (vector of) function to integrate, func(t, ...)
  # T    : upper bound (vector)
  # step : vector of the steps (one row per T)
  # Nstep : vector of the number of steps (T = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # fail : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-vector("numeric", length(T))
  for(i in 1:length(T)){
    # vector of evaluated t
    theT <- (0:Nstep[i])*step[i]
    # matrix of the evaluated functions (nt row, nfunc col)
    FF <- func(theT, i, ...)
    # weights 1 * nt matrix
    w<-intweightsfunc(Nstep[i])
# numerical integration
    res[i] <- crossprod(w , FF)
  }
  res * step 
}

intTD_NC_debug<- function(func=function(x) return(x), T, step, Nstep, intweightsfunc = intweights_CAV_SIM,
                     fail=NULL,
                      debug=FALSE, ...){
  # compute numerical integral of func in [0 , T] following Newton_Cote method
  # func : function to integrate, func(t, ...)
  # T    : upper bound (vector)
  # step : vector of the steps (one row per T)
  # Nstep : vector of the number of steps (T = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # fail : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-vector("numeric", length(T))
  for(i in 1:length(T)){
    # vector of evaluated t
    theT <- (0:Nstep[i])*step[i]
    # vector of the evaluated functions
    FF <- func(theT, i, ...)
    # weights
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i] <- crossprod(w , FF)
  }
  res * step 
}

intTD_SIM3_8 <- function(func=function(x) x, T, step, Nstep, ...){
  # compute numerical integral of func in [0 , T] following cavalieri Simpson method
  # func : function to integrate
  # T    : upper bound (vector)
  # step : vector of the steps (one lig per T)
  # Nstep : vector of the number of steps (T = Nstep * step), Nstep = 3 * k
  # weights are (1 3 3 2 3 3 2 3 3  ... 3 3 2 3 3 2 3 3 1 )* step * 3 / 8
  res<-rep(0, length(T))
           for(i in 1:length(T)){
             res[i] <- 3 * sum(func( (1:(Nstep[i]-1))*step, i, ...)) -
                           sum(func((3*(1:(Nstep[i]/3-1)))*step, i, ...))
                               
           }
           (res + func(0, i, ...) + func(T, i, ...) )* step * 3 / 8
}


intTD_BOOLE <- function(func=function(x) x, T, step, Nstep, ...){
  # compute numerical integral of func in [0 , T] following cavalieri Simpson method
  # func : function to integrate
  # T    : upper bound (vector)
  # step : vector of the steps (one lig per T)
  # Nstep : vector of the number of steps (T = Nstep * step), Nstep = 4 * k
  # weights are (7 32 12 32 14 32 12 32 14 ... 14 32 12 32 14 32 12 32 7 ) * 4 step / 90 
  res<-rep(0, length(T))
           for(i in 1:length(T)){
             res[i] <- 32 * sum(func((1:(Nstep[i]/2-1))*step, ...)) + 
                       12 * sum(func((4*(1:(Nstep[i]/4))-2)*step, ...)) +
                       14 * sum(func((4*(1:(Nstep[i]/4-1)))*step, ...))
           }
           (res + 7 * (func(0, ...) + func(T, ...)) )* step / 90
}

intTD_GL <- function(func=function(x) x, T, step, Nstep, ...){
  # compute numerical integral of func in [0 , T] following Gauss Legendre quadrature
}




intTD_GLM <- function(func=function(x) return(x), T, step, Nstep, intweightsfunc=NULL,
                      fail,
                      debug=FALSE, #Zalphabeta,
                     ...){
  # compute numerical integral of func in [0 , T] for equivalence with the poisson GLM trick
  # func : function to integrate, func(t, ...)
  # T    : upper bound (vector)
  # step : object of class GLMStepParam
  # Nstep : number of complete bands (< T)
  # T is in the Nstep'th band
  # intweightfunc function for computing weights : 
  #     weights are (b_i - b_(i-1))
  #                  with b0=0 and bn = T[j]
  # intweightsfunc=NULL, not used, for compatibility with ind_TD_base_NC
  # fail : statuts at T
  # ... : parameters of func()
  res<-vector("numeric", length(T))
  for(i in 1:length(T)){
    # vector of evaluated t
    if(Nstep[i]>0){
      if( fail[i] !=0 ){
        theT <- c(step@points[1:Nstep[i]] , T[i])
      }
      else {
        theT <- c(step@points[1:Nstep[i]] , (step@cuts[1+Nstep[i]]+T[i])/2)
      }
                                        # vector of the evaluated functions
      FF <- func(theT, i, ...)
                                        # weights
      w<- c(step@steps[1:Nstep[i]], T[i]-step@cuts[1+Nstep[i]])

     # numerical integration
      res[i] <- crossprod(w , FF)

    }
    else {
      if( fail[i] != 0 ){
        res[i] <- T[i]* func(T[i], i,  ...)
      }
      else {
        res[i] <- T[i]* func(T[i]/2, i, ...)
      }
    }
  }

 res  
}


