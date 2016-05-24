intTD_base_base_NC <- function(func=function(x) return(x), T, Spline, step, Nstep, intweightsfunc = intweights_CAV_SIM,
                     ...){
  # compute numerical integral of func*base_i(t)*base_j(t)  in [0 , T] following Newton_Cote method
  # func : (vector of) function to integrate, func(t, ...)
  # T    : upper bound (vector)
  # Spline : Spline parameters
  # step : vector of the lags (one row per T)
  # Nstep : vector of the number of lags (T = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # ... : parameters of func()
#  cat("inintTD_NC\n")
  res<-array(0, dim=c(length(T), Spline@nbase, Spline@nbase))
  for(i in 1:length(T)){
    # vector of evaluated t
    theT <- (0:Nstep[i])*step[i]
    # matrix of the evaluated functions (nt row, nfunc col)
    FF <- func(theT, i, ...)
    TBase <- fevaluate(Spline, theT, intercept=TRUE)
    # weights 1 * nt matrix
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i,,] <- crossprod(TBase, (w*FF)*TBase )
  }
#  cat("outinintTD_NC\n")
  res * step 
}

intTD_base_base_NC_debug<- function(func=function(x) return(x), T, Spline, step, Nstep, intweightsfunc = intweights_CAV_SIM,
                     ...){
  # compute numerical integral of func*base_i(t))*base_j(t) in [0 , T] following Newton_Cote method
  # func : function to integrate, func(t, ...)
  # T    : upper bound (vector)
  # Spline : Spline parameters
  # step : vector of the lags (one row per T)
  # Nstep : vector of the number of lags (T = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # ... : parameters of func()
  cat("inintTD_NC\n")
  cat("lengthT lengthNstep step \n")
  cat(length(T), length(Nstep), length(step))
  cat("\n")
  print(cbind(T,Nstep, step)[1:20,])
  cat("\n")
  res<-array(0, dim=c(length(T), Spline@nbase, Spline@nbase))
  for(i in 1:length(T)){
    # vector of evaluated t
    theT <- (0:Nstep[i])*step[i]
    TBase <- fevaluate(Spline, theT, intercept=TRUE)
    # vector of the evaluated functions
    FF <- func(theT, i, ...)
    # weights
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i,,] <- crossprod(TBase, (w*FF)*TBase )
  }
  cat("outinintTD_NC\n")
  res * step 
}


intTD_base_base_GLM <- function(func=function(x) return(x), T, Spline, step, Nstep, intweightsfunc=NULL,
                     ...){
  # compute numerical integral of func*b_i(t)*base_j(t) in [0 , T] for equivalence with the poisson GLM trick
  # func : function to integrate, func(t, ...)
  # T    : upper bound (vector)
  # Spline : Spline parameters
  # step : object of class GLMLagParam
  # Nstep : number of complete bands (< T)
  # T is in the (Nstep+1)'th band
  # ... : parameters of func()
#  cat("inintTD_NC\n")
  res<-array(0, dim=c(length(T), Spline@nbase, Spline@nbase))
  for(i in 1:length(T)){
    # vector of evaluated t
    if(Nstep[i]>0){
      theT <- c(step@points[1:Nstep[i]] , T[i])
                                        # vector of the evaluated functions
      FF <- func(theT, i, ...)
      # matrix of bases evaluated at theT
      TBase <- fevaluate(Spline, theT, intercept=TRUE)
                                        # weights
      w<- c(step@steps[1:Nstep[i]], T[i]-step@bands[1+Nstep[i]])
#      cat("i length T F W nlag\n")
#      cat(c(i, length(theT), length(FF), length(w), Nstep[i]))
#      cat("\n")
#      if( Nstep[i]==1){
#        cat("w\n")
#        cat(w)
#        cat("\n")
#      }
                                        # numerical integration
      res[i,] <- crossprod( TBase, (w*FF) * TBase)
    }
    else {
# Nstep[i]==0
      res[i,] <- crossprod(fevaluate(Spline, T, intercept=TRUE),
                           T[i]* func(T[i], i, ...)%*% fevaluate(Spline, T, intercept=TRUE))
    }
  }
  res  
}

fastintTD_base_base_GLM <- function(func=function(x) return(x), T, Spline, step, Nstep, intweightsfunc=NULL,
                     ...){
  # compute numerical integral of func*b_i(t)*base_j(t) in [0 , T] for equivalence with the poisson GLM trick
  # func : function to integrate, func(t, ...)
  # T    : upper bound (vector)
  # Spline : Spline parameters
  # step : object of class GLMStepParam
  # Nstep : number of complete bands (< T)
  # T is in the (Nstep+1)'th band
  # ... : parameters of func()
#  cat("inintTD_NC\n")
  res<-array(0, dim=c(length(T), Spline@nbase, Spline@nbase))
                                        # matrix of bases evaluated at theT
  allTBase <- fevaluate(Spline,step@points[1:Nstep[i]] , intercept=TRUE)
  TBaseatT <- fevaluate(Spline,T , intercept=TRUE)
  for(i in 1:length(T)){
    # vector of evaluated t
    if(Nstep[i]>0){
      theT <- c(step@points[1:Nstep[i]] , T[i])
                                        # vector of the evaluated functions
      FF <- func(theT, i, ...)
                                        # weights
      w<- c(step@steps[1:Nstep[i]],  T[i]-step@bands[1+Nstep[i]] )
#      cat("i length T F W nstep\n")
#      cat(c(i, length(theT), length(FF), length(w), Nstep[i]))
#      cat("\n")
#      if( Nstep[i]==1){
#        cat("w\n")
#        cat(w)
#        cat("\n")
#      }
                                        # numerical integration
      res[i,] <- crossprod(allTBase[1:Nstep[i],], (w*FF)*allTBase[1:Nstep[i],] ) +
         crossprod(TBaseatT[i,], (T[i]-step@bands[1+Nstep[i]])*func(T[i], i, ...) *  TBaseatT[i,])
    }
    else {
# Nstep[i]==0
      res[i,] <- crossprod(TBaseatT[i,], T[i]* func(T[i], i, ...) * TBaseatT[i,])
    }
  }
  res  
}


