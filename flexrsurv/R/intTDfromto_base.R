intTDft_base_NC <- function(func=function(x) return(x), fromT, toT,
                            Spline,
                            step, Nstep, intweightsfunc = intweights_CAV_SIM,
                            fail=NULL,
                            debug=TRUE,
                            ...){
  # compute numerical integral of func*base_i(t)  in [fromT , toT] following Newton_Cote method
  # func : (vector of) function to integrate, func(t, ...)
  # fromT : lower bound (vector)
  # toT   : upper bound (vector)
  # Spline : Spline parameters
  # step : vector of the steps (one row per T)
  # Nstep : vector of the number of steps ((toT - fromT) = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # fail : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-matrix(0, nrow = length(toT), ncol = Spline@nbases + Spline@log)
  for(i in 1:length(toT)){
    # vector of evaluated t
    theT <- fromT[i] + (0:Nstep[i])*step[i]
    # matrix of the evaluated functions (nt row, nfunc col)
    FF <- func(theT, i, ...)
    TBase <- fevaluate(Spline, theT, intercept=TRUE)
    # weights 1 * nt matrix
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i,] <- crossprod(w*FF, TBase)
  }
#  cat("outinintTD_NC\n")
  res * step 
}

intTDft_base2_NC <- function(func=function(x) return(x), fromT, toT,
                            Spline,
                            step, Nstep, intweightsfunc = intweights_CAV_SIM,
                            fail=NULL,
                            debug=TRUE,
                            ...){
  #similar to intTDft_base_NC but
  # compute numerical integral of func*base_i(t - fromT)  in [fromT , toT] following Newton_Cote method
  # func : (vector of) function to integrate, func(t, ...)
  # fromT : lower bound (vector)
  # toT   : upper bound (vector)
  # Spline : Spline parameters
  # step : vector of the steps (one row per T)
  # Nstep : vector of the number of steps ((toT - fromT) = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # fail : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  res<-matrix(0, nrow = length(toT), ncol = Spline@nbases + Spline@log)
  for(i in 1:length(toT)){
    # vector of evaluated t
    theT <- fromT[i] + (0:Nstep[i])*step[i]
    # evaluate spline basis at t - fromT
    TBase <- fevaluate(Spline, (0:Nstep[i])*step[i], intercept=TRUE)
    # matrix of the evaluated functions (nt row, nfunc col)
    FF <- func(theT, i, ...)
    # weights 1 * nt matrix
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i,] <- crossprod(w*FF, TBase)
  }
#  cat("outinintTD_NC\n")
  res * step 
}

intTDft_base_NC_debug<- function(func=function(x) return(x),  fromT, toT,
                                 Spline,
                                 step, Nstep, intweightsfunc = intweights_CAV_SIM,
                                 fail=NULL,
                                 debug=TRUE,
                                 ...){
  # compute numerical integral of func*base_i(t)  in [fromT , toT] following Newton_Cote method
  # func : (vector of) function to integrate, func(t, ...)
  # fromT : lower bound (vector)
  # toT   : upper bound (vector)
  # Spline : Spline parameters
  # step : vector of the steps (one row per T)
  # Nstep : vector of the number of steps ((toT - fromT) = Nstep * step), Nstep is even
  # intweightfunc function for computing weights : 
  #     - NC-2 : Cavalieri-Simpson method intweight_CAV_SIM(), Nstep is even
  #     - NC-3 : Simpson 3/8   intweight_SIM_3_8(), Nstep = 3*1
  #     - NC-4 : Boole intweight_BOOLE(), Nstep = 4 I
  # fail : unused but present for compatibility with inTD_GLM
  # ... : parameters of func()
  cat("inintTD_NC_debug\n")
  cat("lengthT lengthNstep step \n")
  cat(length(T), length(Nstep), length(step))
  cat("\n")
  print(cbind(T,Nstep, step)[1:20,])
  cat("\n")
  res<-matrix(0, nrow = length(toT), ncol = Spline@nbases + Spline@log)
  for(i in 1:length(toT)){
    # vector of evaluated t
    theT <- fromT[i] + (0:Nstep[i])*step[i]
    TBase <- fevaluate(Spline, theT, intercept=TRUE)
    # vector of the evaluated functions
    FF <- func(theT, i, ...)
    # weights
    w<-intweightsfunc(Nstep[i])

# numerical integration
    res[i,] <- crossprod(w*FF, TBase)
  }
  cat("outinintTD_NC\n")
  res * step 
}



fastintTDft_base_GLM <- function(func=function(x) return(x), fromT, toT,
                                 Spline,
                               step, Nstep, intweightsfunc=NULL,
                               fail,
                               debug=FALSE,
                               ...){
  # compute numerical integral of func*b_i(t) in [fromT , toT] for equivalence with the poisson GLM trick
  # func : function to integrate, func(t, ...)
  # fromT : lower bound (vector)
  # toT   : upper bound (vector)
  # Spline : Spline parameters
  # step : object of class GLMStepParam
  # Nstep : index of the first and last complete band ( fromT[i] < step@cuts[Nstep[i,1]] <= step@cuts[Nstep[i,2]+1] < toT)
  #                                                    ( fromT[i] < step@points[Nstep[i,1]] <= step@points[Nstep[i,2]] < toT)
  # intweightsfunc=NULL, not used, for compatibility with ind_TD_base_NC
  # fail : statuts at T
  # ... : parameters of func()
  if(debug>200) {
    cat("fastinintTD_base_glm\n")
  }
  res<-matrix(0, nrow = length(fromT), ncol = Spline@nbases+Spline@log)
                     # matrix of bases evaluated at the points and T
  allTBase <- fevaluate(Spline,step@points , intercept=TRUE)
  Tpoints <- ifelse(fail, toT,  (step@cuts[1+Nstep[,2]]+toT)/2) 
  TBaseattoT <- fevaluate(Spline, Tpoints , intercept=TRUE)
  TBaseatfromT <- fevaluate(Spline, (step@cuts[Nstep[,1]]+fromT)/2 , intercept=TRUE)
  for(i in 1:length(toT)){
    # vector of evaluated t
    if(Nstep[i,2]>= Nstep[i,1]){
      # at least one complete step
        theT <- c((step@cuts[Nstep[i,1]]+fromT[i])/2,
                  step@points[Nstep[i,1]:Nstep[i,2]] ,
                  Tpoints[i] )
                                        # vector of the evaluated functions
        FF <- func(theT, i, ...)
                                        # weights
      w<- c(step@cuts[Nstep[i,1]] - fromT[i],
            step@steps[Nstep[i,1]:Nstep[i,2]],
            toT[i]-step@cuts[1+Nstep[i,2]])
                                        # numerical integration of the complete bands
        res[i,] <- crossprod(w*FF,
                             rbind( TBaseatfromT[i,],
                                   allTBase[Nstep[i,1]:Nstep[i,2],, drop=FALSE],
                                   TBaseattoT[i,]))
      }
    else if(Nstep[i,2] - Nstep[i,1] == -1L){
# fromT and toT are in 2 successive bands
      # Nstep[i,2] + 1 =  Nstep[i,1]
        theT <- c((step@cuts[Nstep[i,1]]+fromT[i])/2,
                  Tpoints[i] )
                                        # vector of the evaluated functions
        FF <- func(theT, i, ...)
                                        # weights
      w<- c(step@cuts[Nstep[i,1]] - fromT[i],
            toT[i]-step@cuts[Nstep[i,1]])
                                        # numerical integration of the complete bands
        res[i,] <- crossprod(w*FF,
                             rbind( TBaseatfromT[i,], TBaseattoT[i,]))
      }
    else { #if((Nstep[i,2] - Nstep[i,1]) == -2L){
# fromT and toT are in the same band
      if( fail[i]!=0 ){
        res[i,] <- (toT[i]- fromT[i]) * func(toT[i], i, ...) * TBaseattoT[i,]
      }
      else {
        res[i,] <-  ((toT[i] - fromT[i]) * func((toT[i] + fromT[i])/2, i,  ...)) %*% fevaluate(Spline, (toT[i] + fromT[i])/2 , intercept=TRUE) #[,,drop=TRUE] 
      }
    }
  }
  res  
}



fastintTDft_base2_GLM <- function(func=function(x) return(x), fromT, toT,
                                 Spline,
                               step, Nstep, intweightsfunc=NULL,
                               fail,
                               debug=FALSE,
                               ...){
  #similar to fastintTDft_base2_GLM but
  # compute numerical integral of func*b_i(t-fromT) in [fromT , toT] for equivalence with the poisson GLM trick
  # func : function to integrate, func(t, ...)
  # fromT : lower bound (vector)
  # toT   : upper bound (vector)
  # Spline : Spline parameters
  # step : object of class GLMStepParam
  # Nstep : index of the first and last complete band ( fromT[i] < step@cuts[Nstep[i,1]] <= step@cuts[Nstep[i,2]+1] < toT)
  #                                                    ( fromT[i] < step@points[Nstep[i,1]] <= step@points[Nstep[i,2]] < toT)
  # intweightsfunc=NULL, not used, for compatibility with ind_TD_base_NC
  # fail : statuts at T
  # ... : parameters of func()
  if(debug>200) {
    cat("fastinintTD_base_glm\n")
  }
  res<-matrix(0, nrow = length(fromT), ncol = Spline@nbases+Spline@log)
                     # matrix of bases evaluated at the points and T
  Tpoints <- ifelse(fail, toT,  (step@cuts[1+Nstep[,2]]+toT)/2) 
  TBaseattoT <- fevaluate(Spline, Tpoints - fromT, intercept=TRUE)
  TBaseatfromT <- fevaluate(Spline, (step@cuts[Nstep[,1]] - fromT)/2 , intercept=TRUE)
  for(i in 1:length(toT)){
    # vector of evaluated t
    if(Nstep[i,2]>= Nstep[i,1]){
      # at least one complete step
        theT <- c((step@cuts[Nstep[i,1]]+fromT[i])/2,
                  step@points[Nstep[i,1]:Nstep[i,2]] ,
                  Tpoints[i] )
        #evaluated bases 
    allTBase <- fevaluate(Spline, theT - fromT, intercept=TRUE)
                                        # vector of the evaluated functions
        FF <- func(theT, i, ...)
                                        # weights
      w<- c(step@cuts[Nstep[i,1]] - fromT[i],
            step@steps[Nstep[i,1]:Nstep[i,2]],
            toT[i]-step@cuts[1+Nstep[i,2]])
                                        # numerical integration of the complete bands
        res[i,] <- crossprod(w*FF, allTBase)
      }
    else if(Nstep[i,2] - Nstep[i,1] == -1L){
# fromT and toT are in 2 successive bands
      # Nstep[i,2] + 1 =  Nstep[i,1]
        theT <- c((step@cuts[Nstep[i,1]]+fromT[i])/2,
                  Tpoints[i] )
        #evaluated bases 
        allTBase <- fevaluate(Spline, theT - fromT, intercept=TRUE)
                                        # vector of the evaluated functions
        FF <- func(theT, i, ...)
                                        # weights
      w<- c(step@cuts[Nstep[i,1]] - fromT[i],
            toT[i]-step@cuts[Nstep[i,1]])
                                        # numerical integration of the complete bands
        res[i,] <- crossprod(w*FF, allTBase)
      }
    else { #if((Nstep[i,2] - Nstep[i,1]) == -2L){
# fromT and toT are in the same band
      if( fail[i]!=0 ){
        res[i,] <- (toT[i]- fromT[i]) * func(toT[i], i, ...) %*% fevaluate(Spline, toT[i] - fromT[i], intercept=TRUE) #[,,drop=TRUE]
      }
      else {
        res[i,] <- ((toT[i] - fromT[i]) * func((toT[i] + fromT[i])/2, i,  ...)) %*% fevaluate(Spline, (toT[i] - fromT[i])/2 , intercept=TRUE) #[,,drop=TRUE] 
      }
    }
  }
  res  
}

