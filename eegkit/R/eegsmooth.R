eegsmooth <- 
  function(voltage,space=NULL,time=NULL,nknots=NULL,rparm=NULL,
           lambdas=NULL,skip.iter=TRUE,se.fit=FALSE,rseed=1234){
    ###### Spatial and/or temporal smoothing of EEG data
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: February 16, 2015
    
    ### check voltage
    voltage <- as.matrix(voltage)
    nv <- nrow(voltage)
    if(ncol(voltage)!=1L){stop("Input 'voltage' must be column vector.")}
    
    ### check and fit model
    if(is.null(space[1]) & is.null(time[1])){
      stop("You must input either 'space' or 'time' to smooth.")
    } else if(is.null(space[1])){
      
      # initial checks
      time <- as.matrix(time)
      if(ncol(time)!=1L){stop("Input 'time' must be column vector.")}
      if(nv!=nrow(time)){stop("Inputs 'voltage' and  'time' must have same number of rows.")}
      
      # check knots and rparm
      if(is.null(nknots)){nknots <- 30L} else {nknots <- as.integer(nknots[1])}
      if(is.null(rparm)){rparm <- 0.001}
      
      # fit model
      eegmod <- bigspline(time,voltage,nknots=nknots,rparm=rparm,
                          lambdas=lambdas,se.fit=se.fit,rseed=rseed)
      
    } else if(is.null(time[1])){
      
      # initial checks
      space <- as.matrix(space,rownames=0)
      if(ncol(space)!=3L){stop("Input 'space' must be 3-column matrix of spatial coordinates.")}
      if(nv!=nrow(space)){stop("Inputs 'voltage' and  'space' must have same number of rows.")}
      
      # check knots and rparm
      if(is.null(nknots)){nknots <- 100L} else {nknots <- as.integer(nknots[1])}
      if(is.null(rparm)){rparm <- 0.1}
      
      # fit model
      eegmod <- bigtps(space,c(voltage),nknots=nknots,rparm=rparm,
                       lambdas=lambdas,se.fit=se.fit,rseed=rseed)
      
    } else {
      
      # initial checks
      time <- as.matrix(time)
      if(ncol(time)!=1L){stop("Input 'time' must be column vector.")}
      if(nv!=nrow(time)){stop("Inputs 'voltage' and  'time' must have same number of rows.")}
      space <- as.matrix(space)
      if(ncol(space)!=3L){stop("Input 'space' must be 3-column matrix of spatial coordinates.")}
      if(nv!=nrow(space)){stop("Inputs 'voltage' and  'space' must have same number of rows.")}
      
      # check knots and rparm
      if(is.null(nknots)){nknots <- 500L} else {nknots <- as.integer(nknots[1])}
      if(is.null(rparm)){rparm <- list(space=0.1,time=0.001)}
      
      # fit model
      type <- list(space="tps",time="cub")
      eegmod <- bigssa(voltage~space*time,nknots=nknots,type=type,rparm=rparm,
                       lambdas=lambdas,skip.iter=skip.iter,se.fit=se.fit,rseed=rseed)
      
    } # end if(is.null(space[1]) & is.null(time[1]))
    
    eegmod
    
  }