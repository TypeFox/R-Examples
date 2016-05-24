eegica <-
  function(X,nc,center=TRUE,maxit=100,tol=1e-6,Rmat=diag(nc),
           type=c("time","space"),method=c("imax","fast","jade"),...){
    ###### Independent Component Analysis for EEG data
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: February 16, 2015
    
    ### initial checks
    X <- as.matrix(X)
    type <- type[1]
    if(any(type==c("time","space"))==FALSE){stop("Input 'type' must be 'time' or 'space'.")}
    if(type=="time"){X <- t(X)}
    method <- method[1]
    if(!any(method==c("imax","fast","jade"))){stop("Input 'method' must be 'imax', 'fast', or 'jade'.")}
    
    ### call ica algorithm
    if(method=="imax"){
      imod <- icaimax(X,nc,center,maxit,tol,Rmat,...)
    } else if(method=="fast"){
      imod <- icafast(X,nc,center,maxit,tol,Rmat,...)
    } else {
      imod <- icajade(X,nc,center,maxit,tol,Rmat)
    }
    return(imod <- c(imod,type=type,method=method))
    
  }