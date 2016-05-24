multRepl <-
  function(X,label=NULL,dl=NULL,delta = 0.65){
    
    if (is.character(dl)) stop("dl must be a numeric vector or matrix")
    if (is.vector(dl)) dl <- matrix(dl,nrow=1)
    
    if (is.character(X)) stop("X is not a valid data matrix or vector.")
    if (is.null(label)) stop("A value for label must be given")
    if (!is.na(label)){
      if (!any(X==label,na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
      if (label!=0 & any(X==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data set")
      if (any(is.na(X))) stop(paste("NA values not labelled as censored values were found in the data set"))
    }
    if (is.na(label)){
      if (any(X==0,na.rm=T)) stop("Zero values not labelled as censored values were found in the data set")
      if (!any(is.na(X),na.rm=T)) stop(paste("Label",label,"was not found in the data set"))
    }
    if (is.vector(X))
      if (ncol(dl)!=ncol(as.data.frame(matrix(X,ncol=length(X))))) stop("The number of columns in X and dl do not agree")
    if (!is.vector(X)){
      if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
      if ((nrow(dl)>1) & (nrow(dl)!=nrow(X))) stop("The number of rows in X and dl do not agree")
    }
    
    nam <- NULL
    if (!is.null(names(X))) nam <- names(X)
    if (is.vector(X)) X <- as.data.frame(matrix(X,ncol=length(X)))
    
    X[X==label] <- NA
    X <- apply(X,2,as.numeric)
    if (is.vector(X)) X <- as.data.frame(matrix(X,ncol=length(X)))
    
    nn <- nrow(X); p <- ncol(X)
    c <- apply(X,1,sum,na.rm=TRUE)
    
    # Check for closure
    closed <- 0
    if (all( abs(c - mean(c)) < .Machine$double.eps^0.5 )) closed <- 1
    
    if (nrow(dl)==1) dl <- matrix(rep(1,nn),ncol=1)%*%dl
    
    Y <- X
    
    for (i in 1:nn){
      if (any(is.na(X[i,]))){
        z <- which(is.na(X[i,]))
        Y[i,z] <- delta*dl[i,z]
        Y[i,-z] <- (1-(sum(Y[i,z]))/c[i])*X[i,-z]
        X[i,z] <- as.numeric((X[i,-z][1]/Y[i,-z][1]))*Y[i,z]
      }
    }        

    if (!is.null(nam)) names(X) <- nam
       
    if (closed==1){
      X <- t(apply(X,1,function(x) x/sum(x)*c[1]))
    } 
    
    return(as.data.frame(X))
  }
