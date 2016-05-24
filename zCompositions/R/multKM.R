multKM <-
  function (X,label=NULL,dl=NULL,n.draws=1000,n.knots=NULL)
  {
    
    if (is.character(dl)) stop("dl must be a numeric vector or matrix")
    if (is.vector(dl)) dl <- matrix(dl,nrow=1)
    
    if ((is.vector(X)) | (nrow(X)==1)) stop("X must be a data matrix")
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
    if (ncol(dl)!=ncol(X)) stop("The number of columns in X and dl do not agree")
    if ((nrow(dl)>1) & (nrow(dl)!=nrow(X))) stop("The number of rows in X and dl do not agree")
    
    if ((!is.null(n.knots)) & (length(n.knots)!=1) & (length(n.knots)!=ncol(X))) stop("The dimensions of n.knots and X do not agree")
    if ((!is.null(n.knots)) & (length(n.knots)==1)) {n.knots <- rep(list(n.knots),ncol(X))}
    
    km.imp <- function(x,dl,...){
      
      who <- is.na(x); w <- which(who)
      
      xcen <- ifelse(who,TRUE,FALSE)
      x[who] <- dl[who]
      
      km.ecdf <- cenfit(x,xcen)
      x.km <- rev(km.ecdf@survfit$time) 
      y.km <- rev(km.ecdf@survfit$surv)
      if (is.null(n.knots.part)) {scdf <- smooth.spline(x.km,y.km)}
      if (!is.null(n.knots.part)) {scdf <- smooth.spline(x.km,y.km,nknots=n.knots.part)}
      scdf.fun <- approxfun(scdf$x,scdf$y)
      inv.scdf <- approxfun(scdf$y,scdf$x)
      
      for (i in 1:length(w)){
        if (dl[w[i]] > min(x[!who])){
          temp <- inv.scdf(runif(n.draws,0,scdf.fun(dl[w[i]])))
          x[w[i]] <- exp(mean(log(temp),na.rm=T))
        }
      }
      return(as.numeric(x))
    }
    
    X[X==label] <- NA
    X <- apply(X,2,as.numeric)
    
    nn <- nrow(X); p <- ncol(X)
    c <- apply(X,1,sum,na.rm=TRUE)
    
    # Check for closure
    closed <- 0
    if (all( abs(c - mean(c)) < .Machine$double.eps^0.5 )) closed <- 1
    
    if (nrow(dl)==1){
      dl <- matrix(rep(1,nn),ncol=1)%*%dl
      est <- dl
    }
    else est <- dl
    
    for (part in 1:p)
    {
      if (any(is.na(X[,part]))) 
      {
        n.knots.part <- n.knots[[part]]
        est[,part] <- km.imp(X[,part],dl[,part],n.draws,n.knots.part)
      }
      else {est[,part] <- 0}
    }
    
    Y <- X
    
    for (i in 1:nn){
      if (any(is.na(X[i,]))){
        z <- which(is.na(X[i,]))
        Y[i,z] <- est[i,z]
        Y[i,-z] <- (1-(sum(Y[i,z]))/c[i])*X[i,-z]
        X[i,z] <- as.numeric((X[i,-z][1]/Y[i,-z][1]))*Y[i,z]
      }
    }   
  
    if (closed==1){
      X <- t(apply(X,1,function(x) x/sum(x)*c[1]))
    }
    
    return(as.data.frame(X))
  }  
    