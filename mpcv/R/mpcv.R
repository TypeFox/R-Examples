mpcv <-
  function(x, indepvar=1, LSL, USL, Target, alpha=0.0027, distance = list( "mahalanobis", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"), 
           n.integr=100,  coef.up, coef.lo)
  {
    n.coef <- 5 # No of coeff.

    if(missing(x)){
      stop("argument 'x' is missing")
    }
    if(!is.matrix(x)){
      stop("argument 'x' is not a matrix")
    }
    if(!is.numeric(x)){
      stop("argument 'x' is not numeric")
    }
    if (missing(indepvar)) {
      IndepId <- 1
    }
    else  {
      if(is.numeric(indepvar) && indepvar<=ncol(x) && indepvar>0) 
        IndepId <- indepvar
      else 
        IndepId <- which(colnames(x)==indepvar)
    }
    if(length(IndepId)==0){
      stop("argument 'indepvar': undefined variable selected")
    }
    if (missing(Target)) {
      Target <- LSL + (USL - LSL)/2
    }
    if (missing(coef.up)) {
      coef.up <- rep(0, ncol(x))
      coef.up[IndepId] <- NA
    }
    if (missing(coef.lo)) {
      coef.lo <- rep(0, ncol(x))
      coef.lo[IndepId] <- NA
    }
    if(any(LSL>=Target) || any(USL<=Target))
      stop("Target is outside the limits")
    if (missing(distance))
      distance <- "mahalanobis"
    else
      distance <- match.arg(distance)    
    
    n.examples <- nrow(x)
    n.features <- ncol(x)
    n.del <- round(n.examples*alpha,0)
    
    if(n.del > 0)
    {
      colM <- apply(x, 2, median) # originally: colM <- colMeans(x)
      if(distance == "mahalanobis") {
        S <- cov(x)
        d <- sqrt(apply(x, 1, function(x) t(x-colM) %*% solve(S) %*% (x-colM)))
      }
      else {
        colM_x <- rbind(colM, x)
        d <- as.matrix(dist(colM_x, distance))[-1,1]
      }
      idx <- which(d %in% sort(unique(d), decreasing=TRUE)[1:n.del])
      x <- x[-idx,]
    }
    
    n.examples <- nrow(x)
    x <- as.matrix(x[order(x[,IndepId]),])
    
    # -------- one sided models  -------------------
    # upper model
    d.u <- rbind(n.examples,-n.examples, sum(x[,IndepId]), -sum(x[,IndepId]), -sum(x[,IndepId]^2))
    A.u <- matrix(1,nrow=n.examples)
    A.u <- cbind(A.u, matrix(-1,nrow=n.examples), matrix(x[,IndepId]), -matrix(x[,IndepId]), -matrix(x[,IndepId]^2))
    cond.u <- matrix(">=",ncol=n.examples)
    # lower model
    d.l <- rbind(n.examples,-n.examples, sum(x[,IndepId]), -sum(x[,IndepId]), sum(x[,IndepId]^2))
    A.l <- matrix(1,nrow=n.examples)
    A.l <- cbind(A.l,-1, matrix(x[,IndepId]), -matrix(x[,IndepId]), matrix(x[,IndepId]^2))
    cond.l <- matrix("<=",ncol=n.examples)
    
    
    add.coef <- c(rep(0, n.coef-1), 1) # last coef have to be > 0
    A.u <- rbind(A.u, add.coef) 
    cond.u <- cbind(cond.u, ">=") 
    A.l <- rbind(A.l, add.coef)
    cond.l <- cbind(cond.l, ">=") 
    
    a.u <- matrix(ncol=n.coef, nrow=n.features, dimnames=list(colnames(x),NULL))
    a.l <- matrix(ncol=n.coef, nrow=n.features, dimnames=list(colnames(x),NULL))
    
    ForIds <- c(1:n.features)
    ForIds <- ForIds[-IndepId]
    for(i in ForIds)
    {
      # coeff of models (up&lo)
      a.u[i,] <- lp("min", d.u, A.u, cond.u, c(x[,i], coef.up[i]))$solution
      a.l[i,] <- lp("max", d.l, A.l, cond.l, c(x[,i], coef.lo[i]))$solution
    }
    
    marginal.median <- apply(x, 2, median)
    
    # integration models
    points.integr <- seq(min(x[,IndepId]), max(x[,IndepId]), length=n.integr)
    models.integr.u <- matrix(ncol=n.features, nrow=length(points.integr))
    models.integr.l <- matrix(ncol=n.features, nrow=length(points.integr))
    for(i in ForIds)
    {
      models.integr.u[,i] <- a.u[i,1]-a.u[i,2]+a.u[i,3]*points.integr-a.u[i,4]*points.integr-a.u[i,5]*points.integr^2
      models.integr.l[,i] <- a.l[i,1]-a.l[i,2]+a.l[i,3]*points.integr-a.l[i,4]*points.integr+a.l[i,5]*points.integr^2
    }
    
    # width of integr. models intervals
    # dependent vars
    R.m <- abs(models.integr.u -models.integr.l)
    # indep. var
    R.m[,IndepId] <- c(abs(diff(points.integr)),0)
    
    # numerical integr.
    integr <- 1
    for(i in 1:n.features) # volumes of intervals
      integr <- integr*R.m[,i]
    # summary volume
    v.m <- sum(integr)
    
    # ULS & LSL volume
    v.T <- prod(USL - LSL)
    
    models.max <- apply(x, 2, max)
    models.min <- apply(x, 2, min)
            
    CpV <- (v.m/v.T)^(1/n.features)* 100
    CpV <- round(CpV)
        
    
    # is the median greater of Target
    median.position <- Target < marginal.median
    
    # lengths of the tolerances on the side on which the median is located
    toler.length <- c()
    max.shift <- max.shift.id <- 0
    for(i in 1:n.features)
    {
      if(median.position[i]) 
        toler.length[i] <- USL[i] - Target[i]
      else
        toler.length[i] <- Target[i] - LSL[i]
      shift <- abs(Target[i]- marginal.median[i])/toler.length[i]
      if(shift > max.shift)
      {
        max.shift <- shift
        max.shift.id <- i
      }
    }
    
    PS <- max.shift * 100
    PS <- as.integer(round(PS))
    PSVar <- colnames(x)[max.shift.id]

    PD.matrix <- rbind((models.max-Target)/(USL-Target), (Target-models.min)/(Target-LSL))
    PD.max <- max(PD.matrix)
    PD.which <- which(PD.matrix == PD.max, arr.ind = TRUE)[2]
    PD <- round(PD.max * 100)
    PDVar <- colnames(x)[PD.which]

    
    res <- list(CpV=CpV, PS=PS, PSvar=PSVar, PD=PD, PDvar=PDVar, coef.lo=a.l, coef.up=a.u, x=x, Target=Target, LSL=LSL, USL=USL, indepvar=IndepId)
    class(res) <-"mpcv"
    return(res)
  }
