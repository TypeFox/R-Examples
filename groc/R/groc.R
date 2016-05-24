### groc.R: new groc modelling functions
###
### $Id: 
###
### The top level user function.  Implements a formula interface and calls the
### correct fit function to do the work.
### The function borrows heavily from lm() and mvr().

groc <- function(...)
  UseMethod("groc")

# Code similar to mvr() from pls package
groc.default <- function(formula, ncomp, data, subset, na.action,plsrob = FALSE,
                          method = c("lm","lo","s","lts"),
                          D=NULL,gamma=0.75, 
                          Nc=10,Ng=20,scale=FALSE,Cpp=TRUE,
                          model = TRUE, x = FALSE, y = FALSE, sp = NULL, ...)
{


  Dname <- match.call()$D
  if (is.null(Dname) || nchar(Dname) == 0) {
    Dname <- "covariance"
    D <- .covariance
  } else {
    if (is.character(Dname)) {Dname <- as.name(Dname) ; D <- eval(Dname)}
    if (Dname == "spearman") D <- .spearman
    if (Dname == "kendall") D <- .kendall
    if (Dname == "covariance") D <- .covariance
    if (plsrob) Dname <- "covrob"
  }
  
  ret.x <- x                          # More useful names
  ret.y <- y

  ## Get the model frame
  mf <- match.call(expand.dots = FALSE) # Contains the 'call' used with the values of arguments
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]                # Retain only the named arguments
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame()) # mf is a data.frame that contains y and X
  if (missing(method)) method <- "lm"
  if (plsrob) method <- "lts"
  
  if (method == "model.frame") return(mf)
  ## Get the terms
  mt <- attr(mf, "terms")        # This is to include the `predvars'
                                 # attribute of the terms
  ## Get the data matrices
  Y <- model.response(mf, "numeric")   # Retrieve the y vector values
  if (is.matrix(Y)) {
    if (is.null(colnames(Y)))
      colnames(Y) <- paste("Y", 1:dim(Y)[2], sep = "")
  } else {
    Y <- as.matrix(Y)
    colnames(Y) <- deparse(formula[[2]])
  }
  X <- delete.intercept(model.matrix(mt, mf))
  
  nobj <- dim(X)[1]  
  npred <- dim(X)[2]
  
  ## model.matrix prepends the term name to the colnames of matrices.
  ## If there is only one predictor term, and the corresponding matrix
  ## has colnames, remove the prepended term name:
  if (length(attr(mt, "term.labels")) == 1 &&
      !is.null(colnames(mf[[attr(mt, "term.labels")]])))
    colnames(X) <- sub(attr(mt, "term.labels"), "", colnames(X))
  
  ## Set or check the number of components:
  if (missing(ncomp)) {
    ncomp <- min(nobj - 1, npred)
    ncompWarn <- FALSE              # Don't warn about changed `ncomp'
  } else {
    if (ncomp < 1 || ncomp > min(nobj, npred))
      stop(paste("Invalid number of components, ncomp should be less than ",min(nobj, npred)))
    ncompWarn <- TRUE
  }
  
  ## Select fit function:
  fitFunc <- groc.fit

  ## Fit the model:
  z <- fitFunc(X, Y, ncomp, D, gamma, method, plsrob, Nc, Ng, scale, Cpp, FALSE, 100, sp, ...)
  
  ## Build and return the object:
  class(z) <- "groc"
  z$na.action <- attr(mf, "na.action")
  z$ncomp <- ncomp
  z$method <- method
  z$scale <- scale
  z$call <- match.call()
  z$terms <- mt
  z$plsrob <- plsrob
  z$D <- Dname
  if (model) z$model <- mf
  if (ret.x) z$x <- X
  if (ret.y) z$y <- Y
  
  return(z)
  
}


.simplegrid <- function(r,U,Y,Ng,Nc,D) { 
  
  D.vect <- rep(NA,Ng)
  grid.theta <- rep(NA,Ng)
  
  p <- length(r)
  e <- diag(rep(1,p))
  
  for (i in 1:Nc) {
    for (k in 1:Ng) grid.theta[k] <- - (pi/(2^(i-2))) * (0.5 - (k-1)/Ng) 
    for (j in 1:p) {
      for (l in 1:Ng) {
        theta <- grid.theta[l]
        vect <- cos(theta) * r + sin(theta) * e[,j]
        vect <- vect/sqrt(sum(vect^2))
        D.vect[l] <- D(U%*%vect,Y)
      }
      theta0 <- (grid.theta)[which.max(D.vect)]
      vect <- cos(theta0) * r + sin(theta0) * e[,j]
      r <- vect/sqrt(sum(vect^2))
    }
  }
  
  return(r) 
}

groc.fit <- function(X,Y,ncomp=min(nrow(X)-1,ncol(X)),D=NULL,gamma=0.75,method=NULL,plsrob=FALSE,Nc=10,Ng=20,scale=FALSE,Cpp=TRUE,stripped=FALSE,maxiter=100,sp=NULL,...) {


  tryCatch.W.E <- function(expr)
    {
      W <- NULL
      w.handler <- function(w){ # warning handler
        W <<- w
        invokeRestart("muffleWarning")
      }
      list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
             warning = w.handler),
           warning = W)
    }

  
  Dname <- match.call()$D
  if (is.null(Dname) || nchar(Dname) == 0) {
    Dname <- "covariance"
    D <- .covariance
  } else {
    if (is.character(Dname)) {Dname <- as.name(Dname) ; D <- eval(Dname)}
    if (Dname == "spearman") D <- .spearman
    if (Dname == "kendall") D <- .kendall
    if (Dname == "covariance") D <- .covariance
  }
  
  if (plsrob) {D <- covrob ; Dname <- "covrob" ; method <- "lts"}
  if (is.null(method)) method <- "lm"
  
  X <- as.matrix(X) # n x p
  Y <- as.matrix(Y) # n x q
  if(!stripped) {
  ## Save dimnames:
    dnX <- dimnames(X)
    dnY <- dimnames(Y)
  }

  ## Remove dimnames during calculation.  (Doesn't seem to make a
  ## difference here (2.3.0).)
  dimnames(X) <- dimnames(Y) <- NULL
  
# 1.

  nobj <- n <- nrow(X)
  npred <- p <- ncol(X)
  nresp <- q <- ncol(Y)
  if (ncomp>min(n-1,p)) stop("ncomp should be <= min(n-1,p)")

  ## Center variables:
  Xmeans <- colMeans(X)                                                  
  #U <- sweep(X, 2, colMeans(X)) # n x p , subtract the column mean from X to obtain a centered X
  U <- scale(X,scale=scale) # I think this one (with scaling) works faster!
  Ymeans <- colMeans(Y)
  Y.cent <- scale(Y,scale=FALSE) # n x q  mean centered

# 2.
  if (n<p) { # Note: the result is not the same as prcomp() since prcomp() does not use Dinv (i.e. U <- U %*% V for the last instruction)
    svd.Xc <- svd(U,nu=0,nv=n)
    V <- svd.Xc$v[,-n] # p x n
    #Dinv <- diag(1/svd.Xc$d) # n x n
    U <- U %*% V #%*% Dinv # n x n
  } else { V <- NULL}

  p <- ncol(U)

# 3.
  pred <- matrix(0,nrow=n,ncol=q)
  
# 4.
  Y0 <- Y.cent
  Y0.df <- data.frame(matrix(Y0,ncol=q))
  colnames(Y0.df) <- paste("y",1:q,sep="")
  Tcomp.df <- as.data.frame(matrix(NA,nrow=n,ncol=ncomp))
  colnames(Tcomp.df) <- paste("t",1:ncomp,sep="")
  gam.df <- cbind(Tcomp.df,Y0.df)

# 5.  
  Fcomp <- array(NA,dim=c(n,q,ncomp))
  residuals <- array(NA,dim=c(n,q,ncomp))
  Gobjects <- as.list(1:ncomp)
  for (h in 1:ncomp) Gobjects[[h]] <- as.list(1:p)
  Hobjects <- as.list(1:ncomp)
  for (h in 1:ncomp) Hobjects[[h]] <- as.list(1:p)
  B <- matrix(NA,nrow=ncomp,ncol=p)
  
# 6.  
  Tcomp <- matrix(NA,nrow=n,ncol=ncomp) # n x ncomp. 
  rcomp <- matrix(NA,nrow=p,ncol=ncomp) # p x ncomp.
  
  e <- diag(rep(1,p))
  f <- diag(rep(1,q))

  formule <- ""
  
# 7.  
  for (h in 1:ncomp) { # We search for ncomp components

    if (method %in% c("lm","lts")) formule <- paste(formule," + t",h,sep="")
    if (method %in% "lo") formule <- paste(formule," + lo(t",h,")",sep="")
    if (method %in% "s") formule <- paste(formule," + s(t",h,")",sep="")
    
    ## Grid Algorithm (Begining)
    D.vect <- rep(NA,Ng)
    grid.theta <- rep(NA,Ng)
    if (Cpp) { # We use the fast C/C++ version
      r <- e[,1,drop=FALSE] 
      s <- f[,1,drop=FALSE]
      out <- .Call("grid",as.integer(Nc),as.integer(Ng),as.double(r),as.double(s),
                   U,Y.cent,as.integer(maxiter),D,new.env(),PACKAGE="groc")
      if (out[1]=="failed") stop(paste("The algorithm failed to converge in",maxiter,"iterations"))
      r <- out
    } else {    # We use the R version (historical and pedagogical)
      if(q==1){
        r <- e[,1,drop=FALSE]
        r <- .simplegrid(r,U,Y.cent,Ng,Nc,D)
      }
      if (q>1) {
        r <- e[,1,drop=FALSE] 
        s <- f[,1,drop=FALSE] 
        error <- 1
        nbiter <- 1
        while(error>1e-4){
          if (nbiter > maxiter) stop(paste("The algorithm failed to converge in",maxiter,"iterations"))
          r.old <- r
          s.old <- s
          ytmp <- Y.cent%*%s  
          r <- .simplegrid(r,U,ytmp,Ng,Nc,D)
          xtmp <- U%*%r
          s <- .simplegrid(s,Y.cent,xtmp,Ng,Nc,D)          
          error <- max(abs(r.old-r),abs(s.old-s))
          nbiter <- nbiter + 1
        }
      }
    }
   ## Grid Algorithm (End)

# (b)
    tvec <- U %*% r

# (c)
    Tcomp[,h] <- tvec
    gam.df[,h] <- tvec
    rcomp[,h] <- r

# (d)   
    for (j in 1:q) {
      if (method %in% c("lm")) if (!is.null(tryCatch.W.E(Gobjects[[h]][[j]] <- gam(formula=eval(parse(text=paste("y",j," ~ -1 ",formule,sep=""))),data=gam.df))$warning)) stop(paste("gam() failed! Use ncomp <=",h-1))
      if (method %in% c("lo","s")) Gobjects[[h]][[j]] <- gam(formula=eval(parse(text=paste("y",j," ~ ",formule,sep=""))),data=gam.df,sp=sp[1:h])
        #if (!is.null(tryCatch.W.E(Gobjects[[h]][[j]] <- gam(formula=eval(parse(text=paste("y",j," ~ ",formule,sep=""))),data=gam.df))$warning)) stop(paste("gam() failed! Use ncomp <= ",h-1))
      if (method %in% c("lts")) Gobjects[[h]][[j]] <- lqs(formula=eval(parse(text=paste("y",j," ~ ",formule,sep=""))),data=gam.df,method="lts",quantile=floor(gamma*n) + floor((h+1)/2))
      
      gt <- fitted(Gobjects[[h]][[j]])
      
# (e)

# (f)
      pred[,j] <- gt
      Fcomp[,j,h] <- pred[,j] + Ymeans[j]
      
# (g)
      Y.cent[,j] <- Y0[,j] - gt
      residuals[,j,h] <- Y.cent[,j]
    }
# (h)
    if (plsrob) {
      for (i in 1:p) {
        d <- within(data.frame(x = 1:length(tvec)),{x <- NULL;u=U[,i];tvec=tvec;rm(x)})
        Hobjects[[h]][[i]] <- lqs(u ~ tvec,data=d,method="lts",quantile=floor(gamma*n) + 1)
        U[,i] <- residuals(Hobjects[[h]][[i]])
      }
    } else {
      B[h,] <- (t(tvec) %*% U)/sum(tvec^2)
      U <- U - tvec %*% B[h,]
    }
    
  } # End of (h in 1:ncomp) loop


  ## Add dimnames:
  objnames <- dnX[[1]]
  if (is.null(objnames)) objnames <- dnY[[1]]
  prednames <- dnX[[2]]
  respnames <- dnY[[2]]
  compnames <- paste("Comp", 1:ncomp)
  nCompnames <- paste(1:ncomp, "comps")
  
  dimnames(Tcomp) <- list(objnames, compnames)
  if (nobj>=npred) dimnames(rcomp) <- list(prednames, compnames)
  dimnames(Fcomp) <-  dimnames(residuals) <- list(objnames, respnames, nCompnames)
  
  if (stripped) {
        ## Return as quickly as possible
    return(list(Y=Y,R=rcomp,Gobjects=Gobjects,Hobjects=Hobjects,B=B,Xmeans = Xmeans, Ymeans = Ymeans,D=Dname,V=V,dnnames=dimnames(Fcomp)))
  } else {
    
    return(list(Y=Y,fitted.values=Fcomp,residuals=residuals,T=Tcomp,R=rcomp,Gobjects=Gobjects,Hobjects=Hobjects,B=B,Xmeans=Xmeans,Ymeans=Ymeans,D=Dname,V=V,dnnames=dimnames(Fcomp)))
    
  }
  
}

# Code similar to predict.mvr()
predict.groc <- function(object, newdata, ncomp = object$ncomp,na.action = na.pass, ...) {
  
  plsrob <- object$plsrob
  method <- object$method
  scale <- object$scale 
  
  if (missing(newdata) || is.null(newdata))
    newX <- model.matrix(object)
  else if (is.matrix(newdata)) {
    ## For matrices, simply check dimension:
    if (ncol(newdata) != length(object$Xmeans))
      stop("'newdata' does not have the correct number of columns")
    newX <- newdata
  } else if (is.data.frame(newdata)) {
    newX <- newdata[attr(object$terms,"term.labels")]
  } else {
    Terms <- delete.response(terms(object))
    m <- model.frame(Terms, newdata, na.action = na.action)
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, m)
    newX <- delete.intercept(model.matrix(Terms, m))
  }
  
  nobs <- dim(newX)[1]
  
 # D <- object$D

  oldX <- model.matrix(object)

# 1.
  X <- if (is.null(dim(newX))) t(as.matrix(newX)) else as.matrix(newX) # n x p
  Y <- as.matrix(object$Y)

  n <- nrow(X)
  
  Ymeans <- colMeans(Y)
  U <- scale(X,center=object$Xmeans,scale=scale)

  dim.oldX <- dim(oldX)
  if (dim.oldX[1] < dim.oldX[2]) U <- U %*% object$V
  
# 2.

  p <- ncol(U)
  q <- ncol(Y)
  
# 3.
  pred <- array(0,dim=c(n,q,ncomp))
# 4.  

  tvec <- U %*% object$R[, 1, drop = FALSE] # martin
  colnames(tvec) <- "t1"
  newdata <- within(data.frame(x = 1:n), {x <- NULL; t1=tvec  # Because data.frame() removes the class "nmatrix.1" of a matrix with 1 column.
                                    rm(x)
                                  })
# 7.  
  for (h in 1:ncomp) { # We search for ncomp components

    for (j in 1:q) {
      gt <- unlist(predict(object$Gobjects[[h]][[j]],newdata=newdata))
      pred[,j,h] <- gt + Ymeans[j]        
    }
    
    if (h < ncomp) {
      if(plsrob) {
        newdata2 <- within(data.frame(x = 1:n), {x <- NULL; tvec=tvec  # Because data.frame() removes the class "nmatrix.1" of a matrix with 1 column.
                                                rm(x)
                                              })
        
        for (i in 1:p) {U[, i] <- U[, i]-predict(object$Hobjects[[h]][[i]],newdata=newdata2)}
      } else {
        U <- U - tvec %*% object$B[h,]
      }
      
      tvec <- U%*%object$R[, h + 1, drop = FALSE]
      
      eval(parse(text = paste("newdata$t", h + 1, " <- tvec", sep = "")))
      
    }
        
  } # End of (h in 1:ncomp) loop

  dnPred <- object$dnnames
  dnPred[1] <- dimnames(newX)[1]

  dimnames(pred) <- dnPred
  
  return(pred)
  
}

# Inspired from crossval:
grocCrossval <- function (object, segments = 10, segment.type = c("random", "consecutive", 
                                               "interleaved"), length.seg, 
                      trace = 15,       
                      ...)                                                                   
{

  dnPRESS <- object$dnnames
  
  if (!inherits(object, "groc"))                                          
    stop("`object' not an groc object.")                                
  fitCall <- object$call                                                 
  data <- eval(fitCall$data, parent.frame())                             
  if (is.null(data))                                                     
    stop("`object' must be fit with a `data' argument.")               
  if (!is.null(fitCall$subset)) {                                        
    data <- data[eval(fitCall$subset, parent.frame()), ]               
    object$call$subset <- NULL                                         
  }                                                                      
  if (is.na(match("na.action", names(fitCall)))) {                       
    mf <- model.frame(formula(object), data = data)                    
  }                                                                      
  else {                                                                 
    mf <- model.frame(formula(object), data = data, na.action = fitCall$na.action)
  }                                                                                 
  if (!is.null(NAs <- attr(mf, "na.action"))) {                                     
    data <- data[-NAs, ]                                                          
  }                                                                                 
  Y <- as.matrix(model.response(mf))                                                
  nresp <- dim(Y)[2]   # It's q                                                             
  npred <- length(object$Xmeans)                                                    
  nobj <- nrow(data)                                                                
  
  ## Set up segments
  if (is.list(segments)) {                                                          
    if (is.null(attr(segments, "type")))                                          
      attr(segments, "type") <- "user supplied"                                 
  }                                                                                 
  else {                                                                            
    if (missing(length.seg)) {                                                    
      segments <- cvsegments(nobj, k = segments, type = segment.type)           
    }                                                                             
    else {                                                                        
      segments <- cvsegments(nobj, length.seg = length.seg,                     
                             type = segment.type)                                                  
    }                                                                             
  }
  ncomp <- object$ncomp                                                             
  if (ncomp > nobj - max(sapply(segments, length)) - 1)
    stop("`ncomp' too large for cross-validation.", "\nPlease refit with `ncomp' less than ",
         nobj - max(sapply(segments, length)))
  cvPred <- array(dim = c(nobj,nresp,ncomp))
  
  
  for (n.seg in 1:length(segments)) {
    if (n.seg == 1)
      trace.time <- proc.time()[3]
    seg <- segments[[n.seg]]
    fit <- update(object, data = data[-seg, ])
    pred <- predict(fit, newdata = data[seg,])
    cvPred[seg,,] <- pred
    
    if (n.seg == 1) {
      if (is.numeric(trace)) {
        trace.time <- proc.time()[3] - trace.time
        trace <- trace.time * length(segments) > trace
      }
      if (trace)
        cat("Segment: ")
    }
    if (trace)
      cat(n.seg, "")
  }
  
  if (trace)
    cat("\n")
  
  PRESS <- colSums((cvPred - c(Y))^2) # Equivalent (but faster) than : apply((cvPred - c(Y))^2,FUN=sum,MARGIN=c(2,3))
  PREMAD <- apply(abs(cvPred-c(Y)),MARGIN=c(2,3),FUN=median)
  
  
  dimnames(PREMAD) <- dimnames(PRESS) <- dnPRESS[-1]
  
  objnames <- rownames(data)
  if (is.null(objnames)) objnames <- rownames(Y)
  dimnames(cvPred) <- c(list(objnames), dimnames(fitted(object))[-1])
  
  object$validation <- list(method = "CV", pred = cvPred, 
                            PRESS = PRESS, 
                            PREMAD=PREMAD,
                            RMSEP = sqrt(PRESS/nobj),
                            segments = segments,
                            ncomp = ncomp)
  return(object)
}


plot.groc <- function(x,h=x$ncomp,cex=0.8,...) {
  object <- x
  if (!inherits(object, "groc")) 
    stop("`object' not an groc object.")
  if (h>object$ncomp) stop(paste("`h' must be less than or equal to",object$ncomp))
  q <- ncol(object$Y)
  crit.T <- sqrt(qchisq(.975,df=h))
  crit.r <- sqrt(qchisq(.975,df=q))
  T <- as.matrix(object$T[,1:h])
  tau.T2 <- apply(T,2,scaleTau2,mu.too=TRUE)
  dist.T2 <- sqrt(rowSums(scale(T,center=tau.T2[1,],scale=tau.T2[2,])^2))
  r <- as.matrix(residuals(object)[,,h])
  if (q==1){
    tau.r2 <- apply(r,2,scaleTau2,mu.too=TRUE)
    dist.r2 <- sqrt(rowSums(scale(r,center=tau.r2[1,],scale=tau.r2[2,])^2))
  } else {
    dist.r2 <- sqrt(covRob(r,estim = "pairwiseGK")$dist) 
  }
  index2 <- which(dist.T2 > crit.T | dist.r2 > crit.r)
  lindex2 <- length(index2)
  if(!object$plsrob){
    xlim <- c(0,max(crit.T,dist.T2)*1.05)
    ylim <- c(0,max(crit.r,dist.r2)*1.05)
    plot(dist.T2,dist.r2,xlab="Components",ylab="Residuals",xlim=xlim,ylim=ylim,main="Robust distances",pch=20,...)
    abline(v=crit.T,lty=2)
    abline(h=crit.r,lty=2)
    if (lindex2 != 0) text(dist.T2[index2],dist.r2[index2],index2,cex=cex,pos=3)
  } else {
    plsr.out <- plsr(object$call$formula,ncomp=h,scale=object$scale,data=object$model,method= "simpls")
    T1 <- as.matrix(plsr.out$scores[,1:h])
    r1 <- as.matrix(residuals(plsr.out)[,,h])
    s <- apply(T1,2,sd)
    dist.T1 <- sqrt(rowSums(scale(T1,center=FALSE,scale=s)^2))
    dist.r1 <- sqrt(mahalanobis(r1,center=0,var(r1))) 
    par(mfrow=c(1,2))
    xlim <- c(0,max(crit.T,dist.T1)*1.05)
    ylim <- c(0,max(crit.r,dist.r1)*1.05)
    plot(dist.T1,dist.r1,xlab="Components",ylab="Residuals",xlim=xlim,ylim=ylim,main="Classical distances",pch=20,...)
    abline(v=crit.T,lty=2)
    abline(h=crit.r,lty=2)
    index1 <- which(dist.T1>crit.T | dist.r1>crit.r)
    lindex1 <- length(index1)
    if (lindex1 != 0) text(dist.T1[index1], dist.r1[index1],index1,cex=cex,pos=3)
    xlim <- c(0,max(crit.T,dist.T2)*1.05)
    ylim <- c(0,max(crit.r,dist.r2)*1.05)
    plot(dist.T2,dist.r2,xlab="Components",ylab="Residuals",xlim=xlim,ylim=ylim,main="Robust distances",pch=20,...)
    abline(v=crit.T,lty=2)
    abline(h=crit.r,lty=2)
    if (lindex2 != 0) text(dist.T2[index2],dist.r2[index2],index2,cex=cex,pos=3)
    par(mfrow=c(1,1))
  }
}
