cumresMaxSet <- function(m,var,...) {
  if (!require("lava")) stop("package 'lava' not available")
  A <- t(lava::index(m)$A)
  Afix <- A; Afix[t(lava::index(m)$M0)==1] <- 0
  A[A!=0] <- 1
  P <- lava::index(m)$P
  k <- nrow(A)
  I <- diag(k)
  B <- rbind(I,solve(I-A))
  VV <- B%*%P%*%t(B)
  u.var <- lava::index(m)$vars
  V0 <- VV[seq(length(lava::vars(m))),seq(length(lava::vars(m)))+length(lava::vars(m))][,lava::endogenous(m)]
  rownames(V0) <- lava::vars(m)
  names(which(V0[var,]==0))
}



##' Cumulative residual processes for structural equation models
##' 
##' Calculates GoF statistics based on cumulative residual processes for
##' structural equation models fitted with the \code{lava} package.
##' 
##' With \code{y} and \code{x} given as functions the user can decide which
##' variables to use in the prediction of the outcome and predictor (use the
##' \code{predict} method as below).
##' 
##' @param model \code{lvm} object
##' @param y A formula specifying the association to be checked. Alternatively
##' the outcome specified as a function or a string with the name of the outcome
##' in the model.
##' @param x Predictor. A function, vector or character
##' @param full If FALSE the prediction, Pr, of the variable that are ordered
##' after is only calculated based on the conditional distribution given
##' covariates. If TRUE the conditional expectation is based on the largest set
##' of covariates and endogenous variables such that the residual and Pr are
##' uncorrelated.
##' @param data data.frame (default is the model.frame of the model)
##' @param p Optional parameter vector
##' @param R Number of processes to simulate
##' @param b Moving average parameter
##' @param plots Number of processes to save for use with the plot method
##' @param seed Random seed
##' @param \dots Additional arguments parsed on to lower-level functions
##' @method cumres lvmfit
##' @export
##' @return Returns a \code{cumres} object with associated
##' \code{plot},\code{print},\code{confint} methods
##' @author Klaus K. Holst
##' @references B.N. Sanchez and E. A. Houseman and L. M. Ryan (2009)
##' \emph{Residual-Based Diagnostics for Structural Equation Models}. Biometrics
##' Volume 65 (1), pp 104-115.
##' @keywords models regression
##' @examples
##' 
##' \donttest{
##' library(lava)
##' m <- lvm(list(c(y1,y2,y3)~eta,eta~x)); latent(m) <- ~eta
##' ## simulate some data with non-linear covariate effect
##' functional(m,eta~x) <- function(x) 0.3*x^2
##' d <- sim(m,100)
##' 
##' e <- estimate(m,d)
##' ## Checking the functional form of eta on x
##' g <- cumres(e,eta~x,R=1000)
##' plot(g)
##' 
##' x <- function(p) predict(e,x=~y2+y3,p=p)[,"eta"]
##' ## Checking the functional form of y1 on eta
##' cumres(e,y1~eta,R=1000)
##' g <- cumres(e,"y1",x=x,R=1000)
##' plot(g)
##' 
##' }
cumres.lvmfit <- function(model,y,x,full=FALSE,
                          data=model.frame(model),p,
                          R=1000, b=0, plots=min(R,50), seed=round(runif(1,1,1e9)),
                          ...) {
  if (!require("lava")) stop("package 'lava' not available")
  if (!require("numDeriv")) stop("package 'numDeriv' not available")    

  cl <- match.call()
  
  if (class(y)[1]=="formula") {
    y <- lava::getoutcome(y)
    x <- attributes(y)$x
    y <- lava::decomp.specials(y)
  }

  if (is.list(y) && class(y[[1]])[1]=="formula") {
    yy <- y
    y <- list(); x <- list()
    res <- c()
    cl[[1]] <- as.name("cumres")
    for (i in 1:length(yy)) {
      resp <- lava::getoutcome(yy[[i]])
      y <- lava::decomp.specials(resp)
      x <- attributes(resp)$x
      cl$y <- y; cl$x <- x
      res <- c(res, list(eval.parent(cl)))      
    }
    return(res)
  }
  
  cl$y <- y; cl$x <- x
  cl[[1]] <- as.name("cumres")
  if ((is.character(y) | is.list(y)) & length(y)>1) {
    res <- c()
    iy <- 0;
    for (i in y) {
      iy <- iy+1
      cl$y <- i
      yname <- ifelse(is.character(i),i,paste("y",iy,sep=""))
      if (is.character(x) | is.list(x)) {
        ix <- 0
        for (j in x) {
          ix <- ix+1
          cl$x <- j
          xname <- ifelse(is.character(j),j,paste("x",ix,sep=""))          
          newres <- list(eval.parent(cl))
          names(newres) <- paste(yname,xname,sep="<-")
          res <- c(res,newres)
        }
      } else {
        xname <- ifelse(is.character(j),j,"x1")          
        newres <- list(eval.parent(cl))
        names(newres) <- paste(yname,xname,sep="<-")
        res <- c(res,newres)
      }
    }
    return(res)
  }
  
  if ((is.character(x) | is.list(x)) & length(x)>1) {
    res <- c()
    ix <- 0
    yname <- ifelse(is.character(y),y,"y1")
    for (j in x) {
      ix <- ix+1
      cl$x <- j
      xname <- ifelse(is.character(j),j,paste("x",ix,sep=""))          
      newres <- list(eval.parent(cl))
      names(newres) <- paste(yname,xname,sep="<-")
      res <- c(res,newres)
    }
    return(res)
  }

  if (missing(p))
    p <- lava::pars(model)
  x0 <- x
  if (is.function(x))
    x0 <- x(p)
  if (is.character(x)) {
    if (full) {
      predictby <- cumresMaxSet(model,y,...)
      x0 <- predict(model,predictby)[,x]
    } else {
      if (x %in% lava::latent(model)) {
        x0 <- predict(model,~1)[,x]
      } else {
        x0 <- data[,x]
      }
    }
  }
  n <- length(x0)

  hatW.MC <- function(x) {
    ord <- order(x)
    output <- .C("W2",
                 R=as.integer(R), ## Number of realizations
                 b=as.double(b), ## Moving average parameter
                 n=as.integer(n), ## Number of observations
                 npar=as.integer(length(p)), ## Number of parameters (columns in design)
                 xdata=as.double(x), ## Observations to cummulate after 
                 rdata=as.double(r), ## Residuals (one-dimensional)
                 betaiiddata=as.double(beta.iid), ## Score-process
                 etarawdata=as.double(grad), ## Eta (derivative of terms in cummulated process W)
                 plotnum=as.integer(plots), ## Number of realizations to save (for later plot)
                 seed=as.integer(seed), ## Seed (will probably rely on R's randgen in future version)
                 KS=as.double(0), ## Return: Kolmogorov Smirnov statistics for each realization
                 CvM=as.double(0), ## Return: Cramer von Mises statistics for each realization
                 Wsd=as.double(numeric(n)), ## Return: Pointwise variance of W(x)
                 cvalues=as.double(numeric(R)), ## Return: value for each realization s.t.  +/- cvalue * Wsd contains W*
                 Ws=as.double(numeric(plots*n)), ## Return: Saved realizations (for plotting function)
                 Wobs=as.double(numeric(n)) ## Observed process
                 , PACKAGE="gof")
    return(list(output=output,x=x[ord]))
  }
  
  if (is.function(y)) {
    r <- y(p)
    grad <- attributes(r)$grad
    if (is.null(grad))
      grad <- -numDeriv::jacobian(y,p,method=lava::lava.options()$Dmethod)
  } else {
    myres <- function(p) {
      rr <- predict(model,vars(model),residual=TRUE,p=p)
      if (is.matrix(rr)) return(rr[,y]) else return(rr)          
    }
    r <- myres(p)
    grad <- -numDeriv::jacobian(myres,p,method=lava::lava.options()$Dmethod)
  }     

  ord <- order(x0)  
  x0 <- x0[ord]
  grad <- grad[ord,]
  r <- r[ord]
  ##     Ii <- solve(information(model,p), num=TRUE)
  ## Ii <- model$vcov
  ## Score <- lava::score(model,data=data,p=p,indiv=TRUE,weight=lava::Weight(model))[ord,]
  ## beta.iid <- Ii%*%t(Score)
  ## beta.iid[is.na(beta.iid)] <- 0
  beta.iid <- t(iid(model)[ord,])/nrow(data);
  
  onesim <- hatW.MC(x0)
  What <- matrix(onesim$output$Ws,nrow=n);

  if (!is.character(y))  y <- NULL
  res <- with(onesim$output,
              list(W=cbind(Wobs), What=list(What),
                   x=cbind(x0),
                   KS=KS, CvM=CvM,
                   R=R, n=n, sd=cbind(Wsd), 
                   cvalues=cbind(cvalues), variable=ifelse(is.character(x),x,"x"),
                   response=y,
                   type="sem",
                   model=class(model)[1])
              )
  class(res) <- "cumres"
  
  return(res)    
}

  

################################################################################

