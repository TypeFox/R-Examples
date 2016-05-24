predict.Mort2Dsmooth <-
function(object, newdata = NULL,
                                 type = c("link", "response"),
                                 se.fit = FALSE, ...){
  ## Input:
  ## object: a Mort2Dsmooth object
  ## newdata: an optional list in which
  ##          to look for x and/or y
  ##          variables with which to predict. If omitted,
  ##          the fitted values are used.
  ## type: the type of prediction required.
  ##       The default ("link") is
  ##       on the scale of the linear predictors;
  ##       the alternative "response" is on the
  ##       scale of the response variable.
  ## se.fit: logical switch indicating if
  ##         standard errors are required
  
  ## Output: a list with components 
  ## fit: Predictions
  ## se.fit: Estimated standard errors
  type <- match.arg(type)
  
  ## NO standard errors ## NO standard errors
  if (!se.fit) {
    ## NO newdata ## NO newdata
    if(is.null(newdata)) {
      pred <- switch(type,
                     link = object$linear.predictors -
                     object$offset,
                     response = object$fitted.values)
    }else{
      ## YES newdata ## YES newdata
      ## check newdata as a list
      if(!is.list(newdata))
        stop("newdata must be a data frame")
      ## check newdata names and variables
      NEWdata <- list(x=object$x, y=object$y)
      namesNEW <- names(NEWdata)
      NEWdata[(namesnew <- names(newdata))] <- newdata
      test <- noNms <- namesnew[!namesnew %in% namesNEW]
      if(length(test) > 0){
        stop("unknown names in newdata: ",
             paste(noNms, collapse =      ", "))
      }
      ## original data
      x <- object$x
      y <- object$y
      Z <- object$Z
      offset <- object$offset      
      ## new data
      new.x <- sort(unique(c(x, newdata$x)))
      new.y <- sort(unique(c(y, newdata$y)))
      ## new dimensions
      new.m <- length(new.x)
      new.n <- length(new.y)
      ## new weight matrix, where old and new coincide
      new.W1 <- matrix(0, new.m, new.n)
      new.W2 <- matrix(0, new.m, new.n)
      whi.x <- which(new.x%in%x)
      whi.y <- which(new.y%in%y)
      new.W1[,whi.y] <- 10
      new.W2[whi.x,] <- 1
      whi <- (new.W1 - new.W2)==9
      new.W <- matrix(0, new.m, new.n)
      new.W[whi] <- 1
      ## new death and offset matrices
      new.Z <- matrix(0, new.m, new.n)
      new.Z[whi] <- Z
      new.offset <- matrix(100, new.m, new.n)
      new.offset[whi] <- offset

      ## new B-spline basis
      ## over x
      if(min(new.x)>=min(x) & max(new.x)<=max(x)){
        ## only interpolation
        new.ndx <- object$ndx[1]
        xl <- min(new.x)
        xr <- max(new.x)
        xmax <- xr + 0.01 * (xr - xl)
        xmin <- xl - 0.01 * (xr - xl)
        new.Bx <- MortSmooth_bbase(x=new.x,
                                   xl=xmin, xr=xmax,
                                   ndx=new.ndx,
                                   deg=object$deg[1])
      }else{
        ## extrapolation + eventual interpolation
        deg <- object$deg[1]
        xl <- min(x)
        xr <- max(x)
        xmax <- xr + 0.01 * (xr - xl)
        xmin <- xl - 0.01 * (xr - xl)
        dx <- (xmax - xmin)/object$ndx[1]
        xl1 <- min(new.x)
        xr1 <- max(new.x)
        minK0 <- xmin-deg:(deg+100)*dx
        minK <- minK0[which(minK0<=xl1)[deg]]
        maxK0 <- xmax+deg:(deg+100)*dx
        maxK <- maxK0[which(maxK0>=xr1)[deg+1]]
        knots <- seq(minK, maxK, by=dx)
        PP <- outer(new.x, knots, MortSmooth_tpower, deg)
        nn <- dim(PP)[2]
        DD <-diff(diag(nn),diff=deg+1)/(gamma(deg+1)*dx^deg)
        new.Bx <- (-1)^(deg + 1) * PP %*% t(DD)
      }
      nbx <- ncol(new.Bx)

      ## matplot(x, object$Bx, t="l", col=1,
      ##         lty=1, xlim=range(new.x), lwd=2)
      ## matlines(new.x, new.Bx, t="l", col=2, lty=2, lwd=3)

      ## over y
      if(min(new.y)>=min(y) & max(new.y)<=max(y)){
        ## only interpolation
        new.ndx <- object$ndx[2]
        yl <- min(new.y)
        yr <- max(new.y)
        ymax <- yr + 0.01 * (yr - yl)
        ymin <- yl - 0.01 * (yr - yl)
        new.By <- MortSmooth_bbase(x=new.y,
                                   xl=ymin, xr=ymax,
                                   ndx=new.ndx,
                                   deg=object$deg[2])
      }else{
        deg <- object$deg[2]
        yl <- min(y)
        yr <- max(y)
        ymax <- yr + 0.01 * (yr - yl)
        ymin <- yl - 0.01 * (yr - yl)
        dx <- (ymax - ymin)/object$ndx[2]
        yl1 <- min(new.y)
        yr1 <- max(new.y)
        minK0 <- ymin-deg:(deg+100)*dx
        minK <- minK0[which(minK0<=yl1)[deg]]
        maxK0 <- ymax+deg:(deg+100)*dx
        maxK <- maxK0[which(maxK0>=yr1)[deg+1]]
        knots <- seq(minK, maxK, by=dx)
        PP <- outer(new.y, knots, MortSmooth_tpower, deg)
        nn <- dim(PP)[2]
        DD <-diff(diag(nn),diff=deg+1)/(gamma(deg+1)*dx^deg)
        new.By <- (-1)^(deg + 1) * PP %*% t(DD)
      }
      nby <- ncol(new.By)

      ## matplot(y, object$By, t="l", col=1,
      ##         lty=1, xlim=range(new.y), lwd=2)
      ## matlines(new.y, new.By, t="l", col=2, lty=2, lwd=3)
      
      ## Row tensors of B-splines basis
      Bx1 <- kronecker(matrix(1, ncol=nbx, nrow=1), new.Bx)
      Bx2 <- kronecker(new.Bx, matrix(1, ncol=nbx,nrow=1))
      RTBx <- Bx1*Bx2
      By1 <- kronecker(matrix(1, ncol=nby, nrow=1), new.By)
      By2 <- kronecker(new.By, matrix(1, ncol=nby, nrow=1))
      RTBy <- By1*By2
      ## penalty stuff
      Dx <- diff(diag(nbx), diff=object$pord[1])
      Dy <- diff(diag(nby), diff=object$pord[2])
      Px <- kronecker(diag(nby), t(Dx)%*%Dx)
      Py <- kronecker(t(Dy)%*%Dy, diag(nbx))

      ## initialize
      new.Z[is.na(new.Z)] <- 0
      ## 0) simple poisson-GLM with ages and years
      ##    only for the interpolation cases
      new.xx <- rep(new.x, new.n)
      new.yy <- rep(new.y, each=new.m)
      fit0 <- glm(round(c(new.Z)) ~ new.xx * new.yy +
                  offset(c(new.offset)),
                  family=poisson, weights=c(new.W))
      etaGLM <- matrix(log(fit0$fitted) - c(new.offset),
                       new.m, new.n)
      eta0 <- log((new.Z + 1)) - new.offset
      eta0[new.W==0] <- etaGLM[new.W==0]
      ## simple ridge penalized regression
      BBx <- solve(t(new.Bx)%*%new.Bx + diag(nbx) * 1e-6,
                   t(new.Bx))
      BBy <- solve(t(new.By)%*%new.By + diag(nby) * 1e-6,
                   t(new.By))
      a.init <- MortSmooth_BcoefB(BBx, BBy, eta0)

      ## estimating
      new.fit <- Mort2Dsmooth_estimate(x=new.x,
                                       y=new.y,
                                       Z=new.Z,
                                       offset=new.offset,
                                       wei=new.W,
                                       psi2=1,
                                       Bx=new.Bx,
                                       By=new.By,
                                       nbx=nbx,nby=nby,
                                       RTBx=RTBx, RTBy=RTBy,
                                     lambdas=object$lambdas,
                                       Px=Px,Py=Py,
                                       a.init=a.init,
                                       MON=FALSE,
                                       TOL1=1e-6,
                                       MAX.IT=50)

      new.eta0 <- matrix(MortSmooth_BcoefB(new.Bx,
                                           new.By,
                                           new.fit$a),
                         new.m,new.n)
      colnames(new.eta0) <- new.y
      rownames(new.eta0) <- new.x
      
      ## matplot(x, log(Z/E), t="p", cex=0.5, pch=1,
      ##         xlim=range(new.x), ylim=range(new.eta0))
      ## matlines(new.x, new.eta0, t="l", lty=2)

      ## matplot(y, t(log(Z/E)), t="p", cex=0.5, pch=1,
      ##         xlim=range(new.y), ylim=range(new.eta0))
      ## matlines(new.y, t(new.eta0), t="l", lty=2)

      ## for(i in 1:new.m){
      ##   bla <- new.eta0[i,]
      ##   if(i%in%whi.x){
      ##     bla1 <- log(Z[whi.x==i,]/E[whi.x==i,])
      ##   }else{
      ##     bla1 <- rep(NA,ncol(Z))
      ##   }
      ##   ran <- range(bla,bla1, na.rm=T)
      ##   plot(new.y, t(bla), t="l", lty=2,
      ##        xlim=range(new.y), ylim=ran,
      ##        main=paste(new.x[i]))
      ##   points(y, t(bla1))
      ##   locator(1)
      ## }
      
      ## for(i in 1:new.n){
      ##   bla <- new.eta0[,i]
      ##   if(i%in%whi.y){
      ##     bla1 <- log(Z[,whi.y==i]/E[,whi.y==i])
      ##   }else{
      ##     bla1 <- rep(NA,nrow(Z))
      ##   }
      ##   ran <- range(new.eta0)#bla,bla1, na.rm=T)
      ##   plot(new.x, bla, t="l", lty=2,
      ##        xlim=range(new.x), ylim=ran,
      ##        main=paste(new.y[i]))
      ##   points(x, t(bla1))
      ##   locator(1)
      ## }
      
      ## in case either x or y is not provided in newdata
      if(is.null(newdata$x)){
        newdata$x <- x}
      if(is.null(newdata$y)){
        newdata$y <- y}
      ## 
      new.eta <- new.eta0[new.x%in%newdata$x,
                          new.y%in%newdata$y]
      new.fitted0 <- exp(new.eta0 + new.offset)
      new.fitted0[!whi] <- NA
      new.fitted <- new.fitted0[new.x%in%newdata$x,
                                new.y%in%newdata$y]
      pred <- switch(type,
                     link = new.eta,
                     response = new.fitted)
    }
  ## YES standard errors ## YES standard errors
  }else{
    ## NO newdata ## NO newdata
    if(is.null(newdata)){
      Bx <- object$Bx
      By <- object$By
      nbx <- ncol(Bx)
      nby <- ncol(By)
      ## Row tensors of B-splines basis
      Bx1 <- kronecker(matrix(1, ncol=nbx, nrow=1), Bx)
      Bx2 <- kronecker(Bx, matrix(1, ncol=nbx,nrow=1))
      RTBx <- Bx1*Bx2
      By1 <- kronecker(matrix(1, ncol=nby, nrow=1), By)
      By2 <- kronecker(By, matrix(1, ncol=nby, nrow=1))
      RTBy <- By1*By2
      ## penalty stuff
      Dx <- diff(diag(nbx), diff=object$pord[1])
      Dy <- diff(diag(nby), diff=object$pord[2])
      Px <- kronecker(diag(nby), t(Dx)%*%Dx)
      Py <- kronecker(t(Dy)%*%Dy, diag(nbx))
      P <- (object$lambdas[1]*Px) + (object$lambdas[2]*Py)
      ## 
      eta <- MortSmooth_BcoefB(Bx, By, object$coef)
      mu <- exp(object$offset + eta)
      W <- mu
      WW <- object$W*W
      BWB <- MortSmooth_BWB(RTBx, RTBy, nbx, nby, WW)
      ## 
      BWB.P1 <- solve(BWB + P)
      se <- matrix(Mort2Dsmooth_se(RTBx=RTBx, RTBy=RTBy,
                                   nbx=nbx, nby=nby,
                                   BWB.P1=BWB.P1),
                   object$m, object$n,
                   dimnames=list(object$x, object$y))
      ## check usage of overdispersion
      check.over <- eval(object$call$overdispersion)
      if(is.null(check.over)){
        check.over <- FALSE
      }
      if(check.over){
        se.inflate <- sqrt(object$psi2)
        se <- se * se.inflate
      }
      fit <- switch(type,
                    link = object$linear.predictors - object$offset,
                    response = object$fitted.values)
      se.fit <- switch(type,
                       link = se,
                       response = object$fitted.values * (exp(se) -1))
    }else{
      ## YES newdata ## YES newdata
      ## check newdata as a list
      if(!is.list(newdata))
        stop("newdata must be a data frame")
      ## check newdata names and variables
      NEWdata <- list(x=object$x, y=object$y)
      namesNEW <- names(NEWdata)
      NEWdata[(namesnew <- names(newdata))] <- newdata
      test <- noNms <- namesnew[!namesnew %in% namesNEW]
      if(length(test) > 0){
        stop("unknown names in newdata: ",
             paste(noNms, collapse =      ", "))
      }
      ## original data
      x <- object$x
      y <- object$y
      Z <- object$Z
      offset <- object$offset
      ## new data
      new.x <- sort(unique(c(x, newdata$x)))
      new.y <- sort(unique(c(y, newdata$y)))
      ## new dimensions
      new.m <- length(new.x)
      new.n <- length(new.y)
      ## new weight matrix, where old and new coincide
      new.W1 <- matrix(0, new.m, new.n)
      new.W2 <- matrix(0, new.m, new.n)
      whi.x <- which(new.x%in%x)
      whi.y <- which(new.y%in%y)
      new.W1[,whi.y] <- 10
      new.W2[whi.x,] <- 1
      whi <- (new.W1 - new.W2)==9
      new.W <- matrix(0, new.m, new.n)
      new.W[whi] <- 1
      ## new death and offset matrices
      new.Z <- matrix(0, new.m, new.n)
      new.Z[whi] <- Z
      new.offset <- matrix(100, new.m, new.n)
      new.offset[whi] <- offset
      ## new B-spline basis
      ## over x
      if(min(new.x)>=min(x) & max(new.x)<=max(x)){
        ## only interpolation
        new.ndx <- object$ndx[1]
        xl <- min(new.x)
        xr <- max(new.x)
        xmax <- xr + 0.01 * (xr - xl)
        xmin <- xl - 0.01 * (xr - xl)
        new.Bx <- MortSmooth_bbase(x=new.x,
                                   xl=xmin, xr=xmax,
                                   ndx=new.ndx,
                                   deg=object$deg[1])
      }else{
        ## extrapolation + eventual interpolation
        deg <- object$deg[1]
        xl <- min(x)
        xr <- max(x)
        xmax <- xr + 0.01 * (xr - xl)
        xmin <- xl - 0.01 * (xr - xl)
        dx <- (xmax - xmin)/object$ndx[1]
        xl1 <- min(new.x)
        xr1 <- max(new.x)
        minK0 <- xmin-deg:(deg+100)*dx
        minK <- minK0[which(minK0<=xl1)[deg]]
        maxK0 <- xmax+deg:(deg+100)*dx
        maxK <- maxK0[which(maxK0>=xr1)[deg+1]]
        knots <- seq(minK, maxK, by=dx)
        PP <- outer(new.x, knots, MortSmooth_tpower, deg)
        nn <- dim(PP)[2]
        DD <-diff(diag(nn),diff=deg+1)/(gamma(deg+1)*dx^deg)
        new.Bx <- (-1)^(deg + 1) * PP %*% t(DD)
      }
      nbx <- ncol(new.Bx)

      ## matplot(x, object$Bx, t="l", col=1,
      ##         lty=1, xlim=range(new.x), lwd=2)
      ## matlines(new.x, new.Bx, t="l", col=2, lty=2, lwd=3)

      ## over y
      if(min(new.y)>=min(y) & max(new.y)<=max(y)){
        ## only interpolation
        new.ndx <- object$ndx[2]
        yl <- min(new.y)
        yr <- max(new.y)
        ymax <- yr + 0.01 * (yr - yl)
        ymin <- yl - 0.01 * (yr - yl)
        new.By <- MortSmooth_bbase(x=new.y,
                                   xl=ymin, xr=ymax,
                                   ndx=new.ndx,
                                   deg=object$deg[2])
      }else{
        deg <- object$deg[2]
        yl <- min(y)
        yr <- max(y)
        ymax <- yr + 0.01 * (yr - yl)
        ymin <- yl - 0.01 * (yr - yl)
        dx <- (ymax - ymin)/object$ndx[2]
        yl1 <- min(new.y)
        yr1 <- max(new.y)
        minK0 <- ymin-deg:(deg+100)*dx
        minK <- minK0[which(minK0<=yl1)[deg]]
        maxK0 <- ymax+deg:(deg+100)*dx
        maxK <- maxK0[which(maxK0>=yr1)[deg+1]]
        knots <- seq(minK, maxK, by=dx)
        PP <- outer(new.y, knots, MortSmooth_tpower, deg)
        nn <- dim(PP)[2]
        DD <-diff(diag(nn),diff=deg+1)/(gamma(deg+1)*dx^deg)
        new.By <- (-1)^(deg + 1) * PP %*% t(DD)
      }
      nby <- ncol(new.By)

      ## matplot(y, object$By, t="l", col=1,
      ##         lty=1, xlim=range(new.y), lwd=2)
      ## matlines(new.y, new.By, t="l", col=2, lty=2, lwd=3)
      

      ## Row tensors of B-splines basis
      Bx1 <- kronecker(matrix(1, ncol=nbx, nrow=1), new.Bx)
      Bx2 <- kronecker(new.Bx, matrix(1, ncol=nbx,nrow=1))
      RTBx <- Bx1*Bx2
      By1 <- kronecker(matrix(1, ncol=nby, nrow=1), new.By)
      By2 <- kronecker(new.By, matrix(1, ncol=nby, nrow=1))
      RTBy <- By1*By2
      ## penalty stuff
      Dx <- diff(diag(nbx), diff=object$pord[1])
      Dy <- diff(diag(nby), diff=object$pord[2])
      Px <- kronecker(diag(nby), t(Dx)%*%Dx)
      Py <- kronecker(t(Dy)%*%Dy, diag(nbx))

      ## initialize
      new.Z[is.na(new.Z)] <- 0
      ## 0) simple poisson-GLM with ages and years
      ##    only for the interpolation cases
      new.xx <- rep(new.x, new.n)
      new.yy <- rep(new.y, each=new.m)
      fit0 <- glm(round(c(new.Z)) ~ new.xx * new.yy +
                  offset(c(new.offset)),
                  family=poisson, weights=c(new.W))
      etaGLM <- matrix(log(fit0$fitted) - c(new.offset),
                       new.m, new.n)
      eta0 <- log((new.Z + 1)) - new.offset
      eta0[new.W==0] <- etaGLM[new.W==0]
      ## simple ridge penalized regression
      BBx <- solve(t(new.Bx)%*%new.Bx + diag(nbx) * 1e-6,
                   t(new.Bx))
      BBy <- solve(t(new.By)%*%new.By + diag(nby) * 1e-6,
                   t(new.By))
      a.init <- MortSmooth_BcoefB(BBx, BBy, eta0)
      ## estimating
      new.fit <- Mort2Dsmooth_estimate(x=new.x,
                                       y=new.y,
                                       Z=new.Z,
                                       offset=new.offset,
                                       wei=new.W,
                                       psi2=1,
                                       Bx=new.Bx,
                                       By=new.By,
                                       nbx=nbx,nby=nby,
                                       RTBx=RTBx, RTBy=RTBy,
                                     lambdas=object$lambdas,
                                       Px=Px,Py=Py,
                                       a.init=a.init,                                       
                                       MON=FALSE,
                                       TOL1=1e-6,
                                       MAX.IT=50)

      new.eta0 <- matrix(MortSmooth_BcoefB(new.Bx,
                                           new.By,
                                           new.fit$a),
                         new.m,new.n)
      colnames(new.eta0) <- new.y
      rownames(new.eta0) <- new.x
      ## in case either x or y is not provided in newdata
      if(is.null(newdata$x)){
        newdata$x <- x}
      if(is.null(newdata$y)){
        newdata$y <- y}
      ## 
      new.eta <- new.eta0[new.x%in%newdata$x,
                          new.y%in%newdata$y]
      new.fitted0 <- exp(new.eta0 + new.offset)
      new.fitted0[!whi] <- NA
      new.fitted <- new.fitted0[new.x%in%newdata$x,
                                new.y%in%newdata$y]
      ## standard errors
      BWB.P1 <- solve(new.fit$BWB + new.fit$P)
      new.se0 <- matrix(Mort2Dsmooth_se(RTBx=RTBx,
                                        RTBy=RTBy,
                                        nbx=nbx, nby=nby,
                                        BWB.P1=BWB.P1),
                        new.m, new.n,
                        dimnames=list(new.x, new.y))
      new.se <- new.se0[new.x%in%newdata$x,
                        new.y%in%newdata$y]
      ## check usage of overdispersion
      check.over <- eval(object$call$overdispersion)
      if(is.null(check.over)){
        check.over <- FALSE
      }
      if(check.over){
        se.inflate <- sqrt(object$psi2)
        new.se <- new.se * se.inflate
      }
      fit <- switch(type,
                    link = new.eta,
                    response = new.fitted)   
      
      se.fit <- switch(type,
                       link = new.se,
                       response =new.fitted*(exp(new.se)-1))
    }
    pred <- list(fit=fit, se.fit=se.fit)
  }
  return(pred)
}
