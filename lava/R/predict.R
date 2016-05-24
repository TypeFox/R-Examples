##' @export
predict.lvmfit <- function(object,x=NULL,y=NULL,data=model.frame(object),p=coef(object),...) {
  predict(Model(object),x=x,y=y,p=p,data=data,...)
}


##' @export
predict.lvm.missing <- function(object,x=NULL,y=NULL,data=model.frame(object),p=coef(object),...) {
    idx <- match(coef(Model(object)),names(coef(object)))
    xx <- exogenous(object)
    p <- p[idx]
    if (!is.null(x)) {
        if (inherits(x,"formula")) {
            xy <- getoutcome(x)
            if (length(xy)>0) {
                if (is.null(y)) y <- decomp.specials(xy)
            }
            x <- attributes(xy)$x
      }
      x <- intersect(x,endogenous(object))
      if (is.null(y))
          y <- setdiff(vars(object),c(x,xx))
    } 
    obs0 <- !is.na(data[,x,drop=FALSE])
    data[,xx][which(is.na(data[,xx]),arr.ind=TRUE)] <- 0
    pp <- predict.lvmfit(object,x=x,y=y,data=data,p=p,...)
    if (all(obs0)) return(pp)
    
    if (!requireNamespace("mets",quietly=TRUE)) stop("Requires 'mets'")
    obs <- mets::fast.pattern(obs0)
    res <- matrix(nrow=nrow(data),ncol=NCOL(pp))
    for (i in seq(nrow(obs$pattern))) {
        jj <- which(obs$pattern[i,]==1)
        ii <- which(obs$group==i-1)
        res[ii,] <- predict.lvmfit(object,...,p=p,x=x[jj],y=y,data=data[ii,,drop=FALSE])[,colnames(pp),drop=FALSE]
    }
    attributes(res) <- attributes(pp)
    return(res)
}


##' Prediction in structural equation models
##'
##' Prediction in structural equation models
##' @param object Model object
##' @param x optional list of (endogenous) variables to condition on
##' @param y optional subset of variables to predict
##' @param residual If true the residuals are predicted
##' @param p Parameter vector
##' @param data Data to use in prediction
##' @param path Path prediction
##' @param quick If TRUE the conditional mean and variance given covariates are returned (and all other calculations skipped)
##' @param \dots Additional arguments to lower level function
##' @examples
##' m <- lvm(list(c(y1,y2,y3)~u,u~x)); latent(m) <- ~u
##' d <- sim(m,100)
##' e <- estimate(m,d)
##'
##' ## Conditional mean (and variance as attribute) given covariates
##' r <- predict(e)
##' ## Best linear unbiased predictor (BLUP)
##' r <- predict(e,vars(e))
##' ##  Conditional mean of y3 giving covariates and y1,y2
##' r <- predict(e,y3~y1+y2)
##' ##  Conditional mean  gives covariates and y1
##' r <- predict(e,~y1+y2)
##' ##  Predicted residuals (conditional on all observed variables)
##' r <- predict(e,vars(e),residual=TRUE)
##'
##' @method predict lvm
##' @aliases predict.lvmfit
##' @export
predict.lvm <- function(object,x=NULL,y=NULL,residual=FALSE,p,data,path=FALSE,quick=is.null(x)&!(residual|path),...) {
  ## data = data.frame of exogenous variables

  if (!quick && !all(exogenous(object)%in%colnames(data))) stop("data.frame should contain exogenous variables")
  m <- moments(object,p,data=data,...)
  if (quick) { ## Only conditional moments given covariates
      ii <- index(object)
      P.x <- m$P; P.x[ii$exo.idx, ii$exo.idx] <- 0
      Cy.x <- (m$IAi%*% tcrossprod(P.x,m$IAi))[ii$endo.idx,ii$endo.idx,drop=FALSE]
      X <- ii$exogenous
      mu.0 <- m$v; mu.0[ii$exo.idx] <- 0
      if (length(X)>0) {
          mu.x <- matrix(0,ncol=nrow(data),nrow=length(mu.0))
          mu.x[ii$exo.idx,] <- t(data[,X,drop=FALSE])
          xi.x <- (m$IAi[ii$endo.obsidx,,drop=FALSE]%*%(mu.0 + mu.x))
      } else {
          xi.x <- m$xi%*%rep(1,nrow(data))
          rownames(xi.x) <- ii$endogenous
          ##xi.x <- matrix(as.vector(m$IAi[ii$endo.obsidx,]%*%mu.0),ncol=nrow(data),nrow=length(mu.0))
          ##rownames(xi.x) <- names(mu.0)
      }
      return(structure(t(xi.x),cond.var=Cy.x,
                       p=m$p,
                       e=m$e))
  }

    
  X <- exogenous(object)
  Y <- setdiff(manifest(object), X)
  if (path) {
    X <- colnames(data)
    Y <- setdiff(Y,X)
    idx <- which(vars(object)%in%X)
    if (length(Y)==0) stop("New data set should only contain exogenous variables and a true subset of the endogenous variables for 'path' prediction.")
    A <- t(m$A)
    A[,idx] <- 0 ## i.e., A <- A%*%J
    IAi <- solve(diag(nrow=nrow(A))-t(A))
    mu.0 <- m$v;
    mu.0[X] <- 0
    mu.x <- matrix(0,ncol=nrow(data),nrow=length(mu.0))
    mu.x[idx,] <- t(data[,vars(object)[idx],drop=FALSE])
    pred <- t(IAi%*%(mu.0 + mu.x))
    return(pred)
    ##  Y <- endogenous(object,top=TRUE)
    ##  X <- setdiff(manifest(object),Y)
  }


  IAi <- m$IAi
  P <- m$P
  X.idx <- match(X,manifest(object))
  eta.idx <- match(latent(object),vars(object))
  obs.idx <- match(manifest(object),vars(object))
  X.idx.all <- match(X, vars(object))
  Y.idx.all <- match(Y, vars(object))

  ## Calculation of conditional variance given X=x
  P.x <- m$P; P.x[X.idx.all, X.idx.all] <- 0
  C.x <- (IAi%*% P.x %*%t(IAi))
  Cy.x <- C.x[Y.idx.all,Y.idx.all,drop=FALSE]
  ## Calculation of conditional mean given X=x
  G <- m$J%*%IAi
  mu.0 <- m$v; mu.0[X.idx.all] <- 0
  if (length(X)>0) {
    xs <- data[,X,drop=FALSE]
    mu.x <- apply(xs, 1, FUN=function(i) {res <- rep(0,length(mu.0)); res[X.idx.all] <- i; res})
    xi.x <- (IAi%*%(mu.0 + mu.x))
  } else {
    xi.x <- matrix(as.vector(IAi%*%mu.0),ncol=nrow(data),nrow=length(mu.0))
    rownames(xi.x) <- names(mu.0)
  }

  attr(xi.x,"cond.var") <- Cy.x
  if (path) return(t(xi.x))
  Ey.x <- xi.x[Y.idx.all,,drop=FALSE]
  Eeta.x <- xi.x[eta.idx,,drop=FALSE]
  Cy.epsilon <- P.x%*%t(IAi) ## Covariance y,residual
  Czeta.y <- Cy.epsilon[eta.idx,index(object)$endo.idx]
  A <- m$A
  IA <- diag(nrow=nrow(A))-t(A)

  y0 <- intersect(Y,colnames(data))
  ys <- data[,y0,drop=FALSE]
  y0.idx <- match(y0,Y)
  ry <- t(ys)-Ey.x[y0.idx,,drop=FALSE]

  if (!is.null(x)) {
      if (inherits(x,"formula")) {
          xy <- getoutcome(x)
          if (length(xy)>0) {
              if (is.null(y)) y <- decomp.specials(xy)
          }
          x <- attributes(xy)$x
      }
      if (length(x)==0) {
          if (!is.null(y)) {
              xi.x <- xi.x[y,,drop=FALSE]
              attr(xi.x,"cond.var") <- Cy.x[y,y,drop=FALSE]
          }
          return(t(xi.x))
      }
      x <- intersect(x,endogenous(object))
      if (is.null(y))
          y <- setdiff(vars(object),c(x,exogenous(object)))

      E.x <- xi.x[y,,drop=FALSE] + C.x[y,x]%*%solve(C.x[x,x])%*%ry[x,,drop=FALSE]
      if (residual) {
          Vhat <- matrix(0, nrow(data), length(vars(object))); colnames(Vhat) <- vars(object)
          Vhat[,obs.idx] <- as.matrix(data[,manifest(object),drop=FALSE])
          Vhat[,y] <- t(E.x)
          return(t((IA%*%t(Vhat)-m$v)))
      }
      res <- t(E.x); colnames(res) <- y
      attr(res,"cond.var") <-
          C.x[y,y,drop=FALSE]-C.x[y,x,drop=FALSE]%*%solve(C.x[x,x,drop=FALSE])%*%C.x[x,y,drop=FALSE]
      return(res)
  }

  ys <- data[,Y,drop=FALSE]
  ry <- t(ys)-Ey.x

  if (length(eta.idx)>0) {
    Ceta.x <- C.x[eta.idx,eta.idx]
    Lambda <- A[Y.idx.all,eta.idx,drop=FALSE] ##, ncol=length(eta.idx))
    Cetay.x <- Ceta.x%*%t(Lambda)
    KK <- Cetay.x %*% solve(Cy.x)
    Eeta.y <- Eeta.x + KK %*% ry

    Ceta.y <- Ceta.x - KK%*% t(Cetay.x)
  } else {
    Eeta.y <- NA
    Ceta.y <- NA
  }

  Vhat <- matrix(0, nrow(data), length(vars(object))); colnames(Vhat) <- vars(object)
  Vhat[,obs.idx] <- as.matrix(data[,manifest(object)])
  if (length(eta.idx)>0)
    Vhat[,latent(object)] <- t(Eeta.y)
  I <- diag(nrow=nrow(A));
  epsilonhat <- (t( IA%*%t(Vhat) - m$v ))[,c(endogenous(object),latent(object)),drop=FALSE]
  if (residual) {
    return(epsilonhat)
  }

  mydata <- matrix(0,ncol=ncol(A),nrow=nrow(data)); colnames(mydata) <- vars(object)
  mydata[,manifest(object)] <- as.matrix(data[,manifest(object)])
  for (i in latent(object))
    mydata[,i] <- m$v[i]

  Yhat <- t(mydata%*%t(A)) + (m$v)
  res <- cbind(t(Ey.x)) ## Conditional mean

  attr(res, "cond.var") <- Cy.x
  attr(res, "blup") <- t(Eeta.y)
  attr(res, "var.blup") <- Ceta.y
  attr(res, "Ey.x") <- Ey.x
  attr(res, "eta.x") <- Eeta.x
  attr(res, "epsilon.y") <- epsilonhat
  attr(res, "p") <- m$p
  attr(res, "e") <- m$e
  class(res) <- c("lvm.predict","matrix")
  return(res)
}

##' @export
print.lvm.predict <- function(x,...) print(x[,])
