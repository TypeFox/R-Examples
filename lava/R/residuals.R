Isqrt <- function(X) {
    eX <- eigen(X);
    with(eX, vectors %*% diag(1/sqrt(values),nrow=length(values)) %*% t(vectors))
}


##' @export
residuals.multigroupfit <- function(object,data=model.frame(object),p=coef(object), k, ...) {
  pp <- modelPar(object,p,...)
  if (!missing(k)) return(residuals(object$model$lvm[[k]],data=data[[k]],p=pp$p[[k]],...))
  res <- c()
  for (i in seq(length(pp$p))) {
    res <- c(res, list(residuals(object$model$lvm[[i]],data=data[[i]],p=pp$p[[i]],...)))
  }
  return(res)
}


##' @export
residuals.lvmfit <- function(object,data=model.frame(object),p=coef(object),...) {
  residuals(Model(object), data=data, p=p, ...)
}

##' @export
residuals.lvm <- function(object,data=model.frame(object),std=FALSE,p=coef(object),...) {
  Y <- setdiff(manifest(object), X <- exogenous(object))
  Pr <- predict(object,p=p,data=data)
  
  PrY <- Pr[,Y,drop=FALSE]
  ##  y <- endogenous(object)[match(endogenous(object),manifest(object))]
  r <- as.matrix(data[,Y,drop=FALSE]-(PrY))
  res <- r

  if (std) {
    S <- attributes(Pr)$cond.var;
    if (length(Y)>1) {
      res <- r%*%Isqrt(S)
    } else res <- 1/sqrt(S[1,1])*r
  }
  colnames(res) <- colnames(r)
  res
}

gradpredict <- function(p,obj,data=model.frame(obj)) {
##  res <- residuals.lvmfit(object=obj,data=data,std=FALSE,p=p)
  mom <- moments(Model(obj),p,conditional=TRUE)
  D <- with(obj, deriv.lvm(model, meanpar=modelPar(model,p)$meanpar, mom=mom))

  X <- with(obj, exogenous(model))
  Y <- with(obj, setdiff(manifest(model), X))
  X.idx <- with(obj, match(X,manifest(model)))
  X.idx.all <- with(obj, match(X, vars(model)))

  mu0 <- mom$v; mu0[X.idx.all] <- 0
  xs <- data[,X,drop=FALSE]
  mu.x <- apply(xs, 1, FUN=function(i) { res  <- rep(0,length(mu0)); res[X.idx.all] <- i; res })
  ##xi.x <- (mu.0 + mu.x)
  mu.0 <- rbind(rep(1,nrow(data))) %x% mu0
  mu. <- mu.0 + mu.x

  K <- nrow(mom$J)
  I <- diag(nrow=K)

  d1 <- (t(mu.) %x% I)%*%D$dvecG
  px <- index(obj)$px
  G <- mom$G ## index(obj)$Jy %*% solve(diag(nrow(mom$A))-t(mom$A))
  d2 <-  cbind(rep(1,ncol(mu.))) %x% (G%*% px%*% D$dvecv)
  dvecmui = d1+d2
  ## Reorganize vector
  idx <- (0:(nrow(data)-1))*K+1
  idx. <- c(); for (i in seq_len(K)) idx. <- c(idx., idx+i-1)
  dvecmui <- dvecmui[idx.,]
  ## first n rows first coordinate of predicted values,
  ## next n rows second coordinates and so forth
  dvecmui
}
