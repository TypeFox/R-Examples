coxseifit.ex <- function(dat,par.init,m=2, mit=1000,tr=TRUE,method="L-BFGS-B",
                         lower=c(rep(-Inf,ncol(dat)-3),-Inf,0),
                         upper=rep(Inf,ncol(dat)-3 + 2),...){
  ids <- unique(dat$id)
  ng <- length(ids)
  gs <- as.numeric(table(dat$id))
  gofst <- cumsum(gs) - gs
  C <- dat$Y[!dat$delta]
  Z <- as.matrix(dat[, setdiff(colnames(dat), c("Y", "delta", 
                                                "id"))], nrow = nrow(dat))
  Zs <- array(0, dim = c(dim(Z), ng))
  for (l in 1:ng) {
    for (i in 1:ng) {
      for (j in 1:NCOL(Z)) {
        Zs[dat$id == ids[i], j, l] <- if (diff(range(Z[dat$id == 
             ids[l], j])) == 0) {
          Z[dat$id == ids[l], j][1]
        }
        else {
          approxfun(dat$Y[dat$id == ids[l]], Z[dat$id == 
                              ids[l], j], method = "constant", f = 1, rule = 2)(dat$Y[dat$id == 
                                                                        ids[i]])
        }
      }
    }
  }
  loglik <- function(param) {
    .Call("ll2", dat$Y, as.double(as.vector(Z)), as.double(as.vector(Zs)), 
          as.double(C), as.integer(gs), as.integer(gofst), 
          as.double(param), as.integer(m), PACKAGE = "coxsei")
  }
  ##     if(pmatch(exf,c(
  ##                   }
  ret <- optim(par.init, loglik, control = list(trace = tr, maxit = mit),
               hessian = TRUE, method = method, lower=lower, upper=upper,...)
  allpar <- ret$par;
  np <- length(allpar);
  par <- allpar[1:(np-2)];
##   parg <- exp(allpar[np-1:0]);
  parg <- allpar[np-1:0];
  ##   drv <- diag(c(rep(1,np-2),parg));
  ##   vcovmat <- drv %*% solve(ret$hessian, drv);
  vcovmat <- solve(ret$hessian)
  gfun <- function(x, pa) {
    ifelse(x <= 0, 0, pa[1] * exp(-pa[2] * x))
  }
  gfungrd <- function(x, pa) {
    if (length(x) == 0) 
      return(matrix(0, 2, 0))
    rbind(exp(-pa[2] * x), - x* pa[1] * exp(-pa[2] * x))
  }
  js <- function(i) {
    res <- numeric(gs[i] - 1)
    res1 <- matrix(0, nrow = np, ncol = gs[i] - 1)
    for (j in 1:(gs[i] - 1)) {
      posij <- gofst[i] + j
      for (l in (1:ng)[dat$Y[posij] <= C]) {
        yl <- dat$Y[dat$id == l & (dat$delta == 1)]
        R0part <- exp(Zs[posij, , l] %*% par +
                      sum(gfun(dat$Y[posij] -
                               tail(yl[yl < dat$Y[posij]], m), parg)))
        res[j] <- res[j] + R0part
        res1[, j] <- res1[, j] +
          c(Zs[posij, , l],
            rowSums(gfungrd(dat$Y[posij] -
                            tail(yl[yl < dat$Y[posij]],
                                 m), parg))) * R0part
      }
    }
    res1 <- apply(res1, 1, "/", y = res^2)
    dim(res1) <- c(gs[i] - 1, np)
    res1 <- t(res1)
    rbind(1/res, res1)
  }
  x <- dat$Y[dat$delta == 1]
  ox <- order(x)
  jmps <- unlist(lapply((1:ng)[gs > 1], js))
  lj <- length(jmps)
  jo <- jmps[seq(1, lj, by = 1 + np)][ox]
  varestjo <- jmps[-seq(1, lj, by = 1 + np)]
  dim(varestjo) <- c(np, length(x))
  varestjo <- varestjo[, ox]
  varest <- apply(varestjo, 1, cumsum)
  dim(varest) <- c(length(x), np)
  varest <- cumsum(jo^2) + apply(varest, 1, function(x) x %*% (vcovmat %*% x))
  res <- list();
  co <- c(par,parg);
  names(co) <- c(colnames(Z),"alpha","gamma")
  res$coefficients <- co
  res$vcov <- vcovmat
  res$zval <- co/sqrt(diag(vcovmat))
  res$pval <- c(pnorm(abs(res$zval[1:(np-2)]),lower.tail=FALSE)*2,
                pnorm(abs(res$zval[np-1:0]),lower.tail=FALSE))
  res$details.par <- ret
  int <- cumsum(jo);
  res$cintfn <- stepfun(x[ox],c(0,int),right=TRUE)
  res$cintvar <- stepfun(x[ox],c(0,varest),right=TRUE)
  res$details.cint <- list(times=x[ox],cint=int,cintvar=varest)
  res$call <- match.call()
  class(res) <- "coxsei"
  res
}

coxseiexp <- function(Y,delta,id,Z,par.init,m=2,mit=1000,tr=TRUE,
                      method="L-BFGS-B",
                      lower=c(rep(-Inf,ncol(Z)),-Inf,0),
                      upper=rep(Inf,ncol(Z) + 2),...){
  ids <- unique(id)
  ng <- length(ids)
  gs <- as.numeric(table(id))
  gofst <- cumsum(gs) - gs
  C <- Y[!delta]
  Z <- as.matrix(Z)
  Zs <- array(0, dim = c(dim(Z), ng))
  for (l in 1:ng) {
    for (i in 1:ng) {
      for (j in 1:NCOL(Z)) {
        Zs[id == ids[i], j, l] <- if (diff(range(Z[id == 
             ids[l], j])) == 0) {
          Z[id == ids[l], j][1]
        }
        else {
          approxfun(Y[id == ids[l]], Z[id == 
                          ids[l], j], method = "constant", f = 1, rule = 2)(Y[id == 
                                                                    ids[i]])
        }
      }
    }
  }
  loglik <- function(param) {
    .Call("ll2", Y, as.double(as.vector(Z)), as.double(as.vector(Zs)), 
          as.double(C), as.integer(gs), as.integer(gofst), 
          as.double(param), as.integer(m), PACKAGE = "coxsei")
  }
  ret <- optim(par.init, loglik, control = list(trace = tr, maxit = mit),
               hessian = TRUE, method = method, lower=lower,
               upper = upper,...)
  allpar <- ret$par;
  np <- length(allpar);
  par <- allpar[1:(np-2)];
##   parg <- exp(allpar[np-1:0]);
  parg <- allpar[np-1:0];
  ## drv <- diag(c(rep(1,np-2),parg));
##   vcovmat <- if (m>0)drv %*% solve(ret$hessian, drv) else{
##     tmp <- matrix(0,np,np);
##     tmp[1:(np-2),1:(np-2)] <- solve(ret$hessian[1:(np-2),1:(np-2)]);
##     tmp
##   }
  vcovmat <- if (m>0)solve(ret$hessian) else{
    tmp <- matrix(0,np,np);
    tmp[1:(np-2),1:(np-2)] <- solve(ret$hessian[1:(np-2),1:(np-2)]);
    tmp
  }

  gfun <- function(x, pa) {
    ifelse(x <= 0, 0, pa[1] * exp(-pa[2] * x))
  }
  gfungrd <- function(x, pa) {
    if (length(x) == 0) 
      return(matrix(0, 2, 0))
    rbind(exp(-pa[2] * x), - x * pa[1] * exp(-pa[2] * x))
  }
  js <- function(i) {
    res <- numeric(gs[i] - 1)
    res1 <- matrix(0, nrow = np, ncol = gs[i] - 1)
    for (j in 1:(gs[i] - 1)) {
      posij <- gofst[i] + j
      for (l in (1:ng)[Y[posij] <= C]) {
        yl <- Y[id == l & (delta == 1)]
        R0part <- exp(Zs[posij, , l] %*% par +
                      sum(gfun(Y[posij] -
                               tail(yl[yl < Y[posij]], m), parg)))
        res[j] <- res[j] + R0part
        res1[, j] <- res1[, j] +
          c(Zs[posij, , l],
            rowSums(gfungrd(Y[posij] -
                            tail(yl[yl < Y[posij]],
                                 m), parg))) * R0part
      }
    }
    res1 <- apply(res1, 1, "/", y = res^2)
    dim(res1) <- c(gs[i] - 1, np)
    res1 <- t(res1)
    rbind(1/res, res1)
  }
  x <- Y[delta == 1]
  ox <- order(x)
  jmps <- unlist(lapply((1:ng)[gs > 1], js))
  lj <- length(jmps)
  jo <- jmps[seq(1, lj, by = 1 + np)][ox]
  varestjo <- jmps[-seq(1, lj, by = 1 + np)]
  dim(varestjo) <- c(np, length(x))
  varestjo <- varestjo[, ox]
  varest <- apply(varestjo, 1, cumsum)
  dim(varest) <- c(length(x), np)
  varest <- cumsum(jo^2) + apply(varest, 1, function(x) x %*% (vcovmat %*% x))
  res <- list();
  co <- c(par,parg);
  names(co) <- c(colnames(Z),"alpha","gamma")
  res$coefficients <- co
  res$vcov <- vcovmat
  res$zval <- co/sqrt(diag(vcovmat))
  res$pval <- c(pnorm(abs(res$zval[1:(np-1)]),lower.tail=FALSE)*2,
                pnorm(abs(res$zval[np]),lower.tail=FALSE))
  res$details.par <- ret
  int <- cumsum(jo);
  res$call <- match.call()
  res$cintfn <- stepfun(x[ox],c(0,int),right=TRUE)
  res$cintvar <- stepfun(x[ox],c(0,varest),right=TRUE)
  res$details.cint <- list(times=x[ox],cint=int,cintvar=varest)
  class(res) <- "coxsei"
  res
}

coxsei <- function(x,...){UseMethod("coxsei")}
coxsei.default <- function(x,y,delta,id,par.init,m=2,mit=1000,tr=TRUE,
                           method="L-BFGS-B",
                           lower=c(rep(-Inf,ncol(x)),-Inf,0),
                           upper=rep(Inf,ncol(x) + 2),...){
  res <- coxseiexp(Y=y,delta=delta,id=id,Z=x,par.init=par.init,m=m,mit=mit,tr=tr,method=method,lower=lower,upper=upper,...)
  res$call <- match.call()
  class(res) <- "coxsei"
  res
}
print.coxsei <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoeffficients:\n")
  print(x$coefficients)
  cat("\nCumulative baseline intensity function:\n")
  print(x$cintfn)
}
plot.coxsei <- function(x,...){
  plot(x$cintfn,main="cumulative baseline intensity function",...)
}
summary.coxsei <- function(object,...){
  np <- length(coef(object));
  tab <- cbind(Estimate=coef(object),
               StdErr=sqrt(diag(object$vcov)),
               z.value=object$zval,
               p.value=object$pval)
  res <- list(call=object$call,coefficients=tab,cumint=object$cintfn);
  class(res) <- "summary.coxsei";
  res
}
print.summary.coxsei <- function(x,...){
  cat("Call:\n");
  print(x$call);
  cat("\n")
  printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE)
  cat("\nCumulative baseline intensity function:\n")
  print(x$cumint)
}
