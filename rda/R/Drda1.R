Drda1 <- function(x, y, xnew=NULL, ynew=NULL, prior=table(y)/length(y),
                  alpha=seq(0.01, 1, len=10), delta=seq(0, 3, len=10),
                  genelist=FALSE, trace=FALSE){
  this.call <- match.call()

  testdata.x <- (is.null(xnew) || missing(xnew))
  testdata.y <- (is.null(ynew) || missing(ynew))

  mu <- apply(x, 1, mean) 
  x <- x-mu              

  yhat <- array(0, dim=c(length(alpha), length(delta), length(y)))

  if (!(testdata.x)){
    xnew <- xnew-mu
    yhat.new <- array(0, dim=c(length(alpha), length(delta), 
                               ncol(xnew)))
    posterior <- array(0, dim=c(length(alpha), length(delta), 
                                ncol(xnew), length(prior)))
  }

  if (genelist){
    gene.list <- array(0, dim=c(length(alpha), length(delta), nrow(x)))
    dn2 <- list(round(alpha, 3), round(delta, 3), seq(nrow(x)))
    names(dn2) <- c("alpha", "delta", "gene #")
    dimnames(gene.list) <- dn2
  }

  Y <- model.matrix( ~ factor(y)-1)    
  dimnames(Y) <- list(NULL, NULL)

  xbar <- scale(x%*%Y, FALSE, table(y)) 

  xstar <- x-xbar[, unclass(factor(y))]

  DD <- sqrt(apply(xstar^2, 1, sum)/length(y))

  xstar1 <- xstar/DD
  xout <- t(xstar1)%*%xstar1
  sout <- svd(xout)

  svalue <- sout$d
  svaluePos <- seq(svalue)[svalue > 0] 
  svalue <- sqrt(svalue[svaluePos])

  V <- scale(xstar1%*%sout$u[, svaluePos], FALSE, svalue)
  svalue <- svalue/sqrt(length(y))

  xbar1 <- xbar/DD
  Vtxbar <- t(V)%*%xbar1

  error <- ngene <- testerror <- matrix(0, length(alpha), length(delta))
  if (!(testdata.y)){
    testerror <- error
  }

  for(i in seq(along=alpha)){
    dd <- 1/(svalue^2*alpha[i]+(1-alpha[i]))-1/(1-alpha[i])
    temp <- Vtxbar*dd

    coefmat <- (V%*%temp+xbar1/(1-alpha[i]))/DD
    for(j in seq(along=delta)){
      coefmat1 <- I(abs(coefmat)>delta[j])*(abs(coefmat)-
                    delta[j])*sign(coefmat)

      if(genelist){
        gene.list[i, j, ] <- I(apply(I(abs(coefmat1) > 0), 1, sum) > 0)
      }

      ngene[i, j] <- sum(apply(I(abs(coefmat1) > 0), 1, sum) > 0)
      

      coefmat2 <- coefmat1*DD
      ystar <- t(V)%*%coefmat2
      intvec <- -(alpha[i]*apply(ystar^2*svalue^2, 2, sum)+
               (1-alpha[i])*apply(coefmat2^2, 2, sum))/2+log(prior)

      dhat <- scale(t(x) %*% coefmat1, -intvec, FALSE)
      dimnames(dhat) <- list(NULL, names(prior))
      yhat[i, j, ] <- softmax(dhat)

      error[i, j] <- sum(y != yhat[i, j, ])

      if(!(testdata.x)){
        dhat <- scale(t(xnew)%*%coefmat1, -intvec, FALSE)
        posterior[i, j, , ] <- dhat
        dimnames(dhat) <- list(NULL, names(prior))
        yhat.new[i, j, ] <- softmax(dhat)
        if(!(testdata.y)){
          testerror[i, j] <- sum(ynew != yhat.new[i, j, ])
        }
      }
 
      if(trace){
        if(!(testdata.x || testdata.y)){
          cat("alpha=", alpha[i], "delta=", delta[j], "test error=",
              format(round(testerror[i, j], 3)), "\n")
        }
        else{
          cat("alpha=", alpha[i], "delta=", delta[j], "error=",
              format(round(error[i, j], 3)), "\n")
        }
      }
    }
  }
  dn <- list(round(alpha, 3), round(delta, 3))
  names(dn) <- c("alpha", "delta")
  dimnames(error) <- dimnames(ngene) <- dn

  dn3 <- list(round(alpha, 3), round(delta, 3), seq(ncol(x)))
  names(dn3) <- c("alpha", "delta", "sample")
  dimnames(yhat) <- dn3

  obj <- list(alpha = alpha, delta = delta, prior=prior,
              error=error, ngene=ngene, centroids=xbar,
              centroid.overall=mu, call=this.call, yhat=yhat)
  if(!(testdata.x)) {
    dn4 <- list(round(alpha, 3), round(delta, 3), seq(ncol(xnew)))
    names(dn4) <- c("alpha", "delta", "sample")
    dimnames(yhat.new) <- dn4
    obj$yhat.new <- yhat.new

    dn1 <- list(round(alpha, 3), round(delta, 3), seq(ncol(xnew)),
                names(prior))
    names(dn1) <- c("alpha", "delta", "sample #", "class")
    dimnames(posterior) <- dn1
    obj$posterior <- posterior

    if(!(testdata.y)) {
      dimnames(testerror) <- dn
      obj$testerror <- testerror
    }
  }
  if(genelist) {
    obj$gene.list <- gene.list
  }
  obj
}    

