"expand" <-
function(x,k){
  xx <- x
  for ( i in 2:k) xx <- rbind(x,xx)
  xx}

gqz <- function(numnodes=20, minweight=0.000001){
#  Calculate Gaussian Quadrature points for the Normal distribution
#  using the abscissas and weights for Hermite integration. The
#  conversion of the locations and weights is given in Lindsey (1992,
#  page 169:3) and Skrondal & Rabe-Hesketh (2004, page 165:1).
#  The argument numnodes is the theoretical number of quadrature points,
#  locations with weights that are less than the argument minweight will
#  be omitted.
#  The default vale of minweight=0.000001 returns 14 masspoints for the
#  default numnodes=20 as in Aitkin, Francis & Hinde (2005).
    out <- gauss.quad(numnodes, "hermite")
    h <- rbind(out$nodes*sqrt(2), out$weights/sum(out$weights))
#  Sort the locations and weights into columns in decending order of the
#  location vector.
    ord<-order(h[1,], decreasing = TRUE)
    h <- h[,ord]
    h <- cbind(h[1,], h[2,])
    h <- subset(as.data.frame(h), (h[,2] >= minweight))
    names(h) <- c("location","weight")
    h
    }
gauss.quad <- function (n, kind = "legendre", alpha = 0, beta = 0)
{
    n <- as.integer(n)
    if (n < 0)
        stop("need non-negative number of nodes")
    if (n == 0)
        return(list(nodes = numeric(0), weights = numeric(0)))
    kind <- match.arg(kind, c("legendre", "chebyshev1", "chebyshev2",
        "hermite", "jacobi", "laguerre"))
    i <- 1:n
    i1 <- i[-n]
    switch(kind, legendre = {
        muzero <- 2
        a <- rep(0, n)
        b <- i1/sqrt(4 * i1^2 - 1)
    }, chebyshev1 = {
        muzero <- pi
        a <- rep(0, n)
        b <- rep(0.5, n - 1)
        b[1] <- sqrt(0.5)
    }, chebyshev2 = {
        muzero <- pi/2
        a <- rep(0, n)
        b <- rep(0.5, n - 1)
    }, hermite = {
        muzero <- sqrt(pi)
        a <- rep(0, n)
        b <- sqrt(i1/2)
    }, jacobi = {
        ab <- alpha + beta
        muzero <- 2^(ab + 1) * gamma(alpha + 1) * gamma(beta +
            1)/gamma(ab + 2)
        a <- i
        a[1] <- (beta - alpha)/(ab + 2)
        i2 <- 2:n
        abi <- ab + 2 * i2
        a[i2] <- (beta^2 - alpha^2)/(abi - 2)/abi
        b <- i1
        b[1] <- sqrt(4 * (alpha + 1) * (beta + 1)/(ab + 2)^2/(ab +
            3))
        i2 <- i1[-1]
        abi <- ab + 2 * i2
        b[i2] <- sqrt(4 * i2 * (i2 + alpha) * (i2 + beta) * (i2 +
            ab)/(abi^2 - 1)/abi^2)
    }, laguerre = {
        a <- 2 * i - 1 + alpha
        b <- sqrt(i1 * (i1 + alpha))
        muzero <- gamma(alpha + 1)
    })
    A <- rep(0, n * n)
    A[(n + 1) * (i - 1) + 1] <- a
    A[(n + 1) * (i1 - 1) + 2] <- b
    A[(n + 1) * i1] <- b
    dim(A) <- c(n, n)
    vd <- eigen(A, symmetric = TRUE)
    w <- rev(as.vector(vd$vectors[1, ]))
    w <- muzero * w^2
    x <- rev(vd$values)
    list(nodes = x, weights = w)
}

"weightslogl.calc.w" <-
function(p,fjk,weights){
# amended to rescale posterior probabilities
# R.E.D. 15th Feb 2006
# p is a vector of length K containing the mixture proportions
# fjk is a JXK matrix of log densities
logpf <- t(apply(fjk,1,"+",log(p)))
Mi <- apply(logpf,1,max)
Sik <- logpf-Mi
ifelse(Sik < -760, 0, Sik)
eSik <- exp(Sik)
SeSik <- as.vector(apply(eSik,1,sum))
w <- eSik/SeSik
ML.dev <- -2*sum(weights*log(apply(eSik,1,sum)))-2*sum(weights*Mi)
list(w=w, ML.dev=ML.dev)
}

"summary.pattNPML" <-
function(object,digits=max(3,getOption('digits')-3), ...){
  np   <-  length(object$coefficients)
  if (np > 0){
      m <- seq(1,np)[substr(attr(object$coefficients,'names'),1,4)=='MASS']
      mass.points <- object$coefficients[m]
      cat('\nCall: ',deparse(object$call),'\n', fill=TRUE)
      cat('Coefficients')
      cat(":\n")
      df.r <- object$lastglm$df.residual

      glm.dispersion <- if (any(object$family$family == c("poisson", "binomial")))
            1
            else if (df.r > 0) {
                sum(object$lastglm$weights * object$lastglm$residuals^2, na.rm=TRUE)/df.r
            }
            else Inf
      lastglmsumm <- summary.glm(object$lastglm, dispersion=glm.dispersion)
      fitcoef     <- matrix(lastglmsumm$coeff[,1:3], ncol=3,dimnames= list(dimnames(lastglmsumm$coef)[[1]], c(dimnames(lastglmsumm$coeff)[[2]][1:2], "t value") ))  #03-08-06
      print(fitcoef)
  } else {
      cat('\nCall: ',deparse(object$call),'\n', fill=TRUE)
      cat("No coefficients. \n")
  }

  p <- object$masses
  #names(p) <- paste('MASS',seq(1,ncol(object$post.prob)),sep='') # omitted from 0.42 on
  dispersion <- 1    # now calculate dispersion in the sense of 'dispersion parameter':

  cat('\nMixture proportions:\n')
  print.default(format(p,digits),print.gap=2,quote=FALSE)

  if (object$family$family=='gaussian'){
      dispersion <- (object$sdev$sdev)^2
      if (object$Misc$lambda<=1/length(object$masses)){
          cat('\nComponent distribution - MLE of sigma:\t  ',format(signif(object$sdev$sdev,digits)),'\n')
      } else {
          cat('\nMLE of component standard deviations:\n'); s<- object$sdev$sdevk; names(s)<- names(p);  print.default(format(s,digits),print.gap=2,quote=FALSE); cat ('\n')
      }
  } else if (object$family$family=='Gamma'){
      dispersion <- 1/object$shape$shape
      if (object$Misc$lambda<=1/length(object$masses)){cat('\nComponent distribution - MLE of shape parameter:\t  ',format(signif(object$shape$shape,digits)),'\n')
      } else {
          cat('\n MLE of component shape parameters:\n'); s<- object$shape$shapek; names(s)<- names(p);  print.default(format(s,digits),print.gap=2,quote=FALSE); cat ('\n')
      }
  } else cat('\n')

###  cat('Random effect distribution - standard deviation:\t  ', format(object$rsdev),'\n\n') # 03/09/07
##cat('\n')

  cat('Deviance:\t   ',format(round(object$deviance,digits=2)),"\n")
  cat('-2 log L:\t   ',format(round(object$disparity,digits=2)))
  conv <- object$EMconverged
  if(conv != "none"){
      cnv <- ifelse(conv,"","no")
      if (!is.null(object$post.prob)) cat('   ',cnv,'EM convergence at iteration ',round(object$EMiter,0))
  }
  cat('\n')

  if (np >0){
      invisible(c("coefficients"=list(fitcoef), object[-1],list(dispersion=dispersion), list(lastglmsumm=lastglmsumm)))
  } else {
      invisible(c("coefficients"=list(fitcoef), object[-1],list(dispersion=dispersion)) )
  }
}

"print.pattNPML" <-
function(x,digits=max(3,getOption('digits')-3), ...){
  np <- length(x$coefficients)
  # print(names(x))
  if (np > 0){
    m <- seq(1,np)[substr(attr(x$coefficients,'names'),1,4)=='MASS']
    mass.points <- x$coefficients[m]
    cat('\nCall: ',deparse(x$call),'\n', fill=TRUE)
    cat('Coefficients')
    cat(":\n")
    print.default(format(x$coefficients[1:np],digits = digits), print.gap = 2,quote = FALSE)#;cat('\n')
  } else {
    cat('\nCall: ',deparse(x$call),'\n', fill=TRUE)
    cat("No coefficients. \n")
  }

  if (x$family$family=='gaussian' && x$Misc$lambda<=1/length(x$masses)){ # print sigma only if it is constant over components
        cat('Component distribution - MLE of sigma:\t  ',                # 03/09/07
        format(signif(x$sdev$sdev,digits)),'\n')
    }
  if (x$family$family=='Gamma'&& x$Misc$lambda<=1/length(x$masses)){ # print shape only if it is constant over components
       cat('Component distribution - MLE of shape parameter:\t  ',format(signif(x$shape$shape,digits)),'\n')# 03/09/07
    }

###  cat('Random effect distribution - standard deviation:\t  ', format(x$rsdev),'\n\n') # 03/09/07
cat('\n')

  if (!is.null(x$post.prob)){
      p <- x$masses
      # names(p) <- paste('MASS',seq(1,ncol(x$post.prob)),sep='') # omitted from 0.42
      cat('Mixture proportions')
      cat(":\n")
      print.default(format(p,digits),print.gap=2,quote=FALSE)
      cat('\n')
   }
  cat('Deviance:\t   ',format(round(x$deviance,digits=2)),"\n")
  cat('-2 log L:\t   ',format(round(x$disparity,digits=2)))
  conv <- x$EMconverged
  if(conv != "none"){
      cnv <- ifelse(conv,"","no")
      if (!is.null(x$post.prob)) cat('   ',cnv,'EM convergence at iteration ',round(x$EMiter,0))
  }
  cat('\n')
  invisible(x)
}
# obsolete
#"BIC.pattNPML"<-
#function(object, ...){
#object$disparity + (length(c(na.omit(object$coefficients),object$masses))-1) * log(sum(object$data$y))
#}
# Achim's proposal - thanks Achim
# allows for using stats AIC(), BIC()
logLik.pattNPML <- function(object, ...)
    structure(-0.5 * object$disparity,
      df = length(c(na.omit(object$coefficients), object$masses)) - 1,
      class = "logLik")

nobs.pattNPML <- function(object, ...) sum(object$data$y)
