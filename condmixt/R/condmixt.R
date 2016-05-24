# library of functions for unconditional and conditional hybrid Pareto mixture model and variants

softplus <- function(x){
  y <- rep(NaN,length(x))
  if(any(x>0)){
    y[x>0] <- log(1+exp(-x[x>0]))
    if (any(is.finite(y)))
      y[is.finite(y)] <- y[is.finite(y)]+x[is.finite(y)]
     }
    if(any(x<=0)){
      y[x<=0]=log(1+exp(x[x<=0]))
    if (any(!is.finite(y)))
      y[!is.finite(y)] <- 0
    }
    y
}

softplusinv <- function(y){
  log(exp(y)-1)
}

lambertw <- function(z){
  w0 <- .C("lambertwR", as.double(z), results=double(1),PACKAGE="condmixt")
  w0[["results"]]
}

hpareto.nll <- function(theta,x){
  # assume that the 3rd parameter is encoded as log(sigma)
  if (length(theta)!=3)
    stop("THETA must have length 3")
  n <- length(x)
  if (n == 0)
    stop("X must have a positive length")
  results <- .C("hpnll",as.double(theta),as.double(x),as.integer(n),nll=double(1),nllgrad=double(3),PACKAGE="condmixt")

  obj <- results[["nll"]]
  attr(obj,"gradient") <-  results[["nllgrad"]]
  obj
}


hpareto.fit <- function(params,x,...){
  theta <- params
  theta[3] <- log(params[3])
  opt <- nlm(hpareto.nll,theta,x,...)
  params.opt <- opt$estimate
  params.opt[3] <- exp(params.opt[3])
  params.opt
}


hpareto.negloglike <- function(params,x){
  # params is a vector of the three hybrid Pareto parameters xi, mu and sigma
  theta <- params
  theta[3] <- log(params[3])
  hpareto.nll(theta,x)
}

hillest <- function(data,k){
  # Computes the Hill estimator of the tail index xi
  # The threshold used is the K+1 order statistic

  xstat <- sort(data, decreasing = TRUE)
  xihat <- mean(log(xstat[1:k])) - log(xstat[k+1])
  xihat
}


gpd.mme <- function(x){
  #  Moments estimator for the GPD
  xbar <- mean(x)
  s2 <- var(x)
  s2 <- s2 + (s2==0)
  xi0 <- -0.5 * (((xbar * xbar)/s2) - 1)
  beta0 <- 0.5 * xbar * (((xbar * xbar)/s2) + 1)
  c(xi0,beta0)
}


hpareto.mme <- function(x,xi0=c(),p=0.99){
  # equivalent to sgpdmmepos.m, initialises the parameters of a hybrid Pareto mixture based on the method of moments and the data in X, no truncation here
  # more doc to be found in that file /Users/julie/diro/evt/Mfiles/mystats/sgpdmmepos.m

  k <- length(x) - round(p*length(x)) + 1  # use the 10% largest data to estimate xi

  npos <- sum(x > 0) # ensure the k+1 largest values are positive
  if (npos < k + 1)
    k <- npos -1

  if (length(xi0)==0){ # xi0 is not given externally
    xi0 <- hillest(x, k)
    if (xi0 < 0) # enforce positivity of the tail index
      xi0 <- 0
  }

  z <- (1 + xi0)^2/(2 * pi)
  w <- lambertw(z)
  phiw <- pnorm(sign(1 + xi0) * sqrt(w))

  mu0 <- quantile(x, 1/(2 * (1 + phiw)))
  alpha0 <- quantile(x, phiw/(1 + phiw))
  names(mu0) <- NULL
  names(alpha0) <- NULL
  sigma0 <- (alpha0 - mu0)/ (sqrt(w) * sign(1 + xi0))

  if (sigma0 <= 0)
    if (sum(x > alpha0) == 0)
      sigma0 <- 1
    else{
      theta <- gpd.mme(x[x > alpha0] - alpha0)
      sigma0 <- sqrt(w) * theta[2]/sqrt((1+xi0)^2)
    }
  c(xi0,mu0,sigma0)
}


hpareto.beta <- function(xi,sigma=1){
  sigma * sqrt((1 + xi)^2) / sqrt(lambertw((1 + xi)^2 / (2 * pi)))
}


hpareto.alpha <- function(xi,mu=0,sigma=1){
  mu + sign(1 + xi) * sigma * sqrt(lambertw((1 + xi)^2 / (2 * pi)))
}

hpareto.gamma <- function(xi,mu=0,sigma=1,trunc=T){
  alpha <- hpareto.alpha(xi,mu,sigma)
  if (trunc){
    if (alpha > 0)
      return(1 + pnorm(sign(1 + xi) * sqrt(lambertw((1 + xi)^2 / (2 * pi))))-pnorm(-mu/sigma))
    else{
      beta <- hpareto.beta(xi,sigma)
      return(1-pgpd(-alpha,xi,beta))
    }
  } else
    return(1 + pnorm(sign(1 + xi) * sqrt(lambertw((1 + xi)^2 / (2 * pi)))))
}

dhpareto <- function(y,xi,mu=0,sigma=1,log=FALSE,trunc=TRUE){
  # hybrid Pareto density function
  w <- lambertw((1 + xi)^2 / (2 * pi))
  beta <- sigma * sqrt((1 + xi)^2) / sqrt(w)
  alpha <- sign(1 + xi) * sigma * sqrt(w) + mu

  p <- rep(0,length(y))
  if (trunc){
    if (alpha > 0){  # truncation includes the Gaussian part
      gamma <-  1 + pnorm(sign(1 + xi) * sqrt(lambertw((1 + xi)^2 / (2 * pi))))-pnorm(-mu/sigma)
      if (any(y>0 & y<=alpha))
        p[which(y>0 & y<=alpha)] <- dnorm(y[which(y>0 & y<=alpha)],mu,sigma,log)
    }
    else
      gamma <- 1-pgpd(0,alpha,beta,xi) # no points in the Gaussian part
    if (any(y>alpha & y > 0))
      p[which(y>alpha & y > 0)] <- dgpd(y[which(y>alpha & y > 0)],alpha,beta,xi,log)
  } else { # no truncation
    gamma <-  1 + pnorm(sign(1 + xi) * sqrt(lambertw((1 + xi)^2 / (2 * pi))))
    if (any(y<=alpha))
      p[which(y<=alpha)] <- dnorm(y[which(y<=alpha)],mu,sigma,log)
     if (any(y>alpha))
       p[which(y>alpha)] <- dgpd(y[which(y>alpha)],alpha,beta,xi,log)
  }
  if (log)
    p-log(gamma)
  else
    p/gamma
}


phpareto <- function(q,xi,mu=0,sigma=1,trunc=TRUE){
  # hybrid Pareto distribution function
  w <- lambertw((1 + xi)^2 / (2 * pi))
  beta <- sigma * sqrt((1 + xi)^2) / sqrt(w)
  alpha <- sign(1 + xi) * sigma * sqrt(w) + mu

  p <- rep(0,length(q))
  if (trunc){
    if (alpha > 0){ # truncation includes the Gaussian part
      gamma <-  1 + pnorm(sign(1 + xi) * sqrt(lambertw((1 + xi)^2 / (2 * pi))))-pnorm(-mu/sigma)
      if (any(q>0 & q<=alpha))
        p[which(q>0 & q<=alpha)] <- pnorm(q[which(q>0 & q<=alpha)],mu,sigma)-pnorm(0,mu,sigma)
      if (any(q>alpha))
        p[which(q>alpha)] <- pgpd(q[which(q>alpha)],alpha,beta,xi)+pnorm(alpha,mu,sigma)-pnorm(0,mu,sigma)
    } else{ # if alpha < 0 and trunc = T, then all the points are in the GPD part
      gamma <- 1-pgpd(0,alpha,beta,xi)
      if (any(q>alpha & q > 0))
        p[which(q>alpha & q > 0)] <- pgpd(q[which(q>alpha & q > 0)],alpha,beta,xi)-pgpd(0,alpha,beta,xi)
    }
  } else{ # no truncation
    gamma <-  1 + pnorm(sign(1 + xi) * sqrt(lambertw((1 + xi)^2 / (2 * pi))))
    if (any(q<=alpha))
      p[which(q<=alpha)] <- pnorm(q[which(q<=alpha)],mu,sigma)
    if (any(q>alpha))
      p[which(q>alpha)] <- pgpd(q[which(q>alpha)],alpha,beta,xi)+pnorm(alpha,mu,sigma)
  }

  p/gamma
}


rhpareto <- function(n,xi,mu=0,sigma=1,trunc=TRUE){
  w <- lambertw((1 + xi)^2 / (2 * pi))
  beta <- sigma * sqrt((1 + xi)^2) / sqrt(w)
  alpha <- sign(1 + xi) * sigma * sqrt(w) + mu

  u <- runif(n)
  r <- rep(NaN,n)

  if (trunc){
    if (alpha > 0){  # truncation includes the Gaussian part
      F0 <- pnorm(0,mu,sigma) # probability of being below zero for a Gaussian variable
      Falpha <- pnorm(alpha,mu,sigma) # probability of being below alpha for a Gaussian
      Halpha <- (Falpha-F0)/(1+Falpha-F0) # probability of being in the positive Gaussian part
      if (any(u<=Halpha)) # generate data in the positive Gaussian part
        r[which(u<=Halpha)] <- qnorm(runif(sum(u<=Halpha))*(Falpha-F0)+F0,mu,sigma)
      if (any(u>Halpha)) # generate data in the GPD part
        r[which(u>Halpha)] <- qgpd(runif(sum(u>Halpha)),alpha,beta,xi)
    } else { # no points in the Gaussian part
      Galpha <- pgpd(0,alpha,beta,xi) # prob of being in the negative part of the GPD
      r <- qgpd(u*(1-Galpha)+Galpha,alpha,beta,xi)
    }
  }
  else{  # no truncation
    Falpha <- pnorm(sign(1 + xi) * sqrt(lambertw((1 + xi)^2 / (2 * pi))))
    gamma <-  1 + Falpha
    Halpha <- Falpha/gamma   # probability of being in the Gaussian part
    if (any(u<=Halpha)) # generate data in the Gaussian part
      r[which(u<=Halpha)] <- qnorm(runif(sum(u<=Halpha))*Falpha,mu,sigma)
    if (any(u>Halpha)) # generate data in the GPD part
      r[which(u>Halpha)] <- qgpd(runif(sum(u>Halpha)),alpha,beta,xi)
  }

  r
}


qhpareto <- function(p,xi,mu=0,sigma=1,trunc=TRUE){
  w <- lambertw((1 + xi)^2 / (2 * pi))
  beta <- sigma * sqrt((1 + xi)^2) / sqrt(w)
  alpha <- sign(1 + xi) * sigma * sqrt(w) + mu
  Falpha <- pnorm((alpha-mu)/sigma)

  y <- rep(0,length(p))
  if (trunc)
    F0 <- pnorm(0,mu,sigma)
  else
    F0 <- 0

  Halpha <-  (Falpha-F0)/(1+Falpha-F0)
  if (any(p<=Halpha))
    y[p<=Halpha] <- qnorm((1+Falpha-F0)*p[p<=Halpha]+F0,mu,sigma)

  if (any(p>Halpha))
    y[p>Halpha] <- qgpd((1+Falpha-F0)*p[p>Halpha]-Falpha+F0,0,beta,xi)+alpha
  y
}


hparetomixt.init <- function(m,x,iter.max=20,nstart=10){
  # initialize a mixture of hybrid Paretos
  # equivalent to sgpdmminit.m but do not use a structure MIX, parameters of
  # the model are stored in params, 4 x m matrix

  ndata <- length(x) # assumes the data is univariate
  params.init <- matrix(nrow=4,ncol=m)

  if (m>1){
    clustering <- kmeans(x,m,iter.max,nstart)
    params.init[3,] <- clustering$centers

    #  Set priors depending on number of points in each cluster
    cluster.sizes <- sapply(clustering$size,max,1) # Make sure that no prior is zero
    params.init[1,] <- cluster.sizes/sum(cluster.sizes) # Normalise priors

     # Hybrid Pareto parameters are estimated on each cluster making sure that the xi (tail index parameters) are positive
    minel <- 10 # minimum number of points to perform density estimation on a component
    for (j in 1:m){
      if (clustering$size[j] < minel){ # not enough points in the cluster, initialise parameters randomly
        params.init[2,j] <- runif(1)*0.5 # initialise a small tail index (=thin tail)
        params.init[4,j] <- runif(1)*10+10
      } else {
        theta <- hpareto.mme(x[which(clustering$cluster==j)])
        params.init[2,j] <- theta[1]
        params.init[3,j] <- theta[2]
        params.init[4,j] <- theta[3]
      }

    }
  } else { # there is just one cluster
    theta <- hpareto.mme(x)
    params.init[1,1] <- 1
    params.init[2,1] <- theta[1]
    params.init[3,1] <- theta[2]
    params.init[4,1] <- theta[3]
  }
  if(any(params.init[4,]<0.1))
    params.init[4,params.init[4,]<0.1] <- 0.11
  if(any(params.init[2,]<0.001))
    params.init[2,params.init[2,]<0.001] <- 0.001
  params.init
}






dhparetomixt <- function(params,x,log=FALSE,trunc=TRUE){
  # params is 4 x m matrix
  m <- dim(params)[2]
  p <- rep(0,length(x))
  for (j in 1:m){
    p <- p+params[1,j]*dhpareto(x,params[2,j],params[3,j],params[4,j],log=F,trunc=FALSE)
  }
  if (trunc){
    F0 <- phparetomixt(params,0,trunc=FALSE)
    p <- p/(1-F0)
    if (any(x<=0))
      p[which(x<=0)] <- 0
  }

  if (log)
    p <- log(p)
  p
}

phparetomixt <- function(params,x,trunc=TRUE){
  # params is 4 x m matrix
  m <- dim(params)[2]
  p <- rep(0,length(x))
  for (j in 1:m){
    p <- p+params[1,j]*phpareto(x,params[2,j],params[3,j],params[4,j],trunc=FALSE)
  }
  if (trunc){
    F0 <- phparetomixt(params,0,trunc=FALSE)
    p <- p/(1-F0)
    if (any(x<=0))
      p[which(x<=0)] <- 0
  }
  p
}


hparetomixt.nll <- function(theta,m,x){
  results <- .C("ummhnll",as.double(theta),as.integer(m),as.double(x),as.integer(length(x)),nll=double(1),nllgrad=double(length(theta)),PACKAGE="condmixt")

  obj <- results[["nll"]]
  attr(obj,"gradient") <-  results[["nllgrad"]]
  obj
}

hparetomixt.negloglike <- function(params,x){
  # params is an 4 by m matrix
  m <- dim(params)[2]
  theta <- ummhbwd(params,m)
  hparetomixt.nll(theta,m,x)
}

ummhbwd <- function(params,m){
  cfunction <- .C("ummhbwd",as.double(t(params)), as.integer(m), theta=double(4*m-1),PACKAGE="condmixt")
  theta <- cfunction[["theta"]]
  if (any(is.infinite(theta)))
      theta[is.infinite(theta)] <- -1000
  theta
}

hparetomixt.fwd <- function(theta,m){
  cfunction <- .C("ummhfwd",as.double(theta), as.integer(m), params=double(4*m),PACKAGE="condmixt")
  matrix(cfunction[["params"]],4,m,byrow=TRUE)
}

hparetomixt.fit <- function(params,x,...){
  # params0 is an 4 by m matrix
  m <- dim(params)[2]
  theta0 <- ummhbwd(params,m)
  opt <- nlm(hparetomixt.nll,theta0,m,x,...)
  hparetomixt.fwd(opt$est,m)
}


hparetomixt.cvtrain <- function(m,x,nfold=5,nstart=1,...){
  n <- length(x)
  x <- x[sample(1:n,n)]  # shuffle data


  nk <- floor(n/nfold)  # number of observations per fold
  nlltest <- 0
  for (k in 1:nfold){ # train and test on each fold
    if (k < nfold)
      indk <- (1+(k-1)*nk):(k*nk)
    else
      indk <-  (1+(k-1)*nk):max(k*nk,n)
    cat("Training fold ", k ,"\t")
    nllbest <- Inf
    for (i in 1:nstart){
      params.init <- hparetomixt.init(m,x[-indk])
      theta <- ummhbwd(params.init,m)
      opt <- nlm(hparetomixt.nll,theta,m,x[-indk],...)
      if (opt$min<nllbest){
        nllbest <- opt$min
        thetabest <- opt$est
      } # if
    } # for nstart
    nlltest <- nlltest+hparetomixt.nll(thetabest,m,x[indk])
    cat("Test error is :",nlltest,"\n")
  }  # for nfold
  attr(nlltest,"gradient") <- NULL
  nlltest
}


hparetomixt.negloglike.tailpen <- function(params,lambda,w,beta,mu,sigma,x){
  m <- dim(params)[2]
  theta <- ummhbwd(params,m)
  hparetomixt.nll.bimodal.tailpen(theta,m,lambda,w,beta,mu,sigma,x)
}

hparetomixt.nll.bimodal.tailpen <- function(theta,m,lambda,w,beta,mu,sigma,x){
  results <- .C("ummhnll_bimodal_tailpen",as.double(theta),as.integer(m),as.double(lambda),
                as.double(w),as.double(beta),as.double(mu),as.double(sigma),as.double(x),
                as.integer(length(x)),nll=double(1),nllgrad=double(length(theta)),PACKAGE="condmixt")
  obj <- results[["nll"]]
  attr(obj,"gradient") <-  results[["nllgrad"]]
  obj
}

hparetomixt.fit.tailpen <- function(params,lambda,w,beta,mu,sigma,x,...){
  # params0 is an 4 by m matrix
  m <- dim(params)[2]
  theta0 <- ummhbwd(params,m)
  opt <- nlm(hparetomixt.nll.bimodal.tailpen,theta0,m,lambda,w,beta,mu,sigma,x,...)
  hparetomixt.fwd(opt$estimate,m)
}


hparetomixt.cvtrain.tailpen <- function(m,lambda,w,beta,mu,sigma,x,nfold=5,nstart=1,...){
  n <- length(x)
  nk <- floor(n/nfold)  # number of observations per fold
  x <- x[sample(1:n,n)]  # shuffle data

  nlltest <- 0
  for (k in 1:nfold){ # train and test on each fold
    if (k < nfold)
      indk <- (1+(k-1)*nk):(k*nk)
    else
      indk <-  (1+(k-1)*nk):max(k*nk,n)
    cat("Training fold ", k ,"\t")
    nllbest <- Inf
    for (i in 1:nstart){
      params.init <- hparetomixt.init(m,x[-indk])
      theta <- ummhbwd(params.init,m)
      opt <- nlm(hparetomixt.nll.bimodal.tailpen,theta,m,lambda,w,beta,mu,sigma,x[-indk],...)
      if (opt$min<nllbest){
        nllbest <- opt$min
        thetabest <- opt$est
      } # if
    } # for nstart
    nlltest <- nlltest+hparetomixt.nll(thetabest,m,x[indk])

    cat("Test error is :",nlltest,"\n")

  }
  attr(nlltest,"gradient") <- NULL
  nlltest
}

hparetomixt.disp <- function(params){
  # params0 is an 4 by m matrix
  m <- dim(params)[2]
  for (j in 1:m){
    cat(paste(j, ": &", format(params[1,j],digits=4), " &", format(params[2,j],digits=4), " &", format(params[3,j],digits=4), " &", format(params[4,j],digits=4), "\\\\", "\n"))
  }
}


gaussmixt.init <- function(m,x,iter.max=20,nstart=10){
  # initialize a mixture of Gaussians
  # the model are stored in params, m x 3 matrix

  ndata <- length(x) # assumes the data is univariate
  params.init <- matrix(nrow=3,ncol=m)

  clustering <- kmeans(x,m,iter.max,nstart)
  params.init[2,] <- clustering$centers

  #  Set priors depending on number of points in each cluster
  cluster.sizes <- sapply(clustering$size,max,1) # Make sure that no prior is zero
  params.init[1,] <- cluster.sizes/sum(cluster.sizes) # Normalise priors

# Gaussian standard deviations are estimated on each cluster
  for (j in 1:m){
    params.init[3,j] <- mad(x[which(clustering$cluster==j)])
  }
  if(any(params.init[3,]<0.1))
    params.init[3,which(params.init[3,]<0.1)] <- 0.101
  params.init
}


ummgbwd <- function(params,m){
  cfunction <- .C("ummgbwd",as.double(t(params)), as.integer(m), theta=double(3*m-1),PACKAGE="condmixt")
  cfunction[["theta"]]
}


pgaussmixt <- function(params,x,trunc=TRUE){
  # params is 3 x m matrix
  m <- dim(params)[2]
  p <- rep(0,length(x))
  for (j in 1:m){
    p <- p+params[1,j]*pnorm(x,params[2,j],params[3,j])
  }
  if (trunc){
    F0 <- pgaussmixt(params,0,trunc=FALSE)
    p <- p/(1-F0)
    if (any(x<=0))
      p[which(x<=0)] <- 0
  }
  p
}

dgaussmixt <- function(params,x,log=FALSE,trunc=TRUE){
  # params is 3 x m matrix
  m <- dim(params)[2]
  p <- rep(0,length(x))
  for (j in 1:m){
    p <- p+params[1,j]*dnorm(x,params[2,j],params[3,j])
  }

  if (trunc){
    F0 <- pgaussmixt(params,0,trunc=FALSE)
    p <- p/(1-F0)
    if (any(x<=0))
      p[which(x<=0)] <- 0
  }
  if (log)
    p <- log(p)
  p
}

dlognormixt <- function(params,x,log=FALSE){
  # params is 3 x m matrix
  m <- dim(params)[2]
  p <- rep(0,length(x))
  for (j in 1:m){
    p <- p+params[1,j]*dlnorm(x,params[2,j],params[3,j])
  }
  if (log)
    p <- log(p)
  p
}




condhparetomixt.init <- function(d,h,m,y=NULL){
# D: dimension of input to neural network
# H: number of hidden units
# M: number of components
# Y: targets

  set.seed(as.numeric(Sys.time()))
  nout <- 4*m-1  # number of neural network outputs

  theta <- (runif(nout*(d+1))*1.8 - 0.9)/sqrt(1+d); # linear weights

  if (h >0) {
    for (k in 1:h){
      newneuron <- c((runif(1+d)*1.8 - 0.9)/sqrt(1+d),(runif(nout)*1.8 - 0.9)/sqrt(h+1+d))
      theta <- c(theta, newneuron)
    }
  }
  if(!is.null(y)){ #initialize bias according to unconditional distribution
    if(m>1){
      params <- hparetomixt.init(m,y)
      bias <- ummhbwd(as.vector(t(params)),m)
      ind <-seq(1,nout*(d+1),d+1)
      theta[ind[is.finite(bias)]] <- bias[is.finite(bias)]
    }
    else{
      thetainit <- hpareto.mme(y)
      thetainit[1] <- softplusinv(thetainit[1])
      thetainit[3] <- softplusinv(thetainit[3])
      theta[seq(1,nout*(d+1),d+1)] <- thetainit
    }
  }
  theta
}


condhparetomixt.dirac.init <- function(d,h,m,y=NULL){
# D: dimension of input to neural network
# H: number of hidden units
# M: number of components
# Y: targets

  set.seed(as.numeric(Sys.time()))
  nout <- 4*m  # number of neural network outputs

  theta <- (runif(nout*(d+1))*1.8 - 0.9)/sqrt(1+d); # linear weights

  if (h >0) {
    for (k in 1:h){
      newneuron <- c((runif(1+d)*1.8 - 0.9)/sqrt(1+d),(runif(nout)*1.8 - 0.9)/sqrt(h+1+d))
      theta <- c(theta, newneuron)
    }
  }
  if(!is.null(y)){ #initialize bias according to unconditional distribution
    ypos <- y[y>0]
    theta[1] <- -log(length(y)/length(ypos)-1)#empirical proportion of positive outcomes
    if(m>1){
      params <- hparetomixt.init(m,ypos)
      theta[seq(d+2,nout*(d+1),d+1)] <- ummhbwd(as.vector(t(params)),m)
    }
    else{
      thetainit <- hpareto.mme(ypos)
      thetainit[1] <- softplusinv(thetainit[1])
      thetainit[3] <- softplusinv(thetainit[3])
      theta[seq(d+2,nout*(d+1),d+1)] <- thetainit
    }
  }
  theta
}


condhparetomixt.nll <- function(theta,h,m,x,y){
 # X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  if (length(y) != n)
    stop("Y must have the same number of elements as X")
  nout <- 4*m-1
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmhnllR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x), as.double(y), as.integer(n),nll=double(1),nllgrad=double(length(theta)),PACKAGE="condmixt")
  obj <- results[["nll"]]
  attr(obj,"gradient") <-  results[["nllgrad"]]
  obj
}

condhparetomixt.fit <- function(theta,h,m,x,y,...){
  opt <- nlm(condhparetomixt.nll,theta,h,m,x,y,...)
  opt$estimate
}


condhparetomixt.train <- function(h,m,x,y,nstart=1,...){
# X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  if (length(y) != n)
    stop("Y must have the same number of elements as X")

  nllbest <- Inf
  thetabest <- c()

  for (i in 1:nstart){
    theta0 <- condhparetomixt.init(d,h,m,y) # initialization
    cat("Initial loglike:",condhparetomixt.nll(theta0,h,m,x,y),"\t")
    opt <- try(nlm(condhparetomixt.nll,theta0,h,m,x,y,...))
    if (length(opt) > 1){ # means nlm suceeded
      cat("nstart ", i, ":  ",opt$min, "code ", opt$code, "\n")
      if (opt$min<nllbest){
        nllbest <- opt$min
        thetabest <- opt$est
      }
    } else {
      cat("nstart ", i, "failed\n")
    }
  }
  thetabest
}


condhparetomixt.foldtrain <- function(xtrain,ytrain,xtest,ytest,hp,nstart=1,...){
  # Provide the data on one fold only to parallelize n-fold cv
  # xtrain: d by ntrain matrix of covariates to train the various model
  # ytrain: vector of observations of the dependant variable
  # hp: matrix of hyper-parameters: hp[j,1] = h, hp[j,2]= m
  if(!is.matrix(xtrain))
      stop("XTRAIN is a d x n matrix, even if d = 1")
  d <- dim(xtrain)[1]
  n <- dim(xtrain)[2]

  if (length(ytrain) != n)
    stop("YTRAIN must have the same number of elements as XTRAIN")

  if (dim(xtrain)[1] != dim(xtest)[1])
    stop("Training and test data must have the same dimensions")

  nhp <- nrow(hp)  # number of hyper-parameters
  nlltest <- rep(0,nhp)

  cat("Training fold with ", length(ytrain), " obs \n")

  for (j in 1:nhp){
    cat(j,": h = ",hp[j,1]," m = ",hp[j,2],"\n ")
    thetaopt <- condhparetomixt.train(hp[j,1],hp[j,2],xtrain,ytrain,nstart,...)
    if (is.null(thetaopt)){ # means training failed
      nlltest[j] <- NaN # this set of parameters is not good
    } else {
      nlltest[j] <- condhparetomixt.nll(thetaopt,hp[j,1],hp[j,2],xtest,ytest)
    }
    cat("  Test error is :",nlltest[j],"\n")
  }  # HP loop
  nlltest
}


condhparetomixt.fwd <- function(theta,h,m,x){
 # X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  nout <- 4*m-1
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")
  if (h>0)
    results <- .C("cmmhfwdR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x), as.integer(n),params.mixt=double(4*m*n),a=double((4*m-1)*n),z=double(h*n),PACKAGE="condmixt")
  else
    results <- .C("cmmhfwdR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x), as.integer(n),params.mixt=double(4*m*n),a=double((4*m-1)*n),z=NULL,PACKAGE="condmixt")
  array(results[["params.mixt"]],c(m,4,n))
}


condhparetomixt.nll.tailpen <- function(theta,h,m,x,y,lambda=0,w=0.2,beta=50,mu=0.2,sigma=0.2){
# X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  if (length(y) != n)
    stop("Y must have the same number of elements as X")
  nout <- 4*m-1
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmhnll_bimodal_tailpenR",as.double(theta),as.integer(d),as.integer(h),
                as.integer(m),as.double(x), as.double(y), as.integer(n),as.double(lambda),
                as.double(w),as.double(beta),as.double(mu),as.double(sigma),nll=double(1),
                nllgrad=double(length(theta)),PACKAGE="condmixt")
  obj <- results[["nll"]]
  attr(obj,"gradient") <-  results[["nllgrad"]]
  obj
}

condhparetomixt.fit.tailpen <- function(theta,h,m,x,y,lambda=0,w=0.2,beta=50,mu=0.2,sigma=0.2,...){
  opt <- nlm(condhparetomixt.nll.tailpen,theta,h,m,x,y,lambda,w,beta,mu,sigma,...)
  opt$estimate
}

condhparetomixt.train.tailpen <- function(h,m,x,y,lambda=0,w=0.2,beta=50,mu=0.2,sigma=0.2,nstart=1,...){
# X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  if (length(y) != n)
    stop("Y must have the same number of elements as X")

  nllbest <- Inf
  thetabest <- c()

  for (i in 1:nstart){
    theta0 <- condhparetomixt.init(d,h,m,y) # initialization
    cat("Initial loglike:",condhparetomixt.nll(theta0,h,m,x,y),"\t")
    opt <- try(nlm(condhparetomixt.nll.tailpen,theta0,h,m,x,y,lambda,w,beta,mu,sigma,...))
    if (length(opt) > 1){ # means nlm suceeded
      cat("nstart ", i, ":  ",opt$min, "code ", opt$code, "\n")
      if (opt$min<nllbest){
        nllbest <- opt$min
        thetabest <- opt$est
      }
    } else {
      cat("nstart ", i, "failed\n")
    }
  }
  thetabest
}


condhparetomixt.cvtrain.tailpen <- function(x,y,hp,nfold=5,nstart=1,...){
  # x: d by n matrix of covariates
  # y: vector of observations of the dependant variable
  # hp: matrix of hyper-parameters: hp[j,1] = h, hp[j,2]= m, hp[j,3] = lambda, hp[j,4] = w,
  # hp[j,5] = beta, hp[j,6] = mu and hp[j,7] = sigma

  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  if (length(y) != n)
    stop("Y must have the same number of elements as X")

  nk <- floor(n/nfold)  # number of observations per fold
  sh <- sample(1:n,n)
  x <- matrix(x[,sh],nrow=d,ncol=n)  # shuffle data
  y <- y[sh]

  nhp <- nrow(hp)  # number of hyper-parameters

  nlltest <- rep(0,nhp)
  for (k in 1:nfold){ # train and test on each fold
    if (k < nfold)
      indk <- (1+(k-1)*nk):(k*nk)
    else
      indk <-  (1+(k-1)*nk):max(k*nk,n)
    cat("Training fold ", k, " with ", n-length(indk), " obs \n")

    for (j in 1:nhp){
      cat(j,": h = ",hp[j,1]," m = ",hp[j,2]," lambda = ", hp[j,3], " w = ", hp[j,4],
          " beta = ", hp[j,5]," mu = ", hp[j,6]," sigma = ",hp[j,7],"\n")
      if (d==1){
        thetaopt <- condhparetomixt.train.tailpen(hp[j,1],hp[j,2],t(x[,-indk]),y[-indk],
                                                  hp[j,3],hp[j,4],hp[j,5],hp[j,6],hp[j,7],nstart,...)
      } else {
        thetaopt <- condhparetomixt.train.tailpen(hp[j,1],hp[j,2],x[,-indk],y[-indk],hp[j,3],hp[j,4],
                                                  hp[j,5],hp[j,6],hp[j,7],nstart,...)
      }

      if (is.null(thetaopt)){ # means training failed somehow
        nlltest[j] <- NaN # this set of parameters is not good
      } else {
        if (d>1){
          nlltest[j] <- nlltest[j]+condhparetomixt.nll(thetaopt,hp[j,1],hp[j,2],x[,indk],y[indk])
        } else {
          nlltest[j] <- nlltest[j]+condhparetomixt.nll(thetaopt,hp[j,1],hp[j,2],
                                                       t(x[,indk]),y[indk])
        }
      }
      cat("  Test error is :",nlltest[j],"\n")
    }  # HP loop
    cat("-------------------------------------------------------------- \n")
  }  # k-fold loop

  nlltest
}


condhparetomixt.foldtrain.tailpen <- function(xtrain,ytrain,xtest,ytest,hp,nstart=1,...){
  # Provide the data on one fold only to parallelize n-fold cv
  # xtrain: d by ntrain matrix of covariates to train the various model
  # ytrain: vector of observations of the dependant variable
  # hp: matrix of hyper-parameters: hp[j,1] = h, hp[j,2]= m, hp[j,3] = lambda, hp[j,4] = w,
  # hp[j,5] = beta, hp[j,6] = mu and hp[j,7] = sigma

  if(!is.matrix(xtrain))
      stop("XTRAIN is a d x n matrix, even if d = 1")
  d <- dim(xtrain)[1]
  n <- dim(xtrain)[2]

  if (length(ytrain) != n)
    stop("YTRAIN must have the same number of elements as XTRAIN")

  if (dim(xtrain)[1] != dim(xtest)[1])
    stop("Training and test data must have the same dimensions")

  nhp <- nrow(hp)  # number of hyper-parameters
  nlltest <- rep(0,nhp)

  cat("Training fold with ", length(ytrain), " obs \n")

  for (j in 1:nhp){
      cat(j,": h = ",hp[j,1]," m = ",hp[j,2]," lambda = ", hp[j,3], " w = ", hp[j,4],
          " beta = ", hp[j,5]," mu = ", hp[j,6]," sigma = ",hp[j,7],"\n ")
    thetaopt <- condhparetomixt.train.tailpen(hp[j,1],hp[j,2],xtrain,ytrain,hp[j,3],hp[j,4],
                                              hp[j,5],hp[j,6],hp[j,7],nstart,...)
    if (is.null(thetaopt)){ # means training failed
      nlltest[j] <- NaN # this set of parameters is not good
    } else {
      nlltest[j] <- condhparetomixt.nll(thetaopt,hp[j,1],hp[j,2],xtest,ytest)
    }
    cat("  Test error is :",nlltest[j],"\n")
  }  # HP loop
  nlltest
}


condhparetomixt.dirac.nll.tailpen <- function(theta,h,m,x,y,lambda=0,w=0.2,beta=50,mu=0.2,
                                              sigma=0.2){
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  if (length(y) != n)
    stop("Y must have the same number of elements as X")
  nout <- 4*m
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmhnll_bimodal_tailpen_diracR",as.double(theta),as.integer(d),as.integer(h),
                as.integer(m),as.double(x), as.double(y), as.integer(n),as.double(lambda),
                as.double(w),as.double(beta),as.double(mu),as.double(sigma),nll=double(1),
                nllgrad=double(length(theta)),PACKAGE="condmixt")
  obj <- results[["nll"]]
  attr(obj,"gradient") <-  results[["nllgrad"]]
  obj
}

condhparetomixt.dirac.fit.tailpen <- function(theta,h,m,x,y,lambda=0,w=0.2,beta=50,mu=0.2,sigma=0.2,...){
  opt <- nlm(condhparetomixt.dirac.nll.tailpen,theta,h,m,x,y,lambda,w,beta,mu,sigma,...)
  opt$estimate
}


condhparetomixt.dirac.fwd <- function(theta,h,m,x){
# X has many observations in a d by n matrix
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  nout <- 4*m
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")
  if (h>0)
    results <- .C("cmmhfwd_diracR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x), as.integer(n),params.mixt=double((4*m+1)*n),a=double(nout*n),z=double(h*n),PACKAGE="condmixt")
  else
    results <- .C("cmmhfwd_diracR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x), as.integer(n),params.mixt=double((4*m+1)*n),a=double(nout*n),z=NULL,PACKAGE="condmixt")
  matrix(results[["params.mixt"]],(4*m+1),n)
}

condhparetomixt.dirac.nll <- function(theta,h,m,x,y){
# X has many observations in a d by n matrix
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  if (length(y) != n)
    stop("Y must have the same number of elements as X")
  nout <- 4*m
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmhnll_diracR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x), as.double(y), as.integer(n),nll=double(1),PACKAGE="condmixt")
  results[["nll"]]
}


condhparetomixt.dirac.fit <- function(theta,h,m,x,y,...){
  opt <- nlm(condhparetomixt.dirac.nll,theta,h,m,x,y,...)
  opt$estimate
}


condhparetomixt.dirac.train.tailpen <- function(h,m,x,y,lambda=0,w=0.2,beta=50,mu=0.2,sigma=0.2,
                                                nstart=1,...){
  # X has many observations in a d by n matrix
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  if (length(y) != n)
    stop("Y must have the same number of elements as X")

  nllbest <- Inf
  thetabest <- c()

  for (i in 1:nstart){
    theta0 <- condhparetomixt.dirac.init(d,h,m,y) # initialization
    cat("Initial loglike:",condhparetomixt.dirac.nll(theta0,h,m,x,y),"\t")
    opt <- try(nlm(condhparetomixt.dirac.nll.tailpen,theta0,h,m,x,y,lambda,w,beta,mu,sigma,...))
    if (length(opt) > 1){ # means nlm suceeded
      cat("nstart ", i, ":  ",opt$min, "code ", opt$code, "\n")
      if (opt$min<nllbest){
        nllbest <- opt$min
        thetabest <- opt$est
      }
    } else {
      cat("nstart ", i, "failed\n")
    }
  }
  thetabest
}


condhparetomixt.dirac.foldtrain.tailpen <- function(xtrain,ytrain,xtest,ytest,hp,nstart=1,...){
  # Provide the data on one fold only to parallelize n-fold cv
  # xtrain: d by ntrain matrix of covariates to train the various model
  # ytrain: vector of observations of the dependant variable
  # hp: matrix of hyper-parameters: hp[j,1] = h, hp[j,2]= m, hp[j,3] = lambda, hp[j,4] = w,
  # hp[j,5] = beta, hp[j,6] = mu and hp[j,7] = sigma
  # X has many observations in a d by n matrix
  if(!is.matrix(xtrain))
    stop("XTRAIN is a d x n matrix, even if d = 1")
  d <- dim(xtrain)[1]
  n <- dim(xtrain)[2]
  if (length(ytrain) != n)
    stop("YTRAIN must have the same number of elements as XTRAIN")

  if (dim(xtrain)[1] != dim(xtest)[1])
    stop("Training and test data must have the same dimensions")

  nhp <- nrow(hp)  # number of hyper-parameters
  nlltest <- rep(0,nhp)

  cat("Training fold with ", length(ytrain), " obs \n")

  for (j in 1:nhp){
    cat(j,": h = ",hp[j,1]," m = ",hp[j,2]," lambda = ", hp[j,3], " w = ", hp[j,4],
        " beta = ", hp[j,5]," mu = ", hp[j,6]," sigma = ",hp[j,7],"\n ")

    thetaopt <- condhparetomixt.dirac.train.tailpen(hp[j,1],hp[j,2],xtrain,ytrain,hp[j,3],hp[j,4],
                                                    hp[j,5],hp[j,6],hp[j,7],nstart,...)
    if (is.null(thetaopt)){ # means training failed
      nlltest[j] <- NaN # this set of parameters is not good
    } else {
      nlltest[j] <- condhparetomixt.dirac.nll(thetaopt,hp[j,1],hp[j,2],xtest,ytest)
    }
    cat("  Test error is :",nlltest[j],"\n")
  }  # HP loop
  nlltest
}



condhparetomixt.quant <- function(theta,h,m,x,p,a,b,trunc=TRUE){
 # X has many observations in a d by n matrix
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  nout <- 4*m-1
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  if (trunc)
    results <- .C("cmmhquant_trunc",as.double(theta),as.integer(d),as.integer(h),
                  as.integer(m),as.double(x),as.integer(n),as.double(p),
                  as.integer(length(p)),as.double(a),as.double(b),xq=double(length(p)*n),PACKAGE="condmixt")
  else
    results <- .C("cmmhquant",as.double(theta),as.integer(d),as.integer(h),as.integer(m),
                  as.double(x),as.integer(n),as.double(p),as.integer(length(p)),
                  as.double(a),as.double(b),xq=double(length(p)*n),PACKAGE="condmixt")
  matrix(results[["xq"]],length(p),n)
}


condhparetomixt.dirac.quant <- function(theta,h,m,x,p,a,b){
# X has many observations in a d by n matrix
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  nout <- 4*m
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmhquant_dirac",as.double(theta),as.integer(d),as.integer(h),
                as.integer(m),as.double(x),as.integer(n),as.double(p),
                as.integer(length(p)),as.double(a),as.double(b),xq=double(length(p)*n),PACKAGE="condmixt")
  matrix(results[["xq"]],length(p),n)
}


condhparetomixt.dirac.condquant <- function(theta,h,m,x,p,a,b){
# X has many observations in a d by n matrix
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  nout <- 4*m
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmhcquant_dirac",as.double(theta),as.integer(d),as.integer(h),
                as.integer(m),as.double(x),as.integer(n),as.double(p),
                as.integer(length(p)),as.double(a),as.double(b),xq=double(length(p)*n),PACKAGE="condmixt")
  matrix(results[["xq"]],length(p),n)
}


pcondhparetomixt <- function(params,m,y,trunc=TRUE){
# params is a  m x 4 x n matrix as computed by condhparetmixt.fwd

  n <- dim(params)[3]
  p <- rep(NaN,n)
  for (i in 1:n){
    cfunction <- .C("ummhcdfR",as.double(params[,,i]),as.integer(m),as.double(y[i]),as.integer(1),
                    pdf=double(1),PACKAGE="condmixt")
    p[i] <- cfunction[["pdf"]]
  }
    if (trunc){
      F0 <- pcondhparetomixt(params,m,rep(0,n),trunc=FALSE)
      p <- (p-F0)/(1-F0)
      if (any(y<=0))
        p[which(y<=0)] <- 0
    }
  p
}


pcondhparetomixt.dirac <- function(params,m,y){
# params is a  (4m+1) x n matrix as computed by condhparetomixt.fwd.dirac
  n <- length(y)
  p <- rep(0,n)
  indpos <- which(y>0)
  F0 <- pcondhparetomixt(array(params[2:(4*m+1),indpos],c(m,4,sum(y>0))),m,rep(0,sum(y>0)),trunc=FALSE)
  p[indpos] <- (1-params[1,indpos])+params[1,indpos]*(pcondhparetomixt(array(params[2:(4*m+1),indpos],c(m,4,sum(y>0))),m,y[indpos],trunc=FALSE)-F0)/(1-F0)
  p[which(y<=0)] <- 1-params[1,which(y<=0)]
  p
}

dcondhparetomixt <- function(params,m,y,log=FALSE,trunc=TRUE){
# params is a  m x 4 x n matrix as computed by condhparetmixt.fwd

  n <- dim(params)[3]
  p <- rep(NaN,n)
  for (i in 1:n){
    cfunction <- .C("ummhpdfR",as.double(params[,,i]),as.integer(m),as.double(y[i]),
                    as.integer(1),pdf=double(1),PACKAGE="condmixt")
    p[i] <- cfunction[["pdf"]]
  }
    if (trunc){
      F0 <- pcondhparetomixt(params,m,rep(0,n),trunc=FALSE)
      p <- p/(1-F0)
      if (any(y<=0))
        p[which(y<=0)] <- 0
    }

  if(log)
    return(log(p))
  else
    return(p)

}

condhparetomixt.dirac.negloglike <- function(params,m,y){
# params is a  (4m+1) x n matrix as computed by condhparetomixt.fwd.dirac
  n <- length(y)
  nll <- rep(0,n)
  if (any(y==0))
    nll[y==0] <- -log(1-params[1,which(y==0)])

  if (any(y>0)){
    indpos <- which(y>0)
    nll[indpos] <- -log(params[1,indpos]) - log(dcondhparetomixt(array(params[2:(4*m+1),indpos],c(m,4,sum(y>0))),m,y[indpos],trunc=TRUE))
  }
  nll
}



condgaussmixt.init <- function(d,h,m,y=NULL){
# Conditional Gaussian mixture initialization
# D: dimension of input to neural network
# H: number of hidden units
# M: number of components
# Y: targets

  set.seed(as.numeric(Sys.time()))
  nout <- 3*m-1  # number of neural network outputs

  theta <- (runif(nout*(d+1))*1.8 - 0.9)/sqrt(1+d); # linear weights

  if (h >0) {
    for (k in 1:h){
      newneuron <- c((runif(1+d)*1.8 - 0.9)/sqrt(1+d),(runif(nout)*1.8 - 0.9)/sqrt(h+1+d))
      theta <- c(theta, newneuron)
    }
  }
  if(!is.null(y)){ #initialize bias according to unconditional distribution
    if(m>1){
      params <- gaussmixt.init(m,y)
      theta[seq(1,nout*(d+1),d+1)] <- ummgbwd(as.vector(t(params)),m)
    }
    else{
      thetainit <- c(mean(y),softplusinv(sd(y)))
      theta[seq(1,nout*(d+1),d+1)] <- thetainit
    }
  }
  theta
}


condgaussmixt.nll <- function(theta,h,m,x,y){
 # X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  if (length(y) != n)
    stop("Y must have the same number of elements as X")
  nout <- 3*m-1
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmgnllR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x),
                as.double(y), as.integer(n),nll=double(1),nllgrad=double(length(theta)),PACKAGE="condmixt")
  obj <- results[["nll"]]
  attr(obj,"gradient") <-  results[["nllgrad"]]
  obj
}


condgaussmixt.fit <- function(theta,h,m,x,y,...){
  opt <- nlm(condgaussmixt.nll,theta,h,m,x,y,...)
  opt$estimate
}


condgaussmixt.train <- function(h,m,x,y,nstart=1,...){
# X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  if (length(y) != n)
    stop("Y must have the same number of elements as X")

  nllbest <- Inf
  thetabest <- c()

  for (i in 1:nstart){
    theta0 <- condgaussmixt.init(d,h,m,y) # initialization
    cat("Initial loglike:",condgaussmixt.nll(theta0,h,m,x,y),"\t")
    opt <- try(nlm(condgaussmixt.nll,theta0,h,m,x,y,...))
    if (length(opt) > 1){ # means nlm suceeded
      cat("nstart ", i, ":  ",opt$min, "code ", opt$code, "\n")
      if (opt$min<nllbest){
        nllbest <- opt$min
        thetabest <- opt$est
      }
    } else {
      cat("nstart ", i, "failed\n")
    }
  }
  thetabest
}


condgaussmixt.foldtrain <- function(xtrain,ytrain,xtest,ytest,hp,nstart=1,...){
  # Provide the data on one fold only to parallelize n-fold cv
  # xtrain: d by ntrain matrix of covariates to train the various models
  # ytrain: vector of observations of the dependent variable
  # hp: matrix of hyper-parameters: hp[j,1] = h, hp[j,2]= m
  if(!is.matrix(xtrain))
      stop("XTRAIN is a d x n matrix, even if d = 1")
  d <- dim(xtrain)[1]
  n <- dim(xtrain)[2]

  if (length(ytrain) != n)
    stop("YTRAIN must have the same number of elements as XTRAIN")

  if (dim(xtrain)[1] != dim(xtest)[1])
    stop("Training and test data must have the same dimensions")

  nhp <- nrow(hp)  # number of hyper-parameters
  nlltest <- rep(0,nhp)

  cat("Training fold with ", length(ytrain), " obs \n")

  for (j in 1:nhp){
    cat(j,": h = ",hp[j,1]," m = ",hp[j,2],"\n ")
    thetaopt <- condgaussmixt.train(hp[j,1],hp[j,2],xtrain,ytrain,nstart,...)
    if (is.null(thetaopt)){ # means training failed
      nlltest[j] <- NaN # this set of parameters is not good
    } else {
      nlltest[j] <- condgaussmixt.nll(thetaopt,hp[j,1],hp[j,2],xtest,ytest)
    }
    cat("  Test error is :",nlltest[j],"\n")
  }  # HP loop
  nlltest
}


condgaussmixt.fwd <- function(theta,h,m,x){
# X has many observations in a d by n matrix
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  nout <- 3*m-1
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")
  if (h>0)
    results <- .C("cmmgfwdR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x),
                  as.integer(n),params.mixt=double(3*m*n),a=double(nout*n),z=double(h*n),PACKAGE="condmixt")
  else
    results <- .C("cmmgfwdR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x),
                  as.integer(n),params.mixt=double(3*m*n),a=double(nout*n),z=NULL,PACKAGE="condmixt")
  array(results[["params.mixt"]],c(m,3,n))
}


condgaussmixt.quant <- function(theta,h,m,x,p,a,b,trunc=TRUE){
 # X has many observations in a d by n matrix
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  nout <- 3*m-1
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")
  if (trunc)
    results <- .C("cmmgquant_trunc",as.double(theta),as.integer(d),as.integer(h),
                  as.integer(m),as.double(x),as.integer(n),as.double(p),
                  as.integer(length(p)),as.double(a),as.double(b),xq=double(length(p)*n),PACKAGE="condmixt")
  else
    results <- .C("cmmgquant",as.double(theta),as.integer(d),as.integer(h),as.integer(m),
                  as.double(x),as.integer(n),as.double(p),as.integer(length(p)),
                  as.double(a),as.double(b),xq=double(length(p)*n),PACKAGE="condmixt")
  matrix(results[["xq"]],length(p),n)
}


condgaussmixt.dirac.init <- function(d,h,m,y=NULL){
# D: dimension of input to neural network
# H: number of hidden units
# M: number of components
# Y: targets

  set.seed(as.numeric(Sys.time()))
  nout <- 3*m  # number of neural network outputs

  theta <- (runif(nout*(d+1))*1.8 - 0.9)/sqrt(1+d); # linear weights

  if (h >0) {
    for (k in 1:h){
      newneuron <- c((runif(1+d)*1.8 - 0.9)/sqrt(1+d),(runif(nout)*1.8 - 0.9)/sqrt(h+1+d))
      theta <- c(theta, newneuron)
    }
  }
  if(!is.null(y)){ #initialize bias according to unconditional distribution
    ypos <- y[y>0]
    theta[1] <- -log(length(y)/length(ypos)-1)#empirical proportion of positive outcomes
    if(m>1){
      params <- gaussmixt.init(m,ypos)
      theta[seq(d+2,nout*(d+1),d+1)] <- ummgbwd(as.vector(t(params)),m)
    }
    else{
      if (mad(ypos)>0)
        thetainit <- c(mean(ypos),softplusinv(mad(ypos)))
      else
        thetainit <- c(mean(ypos),-10)
      theta[seq(d+2,nout*(d+1),d+1)] <- thetainit
    }
  }
  theta
}

condgaussmixt.dirac.nll <- function(theta,h,m,x,y){
 # X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  if (length(y) != n)
    stop("Y must have the same number of elements as X")
  nout <- 3*m
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmgnll_diracR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x), as.double(y), as.integer(n),nll=double(1),nllgrad=double(length(theta)),PACKAGE="condmixt")
  obj <- results[["nll"]]
  attr(obj,"gradient") <-  results[["nllgrad"]]
  obj
}


condgaussmixt.dirac.fit <- function(theta,h,m,x,y,...){
  opt <- nlm(condgaussmixt.dirac.nll,theta,h,m,x,y,...)
  opt$estimate
}

condgaussmixt.dirac.train <- function(h,m,x,y,nstart=1,...){
# X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  if (length(y) != n)
    stop("Y must have the same number of elements as X")

  nllbest <- Inf
  thetabest <- c()

  for (i in 1:nstart){
    theta0 <- condgaussmixt.dirac.init(d,h,m,y) # initialization
    cat("Initial loglike:",condgaussmixt.dirac.nll(theta0,h,m,x,y),"\t")
    opt <- try(nlm(condgaussmixt.dirac.nll,theta0,h,m,x,y,...))
    if (length(opt) > 1){ # means nlm suceeded
      cat("nstart ", i, ":  ",opt$min, "code ", opt$code, "\n")
      if (opt$min<nllbest){
        nllbest <- opt$min
        thetabest <- opt$est
      }
    } else {
      cat("nstart ", i, "failed\n")
    }
  }
  thetabest
}

condgaussmixt.dirac.foldtrain <- function(xtrain,ytrain,xtest,ytest,hp,nstart=1,...){
  # Provide the data on one fold only to parallelize n-fold cv
  # xtrain: d by ntrain matrix of covariates to train the various models
  # ytrain: vector of observations of the dependent variable
  # hp: matrix of hyper-parameters: hp[j,1] = h, hp[j,2]= m
  if(!is.matrix(xtrain))
      stop("XTRAIN is a d x n matrix, even if d = 1")
  d <- dim(xtrain)[1]
  n <- dim(xtrain)[2]

  if (length(ytrain) != n)
    stop("YTRAIN must have the same number of elements as XTRAIN")

  if (dim(xtrain)[1] != dim(xtest)[1])
    stop("Training and test data must have the same dimensions")

  nhp <- nrow(hp)  # number of hyper-parameters
  nlltest <- rep(0,nhp)

  cat("Training fold with ", length(ytrain), " obs \n")

  for (j in 1:nhp){
    cat(j,": h = ",hp[j,1]," m = ",hp[j,2],"\n ")
    thetaopt <- condgaussmixt.dirac.train(hp[j,1],hp[j,2],xtrain,ytrain,nstart,...)
    if (is.null(thetaopt)){ # means training failed
      nlltest[j] <- NaN # this set of parameters is not good
    } else {
      nlltest[j] <- condgaussmixt.dirac.nll(thetaopt,hp[j,1],hp[j,2],xtest,ytest)
    }
    cat("  Test error is :",nlltest[j],"\n")
  }  # HP loop
  nlltest
}


condgaussmixt.dirac.fwd <- function(theta,h,m,x){
# X has many observations in a d by n matrix
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  nout <- 3*m
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")
  if (h>0)
    results <- .C("cmmgfwd_diracR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x), as.integer(n),params.mixt=double((3*m+1)*n),a=double(nout*n),z=double(h*n),PACKAGE="condmixt")
  else
    results <- .C("cmmgfwd_diracR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x), as.integer(n),params.mixt=double((3*m+1)*n),a=double(nout*n),z=NULL,PACKAGE="condmixt")
  matrix(results[["params.mixt"]],3*m+1,n)
}


condgaussmixt.dirac.quant <- function(theta,h,m,x,p,a,b){
 # X has many observations in a d by n matrix
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  nout <- 3*m
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmgquant_dirac",as.double(theta),as.integer(d),as.integer(h),
                as.integer(m),as.double(x),as.integer(n),as.double(p),
                as.integer(length(p)),as.double(a),as.double(b),xq=double(length(p)*n),PACKAGE="condmixt")
  matrix(results[["xq"]],length(p),n)
}


condgaussmixt.dirac.condquant <- function(theta,h,m,x,p,a,b){
 # X has many observations in a d by n matrix
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  nout <- 3*m
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmgcquant_dirac",as.double(theta),as.integer(d),as.integer(h),
                as.integer(m),as.double(x),as.integer(n),as.double(p),
                as.integer(length(p)),as.double(a),as.double(b),xq=double(length(p)*n),PACKAGE="condmixt")
  matrix(results[["xq"]],length(p),n)
}


condlognormixt.nll <- function(theta,h,m,x,y){
  condgaussmixt.nll(theta,h,m,x,log(y))+sum(log(y))
}


condlognormixt.init <- function(d,h,m,y=NULL){
    if (is.null(y)){
        theta0 <- condgaussmixt.init(d,h,m)
    } else {
        theta0 <- condgaussmixt.init(d,h,m,log(y))
    }

}

condlognormixt.fit <- function(theta,h,m,x,y,...){
  opt <- nlm(condgaussmixt.nll,theta,h,m,x,log(y),...)
  opt$estimate
}

condlognormixt.train <- function(h,m,x,y,nstart=1,...){
# X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  if (length(y) != n)
    stop("Y must have the same number of elements as X")

  if (any(y<=0))
    stop("Y must contain strictly positive data")

  logy <- log(y)  # the condgaussmixt models log(y)

  nllbest <- Inf
  thetabest <- c()

  for (i in 1:nstart){
    theta0 <- condgaussmixt.init(d,h,m,logy) # initialization
    cat("Initial loglike:",condlognormixt.nll(theta0,h,m,x,y),"\t")
    opt <- try(nlm(condgaussmixt.nll,theta0,h,m,x,logy,...))
    if (length(opt) > 1){ # means nlm suceeded
      cat("nstart ", i, ":  ",opt$min, "code ", opt$code, "\n")
      if (opt$min<nllbest){
        nllbest <- opt$min
        thetabest <- opt$est
      }
    } else {
      cat("nstart ", i, "failed\n")
    }
  }
  thetabest
}


condlognormixt.foldtrain <- function(xtrain,ytrain,xtest,ytest,hp,nstart=1,...){
  # Provide the data on one fold only to parallelize n-fold cv
  # xtrain: d by ntrain matrix of covariates to train the various models
  # ytrain: vector of observations of the dependent variable
  # hp: matrix of hyper-parameters: hp[j,1] = h, hp[j,2]= m
  if(!is.matrix(xtrain))
      stop("XTRAIN is a d x n matrix, even if d = 1")
  d <- dim(xtrain)[1]
  n <- dim(xtrain)[2]

  if (length(ytrain) != n)
    stop("YTRAIN must have the same number of elements as XTRAIN")

  if (dim(xtrain)[1] != dim(xtest)[1])
    stop("Training and test data must have the same dimensions")

  nhp <- nrow(hp)  # number of hyper-parameters
  nlltest <- rep(0,nhp)

  cat("Training fold with ", length(ytrain), " obs \n")

  for (j in 1:nhp){
    cat(j,": h = ",hp[j,1]," m = ",hp[j,2],"\n ")
    thetaopt <- condlognormixt.train(hp[j,1],hp[j,2],xtrain,ytrain,nstart,...)
    if (is.null(thetaopt)){ # means training failed
      nlltest[j] <- NaN # this set of parameters is not good
    } else {
      nlltest[j] <- condlognormixt.nll(thetaopt,hp[j,1],hp[j,2],xtest,ytest)
    }
    cat("  Test error is :",nlltest[j],"\n")
  }  # HP loop
  nlltest
}


condlognormixt.quant <- function(theta,h,m,x,p,a,b){
  z <- condgaussmixt.quant(theta,h,m,x,p,a,b)
  exp(z)
}


condlognormixt.dirac.negloglike <- function(params,m,y){
 # params is a  (3m+1) x n matrix as computed by condgaussmixt.fwd.dirac
  n <- length(y)
  nll <- rep(0,n)
  if (any(y==0))
    nll[which(y==0)] <- -log(1-params[1,which(y==0)])

  if (any(y>0)){
    indpos <- which(y>0)
    nll[indpos] <- -log(params[1,indpos]) - log(dcondgaussmixt(array(params[2:(3*m+1),indpos],c(m,3,sum(y>0))),m,log(y[indpos]),trunc=FALSE)) + sum(log(y[indpos]))
  }
  nll
}



condlognormixt.dirac.nll <- function(theta,h,m,x,y){
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  if (length(y) != n)
    stop("Y must have the same number of elements as X")
  nout <- 3*m
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmlnll_diracR",as.double(theta),as.integer(d),as.integer(h),as.integer(m),as.double(x),
                as.double(y), as.integer(n),nll=double(1),nllgrad=double(length(theta)),PACKAGE="condmixt")
  obj <- results[["nll"]]
  attr(obj,"gradient") <-  results[["nllgrad"]]
  obj
}


condlognormixt.dirac.init <- function(d,h,m,y=NULL){
# D: dimension of input to neural network
# H: number of hidden units
# M: number of components
# Y: targets

  set.seed(as.numeric(Sys.time()))
  nout <- 3*m  # number of neural network outputs

  theta <- (runif(nout*(d+1))*1.8 - 0.9)/sqrt(1+d); # linear weights

  if (h >0) {
    for (k in 1:h){
      newneuron <- c((runif(1+d)*1.8 - 0.9)/sqrt(1+d),(runif(nout)*1.8 - 0.9)/sqrt(h+1+d))
      theta <- c(theta, newneuron)
    }
  }
  if(!is.null(y)){ #initialize bias according to unconditional distribution
    ypos <- y[y>0]
    theta[1] <- -log(length(y)/length(ypos)-1)#empirical proportion of positive outcomes
    if(m>1){
      params <- gaussmixt.init(m,log(ypos))
      theta[seq(d+2,nout*(d+1),d+1)] <- ummgbwd(as.vector(t(params)),m)
    }
    else{
      thetainit <- c(mean(log(ypos)),softplusinv(sd(log(ypos))))
      theta[seq(d+2,nout*(d+1),d+1)] <- thetainit
    }
  }
  theta
}



condlognormixt.dirac.fit <- function(theta,h,m,x,y,...){
  opt <- nlm(condlognormixt.dirac.nll,theta,h,m,x,y,...)
  opt$estimate
}

condlognormixt.dirac.train <- function(h,m,x,y,nstart=1,...){
# X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  if (length(y) != n)
    stop("Y must have the same number of elements as X")

  if (any(y<0))
    stop("Y must contain positive data")

  nllbest <- Inf
  thetabest <- c()

  for (i in 1:nstart){
    theta0 <- condlognormixt.dirac.init(d,h,m,y) # initialization
    cat("Initial negloglike:",condlognormixt.dirac.nll(theta0,h,m,x,y),"\t")
    opt <- try(nlm(condlognormixt.dirac.nll,theta0,h,m,x,y,...))
    if (length(opt) > 1){ # means nlm suceeded
      cat("nstart ", i, ":  ",opt$min, "code ", opt$code, "\n")
      if (opt$min<nllbest){
        nllbest <- opt$min
        thetabest <- opt$est
      }
    } else {
      cat("nstart ", i, "failed\n")
    }
  }
  thetabest
}

condlognormixt.dirac.foldtrain <- function(xtrain,ytrain,xtest,ytest,hp,nstart=1,...){
  # Provide the data on one fold only to parallelize n-fold cv
  # xtrain: d by ntrain matrix of covariates to train the various models
  # ytrain: vector of observations of the dependent variable
  # hp: matrix of hyper-parameters: hp[j,1] = h, hp[j,2]= m
  if(!is.matrix(xtrain))
      stop("XTRAIN is a d x n matrix, even if d = 1")
  d <- dim(xtrain)[1]
  n <- dim(xtrain)[2]

  if (length(ytrain) != n)
    stop("YTRAIN must have the same number of elements as XTRAIN")

  if (dim(xtrain)[1] != dim(xtest)[1])
    stop("Training and test data must have the same dimensions")

  nhp <- nrow(hp)  # number of hyper-parameters
  nlltest <- rep(0,nhp)

  cat("Training fold with ", length(ytrain), " obs \n")

  for (j in 1:nhp){
    cat(j,": h = ",hp[j,1]," m = ",hp[j,2],"\n ")
    thetaopt <- condlognormixt.dirac.train(hp[j,1],hp[j,2],xtrain,ytrain,nstart,...)
    if (is.null(thetaopt)){ # means training failed
      nlltest[j] <- NaN # this set of parameters is not good
    } else {
      nlltest[j] <- condlognormixt.dirac.nll(thetaopt,hp[j,1],hp[j,2],xtest,ytest) + sum(log(ytest[ytest>0]))
    }
    cat("  Test error is :",nlltest[j],"\n")
  }  # HP loop
  nlltest
}


condlognormixt.dirac.quant <- function(theta,h,m,x,p,a,b){
 # X has many observations in a d by n matrix
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  nout <- 3*m
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmlquant_dirac",as.double(theta),as.integer(d),as.integer(h),
                as.integer(m),as.double(x),as.integer(n),as.double(p),
                as.integer(length(p)),as.double(a),as.double(b),xq=double(length(p)*n),PACKAGE="condmixt")
  matrix(results[["xq"]],length(p),n)
}


condlognormixt.dirac.condquant <- function(theta,h,m,x,p,a,b){
  z <- condgaussmixt.dirac.condquant(theta,h,m,x,p,a,b)
  exp(z)
}


kneigh.condquant <- function(x,y,k=10,p=0.9){  # conditional quantile estimated from neighbouring observations
  n <- length(y)
  nq <- length(p)
  yq <- matrix(nrow=nq,ncol=n)
  # find k nearest neighbors for point i
  x2 <- apply(x^2,2,sum)
  x2mat <- matrix(rep(x2,n),n,n)
  distance <- x2mat+t(x2mat)- 2*t(x)%*%x
  xorder <- apply(distance,1,order)
  neigh <-  xorder[2:(k+1),] # k by n matrix
  ysub <- matrix(y[neigh],n,k,byrow=TRUE) # n by k matrix
  yq <- apply(ysub,1,quantile,p)
}



pcondgaussmixt <- function(params,m,y,trunc=TRUE){
# params is a  m x 3 x n matrix as computed by condgaussmixt.fwd

  n <- dim(params)[3]
  p <- rep(NaN,n)
  for (i in 1:n){
    cfunction <- .C("ummgcdfR",as.double(params[,,i]),as.integer(m),as.double(y[i]),as.integer(1),
                    pdf=double(1),PACKAGE="condmixt")
    p[i] <- cfunction[["pdf"]]
  }
    if (trunc){
      F0 <- pcondgaussmixt(params,m,rep(0,n),trunc=FALSE)
      p <- (p-F0)/(1-F0)
      if (any(y<=0))
        p[which(y<=0)] <- 0
    }
  p
}


dcondgaussmixt <- function(params,m,y,log=FALSE,trunc=TRUE){
# params is a  m x 3 x n matrix as computed by condgaussmixt.fwd

  n <- dim(params)[3]
  p <- rep(NaN,n)
  for (i in 1:n){
    cfunction <- .C("ummgpdfR",as.double(params[,,i]),as.integer(m),as.double(y[i]),as.integer(1),
                    pdf=double(1),PACKAGE="condmixt")
    p[i] <- cfunction[["pdf"]]
  }

      if (trunc){
       F0 <- pcondgaussmixt(params,m,rep(0,n),trunc=FALSE)
      p <- p/(1-F0)
      if (any(y<=0))
        p[which(y<=0)] <- 0
    }

  if(log)
    return(log(p))
  else
    return(p)

}

condgaussmixt.dirac.negloglike <- function(params,m,y){
# params is a  (3m+1) x n matrix as computed by condgaussmixt.fwd.dirac
  n <- length(y)
  nll <- rep(0,n)
  if (any(y==0))
    nll[which(y==0)] <- -log(1-params[1,which(y==0)])

  if (any(y>0)){
    indpos <- which(y>0)
    nll[indpos] <- -log(params[1,indpos]) - log(dcondgaussmixt(array(params[2:(3*m+1),indpos],c(m,3,sum(y>0))),m,y[indpos],trunc=TRUE))
  }
  nll
}


# ----------------------------------------------------- ------------------
# ---------------- conditional Bernouilli-Gamma mixture ------------------
# ----------------------------------------------------- ------------------


condbergamixt.init <- function(d,h,y=NULL){
# D: dimension of input to neural network
# H: number of hidden units
# Y: targets

  set.seed(as.numeric(Sys.time()))
  nout <- 3  # number of neural network outputs

  theta <- (runif(nout*(d+1))*1.8 - 0.9)/sqrt(1+d); # linear weights

  if (h >0) {
    for (k in 1:h){
      newneuron <- c((runif(1+d)*1.8 - 0.9)/sqrt(1+d),(runif(nout)*1.8 - 0.9)/sqrt(h+1+d))
      theta <- c(theta, newneuron)
    }
  }
  if(!is.null(y)){ #initialize bias according to unconditional distribution
    ypos <- y[y>0]
    theta[1] <- -log(length(y)/length(ypos)-1)#empirical proportion of positive outcomes
    xbar <- mean(ypos)
    s2 <- var(ypos)
    thetainit <- softplusinv(c(xbar^2/s2, s2/xbar))
    theta[seq(d+2,nout*(d+1),d+1)] <- thetainit
  }
  theta
}


condbergamixt.fwd <- function(theta,h,x){
 # X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  nout <- 3
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  if (h>0)
    results <- .C("cmmbergam_fwdR",as.double(theta),as.integer(d),as.integer(h),
                  as.double(x), as.integer(n),params.mixt=double(nout*n),
                  a=double(nout*n),z=double(h*n),PACKAGE="condmixt")
  else
    results <- .C("cmmbergam_fwdR",as.double(theta),as.integer(d),as.integer(h),
                  as.double(x), as.integer(n),params.mixt=double(nout*n),
                  a=double(nout*n),z=NULL,PACKAGE="condmixt")
  matrix(results[["params.mixt"]],nout,n)
}

condbergamixt.nll <- function(theta,h,x,y){
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]

  if (length(y) != n)
    stop("Y must have the same number of elements as X")
  nout <- 3
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  results <- .C("cmmbergam_nllR",as.double(theta),as.integer(d),as.integer(h),as.double(x), as.double(y),
                as.integer(n), nll=double(1),nllgrad=double(length(theta)),PACKAGE="condmixt")
  obj <- results[["nll"]]
  attr(obj,"gradient") <-  results[["nllgrad"]]
  obj
}


condbergamixt.fit <- function(theta,h,x,y,...){
  opt <- nlm(condbergamixt.nll,theta,h,x,y,...)
  opt$estimate
}

condbergamixt.train <- function(h,x,y,nstart=1,...){
# X is a d x n matrix, even if d = 1
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  if (length(y) != n)
    stop("Y must have the same number of elements as X")

  if (any(y<0))
    stop("Y must contain positive data")

  nllbest <- Inf
  thetabest <- c()

  for (i in 1:nstart){
    theta0 <- condbergamixt.init(d,h,y) # initialization
    cat("Initial negloglike:",condbergamixt.nll(theta0,h,x,y),"\t")
    opt <- try(nlm(condbergamixt.nll,theta0,h,x,y,...))
    if (length(opt) > 1){ # means nlm suceeded
      cat("nstart ", i, ":  ",opt$min, "code ", opt$code, "\n")
      if (opt$min<nllbest){
        nllbest <- opt$min
        thetabest <- opt$est
      }
    } else {
      cat("nstart ", i, "failed\n")
    }
  }
  thetabest
}


condbergamixt.foldtrain <- function(xtrain,ytrain,xtest,ytest,hp,nstart=1,...){
  # Provide the data on one fold only to parallelize n-fold cv
  # xtrain: d by ntrain matrix of covariates to train the various models
  # ytrain: vector of observations of the dependent variable
  # hp: matrix of hyper-parameters: hp[j,1] = h
  if(!is.matrix(xtrain))
      stop("XTRAIN is a d x n matrix, even if d = 1")
  d <- dim(xtrain)[1]
  n <- dim(xtrain)[2]

  if (length(ytrain) != n)
    stop("YTRAIN must have the same number of elements as XTRAIN")

  if (dim(xtrain)[1] != dim(xtest)[1])
    stop("Training and test data must have the same dimensions")

  nhp <- nrow(hp)  # number of hyper-parameters
  nlltest <- rep(0,nhp)

  cat("Training fold with ", length(ytrain), " obs \n")

  for (j in 1:nhp){
    cat(j,": h = ",hp[j,1],"\n ")
    thetaopt <- condbergamixt.train(hp[j,1],xtrain,ytrain,nstart,...)
    if (is.null(thetaopt)){ # means training failed
      nlltest[j] <- NaN # this set of parameters is not good
    } else {
      nlltest[j] <- condbergamixt.nll(thetaopt,hp[j,1],xtest,ytest)
    }
    cat("  Test error is :",nlltest[j],"\n")
  }  # HP loop
  nlltest
}



condbergamixt.quant <- function(theta,h,x,p){
 # X has many observations in a d by n matrix, n is the number of observations
 # params.mixt is a 3 by n matrix of parameters
  if(!is.matrix(x))
    stop("X is a d x n matrix, even if d = 1")
  d <- dim(x)[1]
  n <- dim(x)[2]
  nout <- 3
  if (length(theta) != nout*(d+1)+h*(d+nout+1))
    stop("The number of parameters should be equal to NOUT(D+1)+H(D+NOUT+1)")

  params.mixt <- condbergamixt.fwd(theta,h,x)
  quant <- matrix(0,nrow=length(p),ncol=n)
  for (i in 1:length(p)){
    indi <- which(p[i] > 1-params.mixt[1,])
    quant[i,indi] = qgamma((p[i]+params.mixt[1,indi]-1)/params.mixt[1,indi],
           shape=params.mixt[2,indi],scale=params.mixt[3,indi])
  }
  quant
}


condbergamixt.negloglike <- function(params,y){
# params.mixt is a  3 x n matrix as computed by condbergamixt.fwd
  n <- length(y)
  nll <- rep(0,n)
  if (any(y==0))
    nll[y==0] <- -log(1-params[1,which(y==0)])

  if (any(y>0)){
    indpos <- which(y>0)
    nll[indpos] <- -log(params[1,indpos]) - dgamma(y[indpos],shape=params[2,indpos],
                                                   scale=params[3,indpos],log=TRUE)
  }
  nll
}
