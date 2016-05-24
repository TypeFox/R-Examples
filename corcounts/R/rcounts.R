`rcounts` <-
function(N, margins, mu, phi=rep(NA,length(margins)), omega=rep(NA,length(margins)), psi=rep(NA,length(margins)), corstr, corpar, conv=0.01) {

  mysampling <- function(Ctarget, N, margins, mu, phi, omega, psi, corstr, corpar, conv) {
    error <- FALSE
    dimension <- length(margins)
    Theta <- t(c2pc(Ctarget))
    Theta.c <- Theta


    X <- matrix(NA,N,dimension)
    W <- matrix(NA,N,dimension)
    for (i in 1:dimension) { W[,i] <- runif(N) }
    V <- array(NA,dim=c(N,dimension,dimension-1))

    X[,1] <- W[,1]
    V[,1,1] <- W[,1]

    for (i in 2:dimension) {
      V[,i,1] <- W[,i]
      for (k in (i-1):1) {
          out <- modified.cvine.alg(u1=V[,i,1], u2=V[,k,k],
                               mu.x=mu[i], phi.x=phi[i], omega.x=omega[i],
                               mu.y=mu[k], phi.y=phi[k], omega.y=omega[k],
                               psi.x=psi[i], psi.y=psi[k],
                               rho.target=Theta[i,k], conv=conv,
                               margin.x=margins[i], margin.y=margins[k])
          if (out$fehler) { error <- TRUE }
          V[,i,1] <- out$u1
          Theta.c[i,k] <- out$rhoc

      }
      X[,i] <- V[,i,1]

      if (i<dimension) {
        for (j in 1:(i-1)) {
          V[,i,j+1] <- hfunc(V[,i,j], V[,j,j], Theta.c[i,j])
        }
      }
    }

    Y <- X
    for (i in 1:dimension) {
      if (margins[i]=="Poi")  { Y[,i] <- qpois(X[,i], mu[i]) }
      if (margins[i]=="GP")   { Y[,i] <- pseudoinv.zigp(X[,i], mu[i], phi[i]) }
      if (margins[i]=="ZIP")  { Y[,i] <- pseudoinv.zigp(X[,i], mu[i], 1, omega[i]) }
      if (margins[i]=="ZIGP") { Y[,i] <- pseudoinv.zigp(X[,i], mu[i], phi[i], omega[i]) }
      if (margins[i]=="NB")   { Y[,i] <- qnbinom(X[,i], mu=mu[i], size=psi[i]) }
    }
    return(Y)
  }


  dimension <- length(margins)
  if (corstr=="ex") {
    if (is.matrix(corpar)) stop("'corpar' must be scalar.")
    Ctarget <- R11.exchangeable(corpar,dimension)
  }
  if (corstr=="AR1") {
    if (is.matrix(corpar)) stop("'corpar' must be scalar.")
    Ctarget <- R11.AR1(corpar,dimension)
  }
  if (corstr=="unstr") {
    if (is.matrix(corpar)==FALSE) stop("'corpar' must be a matrix.")
    if (dim(corpar)[1]!=dimension | dim(corpar)[2]!=dimension) stop("Dimension of 'corpar' must be T x T.")
    Ctarget <- corpar
  }
  out <- rep(NA,dimension-1)
  for (i in 2:dimension) {
    out[i-1] <- (det(Ctarget[1:i,1:i])>0)
  }
  if (prod(out)==0) stop("Correlation matrix is not positive definite.")

  # check marginal parameters
  for (i in 1:dimension) {
    if (margins[i]!="ZIGP"&margins[i]!="GP"&margins[i]!="ZIP"&margins[i]!="Poi"&margins[i]!="NB") {
      stop(paste("Invalid parameter token (",margins[i], ") for margin ",i,".",sep=""))
    }
    if (margins[i]=="ZIGP") {
      if (is.na(mu[i])|is.na(phi[i])|is.na(omega[i])) stop(paste("Invalid parameters for margin ",i,".",sep=""))
    }
    if (margins[i]=="GP") {
      if (is.na(mu[i])|is.na(phi[i])) stop(paste("Invalid parameters for margin ",i,".",sep=""))
    }
    if (margins[i]=="ZIP") {
      if (is.na(mu[i])|is.na(omega[i])) stop(paste("Invalid parameters for margin ",i,".",sep=""))
    }
    if (margins[i]=="NB") {
      if (is.na(psi[i])) stop(paste("Invalid parameters for margin ",i,".",sep=""))
    }
  }

  Y <- try(mysampling(Ctarget, N, margins, mu, phi, omega, psi, corstr, corpar, conv),silent=TRUE)
  while (class(Y) == "try-error") {
    Y <- try(mysampling(Ctarget, N, margins, mu, phi, omega, psi, corstr, corpar, conv),silent=TRUE)
  }

  return(Y)
}

