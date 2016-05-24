## Wrappers for C implementations of the EM algorithms

# Simple EM
simpleEM = function(x, s.init=100, niter=100) {
  # Get all of the constants
  a = getConstants()
  R1 = a$R1 
  R2 = a$R2
  R3 = a$R3
  R4 = a$R4
  R5 = a$R5
  R = a$R 
  S = a$S
  
  out = .C("EM", x=as.integer(x), n=as.integer(length(x)),
    s.init=as.double(s.init), s=as.double(rep(0,niter)),
    ll=as.double(rep(0,niter)), niter=as.integer(niter),
    R=as.double(c(R1,R2,R3,R4,R5,R)),
    package="darts")
  
  return(list(s.final=out$s[niter], s.init=s.init, s=out$s,
              loglik=out$ll, niter=niter))
}

# EM with arbitrary covariance
generalEM = function(x, Sig.init=c(10^2,10^2,0.1*10*10), niter=100,
  seed=NULL) {
  # Get all of the constants
  a = getConstants()
  R1 = a$R1 
  R2 = a$R2
  R3 = a$R3
  R4 = a$R4 
  R5 = a$R5
  R = a$R 
  S = a$S

  # Areas of different board regions
  ar = rep(0,6)
  ar[1] = pi*R1^2
  ar[2] = pi*(R2^2-R1^2)
  ar[3] = pi*(R3^2-R2^2)/20
  ar[4] = pi*(R5^2-R4^2)/20
  ar[5] = pi*(R^2-R5^2)/20
  ar[6] = pi*(R4^2-R3^2)/20

  if (!is.null(seed)) set.seed(seed) 
  
  out = .C("EMCov", x=as.integer(x), n=as.integer(length(x)),
    Siginit=as.double(Sig.init), Sig1=as.double(rep(0,niter)),
    Sig2=as.double(rep(0,niter)), Sig3=as.double(rep(0,niter)),
    ll=as.double(rep(0,niter)), niter=as.integer(niter),
    computell=as.integer(FALSE), R=as.double(c(R1,R2,R3,R4,R5,R)),
    ar=as.double(ar), ii=as.integer(order(S)),
    package="darts")
  
  return(list(Sig.final=c(out$Sig1[niter],out$Sig2[niter],out$Sig3[niter]),
              Sig.init=Sig.init, Sig=cbind(out$Sig1,out$Sig2,out$Sig3), 
              loglik=out$ll, niter=niter))
}

