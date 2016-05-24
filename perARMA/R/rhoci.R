rhoci <-
function(rhovec,nvec,alpha)
{    nrho=length(rhovec)
     nn=length(nvec)
     if (nrho != nn)  {stop("nrho != nn")}
      zval<-matrix(1,1,nrho)
      zstd<-matrix(0,1,nrho) 
      zupper<-matrix(0,1,nrho)
      zlower<-matrix(0,1,nrho)
      upper<-matrix(1,1,nrho)
      lower<-matrix(1,1,nrho)

       good=which(rhovec!=1)
       zval[good]=0.5*log((1+rhovec[good])/(1-rhovec[good]))
       zstd[good]=sqrt(1/(nvec[good]-3))

       del=qnorm(1-alpha/2)
       zupper[good]=zval[good]+del*zstd[good]
       zlower[good]=zval[good]-del*zstd[good]

      lower[good]=tanh(zlower[good])
      upper[good]=tanh(zupper[good])

    result = list(lower=lower,upper=upper)
    class(result) = "rhoci"
    result
}

