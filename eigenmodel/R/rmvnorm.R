"rmvnorm" <-
function(mu,Sig2){

## sample from multivariate normal distribution

  R<-t(chol(Sig2))
  t(R%*%(rnorm(length(mu),0,1)) +mu)
                           }

