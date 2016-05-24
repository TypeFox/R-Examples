"rwish" <-
function(S0,nu){       # sample from a Wishart distribution
sS0<-chol(S0)
Z<-matrix(rnorm(nu*dim(S0)[1]),nu,dim(S0)[1])%*%sS0
t(Z)%*%Z                     }

