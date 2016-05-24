robust.factor<-function(data,q)
{
   fac   <- NULL
   S     <- huber.M(data)
   means <- S$loc
   datac <- sweep(data,2,means,"-")
 
   # get matrix square root of S$cov:
   apu   <- eigen(S$cov)
   L     <- apu$values
   P     <- apu$vectors
   z     <- datac %*% P%*%(diag(L^(-1/2)))%*%t(P)
   r     <- sqrt(diag(z%*%t(z)))

   fac[1]<- mean(alpha.fun(r,2,q=pchisq(3,2))^2)/8
   fac[2]<- mean(gamma.fun(r,2,q=pchisq(3,2))^2)/2
   fac
}


alpha.fun<-function(r,k,q) 
{
  c<-qchisq(q,k)
  sig<-pchisq(c,k+2)+(c/k)*(1-q)

  c<-sqrt(c)
  eta<-NULL
  eta[r<=c]<-(r[r^2<=c^2])^2/(2*sig^2) 
  eta[r>c]<-c^2/(4*sig^2) 
  eta <- mean(eta)

  alpha<-NULL
  alpha[r<=c]<-(r[r<=c])^2/(eta*sig^2)
  alpha[r>c]<-c^2/(eta*sig^2)
  alpha
}

gamma.fun<-function(r,k,q) 
{
  c<-sqrt(qchisq(q,k))
  sig<-pchisq(c,k+2)+(c/k)*(1-q)
  c<-sqrt(c)
  eta<-NULL
  eta[r<=c]<-1 
  eta[r>c]<-c*(k-1)/(r[r>c]*k)
  eta<-mean(eta)

  gamma<-NULL
  gamma[r<=c]<-(r[r<=c])/eta
  gamma[r>c]<-c/eta 
  gamma
}
