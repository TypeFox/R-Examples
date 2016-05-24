testOmegaPT <-
function(p,n)
{
   mu1  <- double(n)
   mu   <- double(n)
   theta0 <-double(1)
   w <- double(n)
   w <- double(n)
  for (i in 1: n) 
    { 
      mu[i] <- abs(1/runif(1))
      theta0 <- - abs(1/runif(1))
      w[i]<-omega(p,mu[i],theta0)
      mu1[i]<-moyennePT(p,w[i],theta0)
    }
print("mu = ")
print(mu)
print("mu1 = ") 
print(mu1)
print("mu-mu1 = ") 
n=mu-mu1 
print(n)   
}

