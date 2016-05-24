
#library(rootSolve)      #multiroot
#library(numDeriv)       #jacobian

#######################
#
# function which will return \hat{\eta} and a.var(\hat{\eta})
# defaults are linear mu, constant V, Gaussian (N+1)^{st} gap time
# user also needs to provide initial guess for \eta
#



condGEE <- function(data, start, mu.fn=MU, mu.d=MU.d, var.fn=V, k1=K1.norm, k2=K2.norm, 
  robust=TRUE, asymp.var=TRUE, maxiter=100, 
  rtol=1e-6, atol=1e-8, ctol=1e-8, 
  useFortran=TRUE)
{
 p <- length(start)
 N <- length(data[,1])
 n.covs <- length(data[1,])-3

 uniq <- unique(data[,1])
 n <- length(uniq)
 temp <- c(match(uniq,data[,1]),length(data[,1])+1)
 c <- (temp[2:(n+1)]-temp[1:n]) 

 # default mean function is linear
 MU <- function(theta, covs)
   return(theta %*% t(cbind(rep(1,length(covs[,1])),covs)))
 MU.d <- function(theta, covs)
   return(t(cbind(rep(1,length(covs[,1])),covs)))

 # default V function is 1
 V <- function(theta, covs)
   return(rep(1,length(covs[,1])))

 GEE.mean <- function(theta)
 {
  mu <- mu.fn(theta,matrix(data[,4:(4+n.covs-1)],ncol=n.covs))
  var <- var.fn(theta,matrix(data[,4:(4+n.covs-1)],ncol=n.covs))
  d.mu <- mu.d(theta,matrix(data[,4:(4+n.covs-1)],ncol=n.covs))
  std.gaps <- as.numeric((data[,2]-mu)/sqrt(sig2*var))

  cens <- which(data[,3]==0)
  std.gaps[cens] <- k1(std.gaps[cens])

  tot <- as.vector(d.mu %*% diag(sqrt(1/var),N) %*% std.gaps)

  return(tot/n)
 }

 GEE.var <- function(sig2)
 {
  mu <- mu.fn(theta,matrix(data[,4:(4+n.covs-1)],ncol=n.covs))
  var <- var.fn(theta,matrix(data[,4:(4+n.covs-1)],ncol=n.covs)) 
  std.gaps <- as.numeric((data[,2]-mu)/sqrt(sig2*var))
  std.gaps2 <- std.gaps^2

  cens <- which(data[,3]==0)
  std.gaps2[cens] <- k2(std.gaps[cens])

  tot <- sum(std.gaps2)-N

  return(tot/n)
 }

 U <- function(eta)  
 {
  theta <- eta[1:(p-1)]
  sig2 <- eta[p]

  mu <- mu.fn(theta,matrix(data[,4:(4+n.covs-1)],ncol=n.covs))
  var <- var.fn(theta,matrix(data[,4:(4+n.covs-1)],ncol=n.covs))
  d.mu <- mu.d(theta,matrix(data[,4:(4+n.covs-1)],ncol=n.covs)) 
  std.gaps <- as.numeric((data[,2]-mu)/sqrt(sig2*var))
  std.gaps2 <- std.gaps^2

  cens <- which(data[,3]==0)
  std.gaps2[cens] <- k2(std.gaps[cens])
  std.gaps[cens] <- k1(std.gaps[cens])

  tot.mean <- as.vector(d.mu %*% diag(sqrt(1/var),N) %*% std.gaps)
  tot.var <- sum(std.gaps2)-N

  return(c(tot.mean,tot.var)/n)
 }

 var.S <- function(eta)
 {
  theta <- eta[1:(p-1)]
  sig2 <- eta[p]
  tot <- matrix(rep(0,p^2),p,p)

  for(i in 1:n)
  {
    stop.index <- sum(c[1:i])
    start.index <- stop.index - c[i] + 1
    ind <- (start.index:stop.index)
    len <- length(ind)

    mu <- mu.fn(theta,matrix(data[ind,4:(4+n.covs-1)],ncol=n.covs))
    var <- var.fn(theta,matrix(data[ind,4:(4+n.covs-1)],ncol=n.covs))        
    d.mu <- mu.d(theta,matrix(data[ind,4:(4+n.covs-1)],ncol=n.covs))
   
    std.gaps <- as.numeric((data[ind,2]-mu)/sqrt(sig2*var))
    std.gaps2 <- std.gaps^2

    cens <- which(data[ind,3]==0)
    if(length(cens)>0)
    {
      std.gaps2[cens] <- k2(std.gaps[cens])
      std.gaps[cens] <- k1(std.gaps[cens])
    }

    tot.mean <- as.vector(d.mu %*% diag(sqrt(1/var),len) %*% std.gaps)
    tot.var <- sum(std.gaps2)-len
    temp <- c(tot.mean, tot.var)

    tot <- tot + temp %*% t(temp)
  }
  return(tot/n)
 }

 # doing mean and var one by one is slower but more robust
 if(robust==TRUE)
 {
   conv <- 1
   theta <- start[1:(p-1)]
   sig2 <- start[p]
   old <- rep(0,p-1)
   while(conv>1e-4)
   {
     theta <- multiroot(GEE.mean,theta,maxiter, 
        rtol,atol,ctol,useFortran)$root
     sig2 <- multiroot(GEE.var,sig2,maxiter, 
        rtol,atol,ctol,useFortran)$root
     conv <- sum((theta-old)^2)
     old <- theta
   }
   eta <- c(theta,sig2)
 } else { eta <- multiroot(U,start,maxiter, 
        rtol,atol,ctol,useFortran)$root }

 a.var <- NULL
 if(asymp.var==TRUE)
 {
   bread.inv <- jacobian(U,eta)
   bread <- solve(bread.inv)
   meat <- var.S(eta)

   a.var <- bread %*% meat %*% t(bread) / n
 }

 return(list(eta=eta, a.var=a.var))
}
