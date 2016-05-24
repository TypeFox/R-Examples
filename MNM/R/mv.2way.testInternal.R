### Internal functions for
### mv.2way.test
###


###
###
###


`outer.sign` <- function(x, block, eps, maxiter)
{
  d <- dim(x)[2]
  n <- nlevels(block)
  N <- dim(x)[1]

  #Subtract blockwise spatial medians from the blocks
  mx<-by(x,block,spatial.median,eps=eps,maxiter=maxiter)
  mx<-matrix(unlist(mx),ncol=d,byrow=T)
  z<-x-mx[block,]

  #Calculate spatials signs of the centered observations
  sx<-spatial.sign(z,center=F,shape=F)

  return(sx)
}


###
###
###

`inner.sign` <- function(x, block, eps=1.0e-10, maxiter=10000)
{
  d <- dim(x)[2]
  n <- nlevels(block)
  N <- dim(x)[1]

  G <- diag(d)
  iter<-0
  diff<-abs(eps)+1
  smed<-matrix(0,n,d)

  while((diff>eps)&(iter<maxiter))
    {
      iter<-iter+1
      sqrtG<-mat.sqrt(G)               
      #z <- x%*%t(solve(sqrtG))
      z <- tcrossprod(x, syminv(sqrtG))
        
      #Subtract blockwise spatial means from the blocks
      mx<-by(z,block,spatial.median)
      mx<-matrix(unlist(mx),ncol=d,byrow=T)
      z<-z-mx[block,]   

      sx<-spatial.sign(z,center=F,shape=F)
      Cs<-(1/N)*crossprod(sx,sx)
      G<-d*sqrtG%*%Cs%*%sqrtG
      diff<-sqrt(sum((Cs-(1/d)*diag(d))^2))
    }
  if(iter==maxiter)stop("Inner standardization algorithm did not converge")

  return(list(sx=sx,sqrtG=sqrtG))
}


###
###
###

`outer.rank` <- function(x, block)
{
  n <- nlevels(block)

  rx <- x
  for(i in 1:n)
    {
      rx[block==levels(block)[i], ] <-
        spatial.rank(x[block==levels(block)[i], ], shape=FALSE)
    }
  
  return(rx)
}


###
###
###

`inner.rank` <- function(x, block, eps=1.0e-10, maxiter=1000)
{
  d <- dim(x)[2]
  n <- nlevels(block)
  N <- dim(x)[1]

  G <- diag(d)
  iter<-0
  eps<-abs(eps)
  diff<-eps+1
  rx<-x

  while((diff>eps)&(iter<maxiter))
    {
      iter<-iter+1
      sqrtG<-mat.sqrt(G)
      #xg<-x%*%t(solve(sqrtG))
      xg<-tcrossprod(x, syminv(sqrtG))
      for(i in 1:n)
        {
          rx[block==levels(block)[i],]<-
            spatial.rank(xg[block==levels(block)[i],],shape=FALSE)
        }
      Cr<-(1/N)*crossprod(rx,rx)
      coeff<-mean(apply(rx^2,1,sum))
      G<-(d/coeff)*sqrtG%*%Cr%*%sqrtG
      diff<-sqrt(sum((Cr-(coeff/d)*diag(d))^2))
    }
  if(iter==maxiter)stop("Inner standardization algorithm did not converge")
  return(list(rx=rx,sqrtG=sqrtG))
}


###
###
###

`manova.identity` <- function(x, block, treatment,
                              method=c("approximation","permutation"),
                              nsim = 1000)
{
  METHOD <- "MANOVA test in a randomized complete block design"
  d <- dim(x)[2]
  k <- nlevels(treatment)
  n <- nlevels(block)
  N <- dim(x)[1]
  method <- match.arg(method)

  #Subtract block means from the blocks
  mx<-by(x,block,colMeans)
  mx<-matrix(unlist(mx),ncol=d,byrow=T)
  z<-x-mx[block,]

  #Calculate the value of the test statistic
  C <- (1/N)*crossprod(z)
  #invC <- solve(C)
  invC <- syminv(C)
  rdotj<-by(z,treatment,colSums)
  rcr<-sapply(rdotj, my.quad.from, B.inv = invC, simplify = T)
  w0 <- ((k-1)/(n*k))*sum(rcr)

  #Calculate the p-value
  res1<-
    switch(method,
           "approximation"=
           {
             parameter <- d*(k-1)
             names(parameter) <- "df"
             pasym <- 1-pchisq(w0, parameter)
             statistic<-w0
             names(statistic)<-"Q2"
             list(statistic=statistic, p.value = pasym, parameter = parameter, 
                  method = METHOD)
           },
           "permutation"=
           {
             statistics<-
               replicate(nsim,
                         perm.test(z,
                                   block,
                                   factor(replicate(n,sample(k))),
                                   invC))
             pperm <- mean(statistics > w0)
             parameter <- nsim
             names(parameter) <- "replications"
             statistic<-w0
             names(statistic)<-"Q2"
             list(statistic = statistic, p.value = pperm,
                  parameter = parameter, method = METHOD)
           }
           )

  return(res1)
}


###
###
###

`manova.sign` <- function(x, block, treatment,
                          stand=c("outer","inner"),
                          method=c("approximation","permutation"),
                          nsim=1000,eps=1.0e-10,maxiter=10000)
{
  d <- dim(x)[2]
  k <- nlevels(treatment)
  n <- nlevels(block)
  N <- dim(x)[1]
  method <- match.arg(method)

  #Calculate spatial sign vectors
  switch(stand,
         "inner"=
         {
           METHOD <- "Affine invariant multivariate Friedman test using spatial signs"
           sx<-inner.sign(x=x,block=block,eps=eps,maxiter=maxiter)$sx

         },
         "outer"=
         {
           METHOD <- "Multivariate Friedman test using spatial signs"
           sx<-outer.sign(x=x,block=block,eps=eps,maxiter=maxiter)
         }
         )

  #Calculate the value of the test statistic
  Cs <- (1/N)*crossprod(sx)
  #invCs <- solve(Cs)  
  invCs <- syminv(Cs)  
  rdotj<-by(sx,treatment,colSums)
  rcr<-sapply(rdotj, my.quad.from, B.inv = invCs, simplify = T)
  w0 <- ((k-1)/(n*k))*sum(rcr)

  #Calculate the p-value
  res1<-
    switch(method,
           "approximation"=
           {
             parameter <- d*(k-1)
             names(parameter) <- "df"
             pasym <- 1-pchisq(w0, parameter)
             statistic <- w0
             names(statistic)<-"Q2"
             list(statistic=statistic, p.value = pasym, parameter = parameter, 
                  method = METHOD)
           },
           "permutation"=
           {
             statistics<-
               replicate(nsim,
                         perm.test(sx,
                                   block,
                                   c(replicate(n,sample(k))),
                                   invCs))
             pperm <- mean(statistics > w0)             
             parameter <- nsim
             names(parameter) <- "replications"
             statistic <- w0
             names(statistic) <- "Q2"
             list(statistic = statistic, p.value = pperm, parameter = parameter,
                  method = METHOD)             
           }
           )

  return(res1)
}



###
###
###

`manova.rank` <- function(x, block, treatment,
                          stand=c("outer","inner"),
                          method=c("approximation","permutation"),
                          nsim=1000,eps=1.0e-10,maxiter=1000)
{
  d <- dim(x)[2]
  k <- nlevels(treatment)
  n <- nlevels(block)
  N <- dim(x)[1]
  method <- match.arg(method)

  #Calculate the spatial rank vectors
  switch(stand,
         "inner"=
         {
           METHOD <- "Affine invariant multivariate Friedman test using spatial ranks"
           rx<-inner.rank(x=x,block=block,eps=eps,maxiter=maxiter)$rx
         },
         "outer"=
         {
           METHOD <- "Multivariate Friedman test using spatial ranks"
           rx<-outer.rank(x=x,block=block)
         }
         )
  
  #Calculate the value of the test statistic
  Cr <- (1/N)*crossprod(rx)
  # invCr <- solve(Cr)
  invCr <- syminv(Cr)  
  rdotj<-by(rx,treatment,colSums)
  rcr<-sapply(rdotj, my.quad.from, B.inv = invCr, simplify = T)
  w0 <- ((k-1)/(n*k))*sum(rcr)

  #Calculate the p-value
  res1<-
    switch(method,
           "approximation"=
           {
             parameter <- d*(k-1)
             names(parameter) <- "df"
             pasym <- 1-pchisq(w0, parameter)
             statistic<-w0
             names(statistic)<-"Q2"
             list(statistic=statistic, p.value = pasym, parameter = parameter, 
                  method = METHOD)             
           },
           "permutation"=
           {
             statistics<-
               replicate(nsim,
                         perm.test(rx,
                                   block,
                                   c(replicate(n,sample(k))),
                                   invCr))
             pperm <- mean(statistics > w0)                          
             parameter <- nsim
             names(parameter) <- "replications"
             statistic<-w0
             names(statistic)<-"Q2"
             list(statistic = statistic, p.value = pperm, parameter = parameter,
                  method = METHOD)             
           }
           )

  return(res1)
}



###
###
###

`perm.test`<-function (z, block, treatment, B.inv) 
{
  k<-nlevels(treatment)
  n<-nlevels(block)
    
  rdotj<-by(z,treatment,colSums)
  rcr<-sapply(rdotj, my.quad.from, B.inv = B.inv, simplify = T)
  w <- ((k-1)/(n*k))*sum(rcr)
  w
}
