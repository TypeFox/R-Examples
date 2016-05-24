dirichlet.o <- function(data, ncycles=10, M=1,
                            d=c(.1,.1,0,1000), start=NULL, K=100)
{
  cl <- match.call()
  inv.gam.par <- d[1:2]; norm.par <- d[3:4]
  if (!is.numeric(data) || !is.matrix(data))
    stop("data must be numeric matrix")  
  psi.hat <- data[,1]; se.psi.hat <- data[,2]
  if (!all(se.psi.hat >= 0))
    {
     stop("negative standard error in column 2 of data.\n")
    }
  if (M <= 0)  {
   stop("Dirichlet precision parameter M must be greater than zero.\n")
   }
  if (!is.numeric(d) || length(d) != 4){
   stop("hyperparameter must be numeric of length 4")
   }
  else if (d[1] <= 0 || d[2] <= 0 || d[4] <= 0){
     stop("Gamma shape, Gamma scale, and normal variance hyperparameters
          must all be greater than zero\n")
   }
  names(d) <- c("gam.shape","gam.scale","mu","var.factor") 
  aa <- d[1]; bb <- d[2]; cc <- d[3]; dd <- d[4]
  cycle.number <- 0
  nstudies <- length(psi.hat)    
  if (is.null(start)) {
    start.user <- FALSE
    psi.vec.current <- psi.hat
    psi.current <- mean(psi.vec.current)
    tau.current <- sqrt(var(psi.vec.current))
    start <- c(psi.vec.current,psi.current,tau.current)
  }
  else if (!is.null(start)){
    if(!is.numeric(start) || !(length(start)==nstudies+2)) {
     stop("initial param must be num. vec. length nstudies+2, or omit")
    }
    else {
     start.user <- TRUE
     psi.vec.current <- start[1:nstudies]
     psi.current <- start[nstudies+1]
     tau.current <- start[nstudies+2]
    }
  }  
  vec.to.save  <- c(psi.vec.current, psi.current, tau.current)
  V <- vector(length=K)
  mu <- 0
  output.matrix <- matrix(0, ncycles+1, nstudies+3)
  output.matrix[1,] <- c(vec.to.save, mean(psi.hat))
  for (ic in (1:ncycles))
    {
      for (i in (1:nstudies))
        {
          value <- psi.vec.current[-i]; st.dev <- se.psi.hat[i];
          indices.psi.minus.i <- (1:nstudies)[-i]
          cjs <- dnorm(value, mean=psi.hat[i], sd=se.psi.hat[i])
          term1 <- sum(cjs)
          mu <- (psi.current*st.dev^2 + psi.hat[i]*tau.current^2) /
            (st.dev^2 + tau.current^2)
          sigma.sq <- (st.dev^2*tau.current^2) / (st.dev^2 + tau.current^2)
          sigma <- sqrt(sigma.sq)
          term2 <- ((M*sigma)/(sqrt(2*pi)*tau.current*st.dev))*
            exp(-(psi.hat[i] - psi.current)^2 / (2*(st.dev^2+tau.current^2)))
          A <- term1 + term2
          p0 <- term2/A
          p <- cjs/A
          p.cum <- p0 + cumsum(p)
#          print(c(p0,p.cum))
          runi <- runif(n=1)
          uni.categ <- findInterval(runi,c(0,p0,p.cum)) - 1
#          uni.categ <- cut(runi, c(0,p0,p.cum),labels=FALSE) - 1 #20040109 For R
#          print(uni.categ)
#          uni.categ <- cut(runi, c(0,p0,p.cum)) - 1
          if (uni.categ == 0)
            {psi.vec.current[i] <- rnorm(1, mean=mu, sd=sqrt(sigma.sq))}
#          else if (uni.categ > 0)
          if (uni.categ > 0) #20040109 For R
            {jj <- uni.categ;
             psi.vec.current[i] <- psi.vec.current[indices.psi.minus.i[jj]]}
        }
      psi.vec.unique <- unique(psi.vec.current)
      mstar <- length(psi.vec.unique)
      psi.vec.unique.bar <- mean(psi.vec.unique)
      aa.prime <- aa + (mstar/2)
      ml <- mstar*dd
      bb.prime <- bb + .5 * sum((psi.vec.unique - psi.vec.unique.bar)^2) +
        mstar * (psi.vec.unique.bar - cc)^2 / (2*(1+ml))
      cc.prime <- (cc + mstar * dd * psi.vec.unique.bar)/(ml+1)
      dd.prime <- 1/(mstar + (1/dd))
      tau.current <- 1/sqrt(rgamma(1,aa.prime)/bb.prime)
      psi.current <- rnorm(1, mean=cc.prime,
                           sd=sqrt(dd.prime*tau.current^2))
      cycle.number <- cycle.number + 1 #20110701 deleted from output XXX
      vec.to.save <- c(psi.vec.current, psi.current, tau.current)
      p0 <- M/(M+nstudies)
      p <- rep(1/(M+nstudies),nstudies)
      p.cum <- p0 + cumsum(p)
#      print(c(p0,p.cum))
      runi <- runif(n=K)
      uni.categ <- findInterval(runi,c(0,p0,p.cum)) - 1
#      uni.categ <- cut(runi, c(0,p0,p.cum),labels=FALSE) - 1 #9-1 For R
#      print(uni.categ)
#      uni.categ <- cut(runi,c(0,p0,p.cum)) - 1
      for (ib in (1:K))
        {if(uni.categ[ib]==0)
           V[ib] <- rnorm(1,mean=psi.current,sd=tau.current)
#        else if (uni.categ[ib]>0)
         if (uni.categ[ib] > 0) #20040109 For R
          {jj <- uni.categ[ib];
           V[ib] <- psi.vec.current[jj]}
       }
      B <- rbeta(K,1,M+nstudies)
      prodomB <- cumprod(1-B[-c(K-1,K)])
      factor.right <- c(1,prodomB)
      Pkm1 <- B[-K]*factor.right
      P <- c(Pkm1,1-sum(Pkm1))
      muu <- sum(V*P)
      output.matrix[ic+1,] <- c(vec.to.save, muu)
    }
  z <- list(call = cl, ncycles=ncycles, M=M, prior = d,
            chain = output.matrix,start.user=start.user,start=start)
  class(z) <- c("dir.ord")
  z
}

