check.psrf <-
function(post1=NULL, post2=NULL, post3=NULL, post4=NULL, post5=NULL)
{
  if(is.list(post1)){
    MCMC.list <- post1
  }
  else{
    tmp <-  list(post1,post2,post3,post4,post5)
    count <- 5
    if(is.null(post5)) count <- count-1
    if(is.null(post4)) count <- count-1
    if(is.null(post3)) count <- count-1  
    if(is.null(post2)) 
      stop("at last two Markov Chains are needed to compute psrf")
    
    draw <- vector("list", count)
    for(i in 1:length(tmp))
    {
      if(!is.null(tmp[[i]])) draw[[i]] <- mcmc(tmp[[i]])
      else break
    }
    MCMC.list <- mcmc.list(draw)
  }

  x <- MCMC.list  
  Niter <- niter(x)
  Nchain <- nchain(x)
  Nvar <- nvar(x)
  xnames <- varnames(x)
  x <- lapply(x, as.matrix)
  S2 <- array(sapply(x, var, simplify = TRUE), 
              dim = c(Nvar, Nvar, Nchain))
  W <- apply(S2, c(1, 2), mean)
  PD <- is.positive.definite(W)    
  
  if(!PD)
  { 
    psrf.s <- gelman.diag(MCMC.list, multivariate=FALSE)$psrf
    psrf.m <- NULL
    print("the covariance matrix of the posterior samples is not")  
    print("positive definite, and the multivarite psrf cannot be")
    print("computed") 
  }
  else{
    gelman.plot(MCMC.list)
    psrf.s <- gelman.diag(MCMC.list)[[1]]
    psrf.m <- gelman.diag(MCMC.list)[[2]]
  }
  par(mfrow=c(1,2),mar=c(2,2,1,1)) 
  boxplot(psrf.s[,1]); mtext("psrf",1,cex=1.2)  
  boxplot(psrf.s[,2]); mtext("upper bound of 95% CI",1,cex=1.2)
 
  print(psrf.s); 
  print(psrf.m)
  return(list(psrf.s=psrf.s, psrf.m=psrf.m,
              psrf.s.summ = apply(psrf.s,2,summary))) 
}
