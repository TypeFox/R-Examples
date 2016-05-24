ecm <-function(dat,cvg,eps=1e-6,max.steps=500,eps.bisec=1e-6,ini.m=-7,ini.d=-7){
  m <- dim(dat)
  n <- m[1]
  m <- m[2]
  mu0 <- rep(ini.m,m)
  delta0 <- rep(ini.d,n)
  pm1 <- 0.9
  p0 <- 0.05
  theta0 <- c(mu0,delta0,pm1,p0)
  eout <- estep(mu0,delta0,pm1,p0,dat,cvg)
  mout <- mstep(mu0,delta0,eout[[1]],eout[[2]],eout[[3]],dat,cvg,eps=eps.bisec)
  theta1 <- unlist(mout)
  s <- 1
  while ((max((theta1-theta0)^2)>eps)&(s<=max.steps)){
    theta0 <- theta1
    eout <- estep(mout[[1]],mout[[2]],mout[[3]],mout[[4]],dat,cvg)
    mout <- mstep(mout[[1]],mout[[2]],eout[[1]],eout[[2]],eout[[3]],dat,cvg,eps=eps.bisec)
    theta1 <- unlist(mout)
    s <- s+1
	print(paste0("ecm:",s))
  }
  par.est <- mout
  eout <- estep(mout[[1]],mout[[2]],mout[[3]],mout[[4]],dat,cvg)
  post.probs <- eout
  eoutall <- cbind(c(eout[[1]]),c(eout[[2]]),c(eout[[3]]))
  geno.est <- apply(eoutall,1,which.max)-1
  return(list(par.est=par.est,post.probs=post.probs,steps=s,geno.est=geno.est))
}
