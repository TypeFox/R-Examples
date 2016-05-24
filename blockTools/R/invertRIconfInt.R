invertRIconfInt <- function(dat, outcome.var, tr.var, tau.abs.min = -1, tau.abs.max = 1, tau.length = 10, n.sb.p = 100, id.vars, id.vals, exact.vars = NULL, exact.vals = NULL, exact.restr = NULL, exact.alg = "single", covar.vars = NULL, covar.vals = NULL, covar.restr = NULL, covars.ord = NULL, n.tr = 2, tr.names = NULL, assg.prob = NULL, seed = NULL, seed.dist, assg.prob.stat = NULL, trim = NULL, assg.prob.method = NULL, assg.prob.kfac = NULL, distance = "mahalanobis", file.name = "sbout.RData", query = FALSE, verbose = TRUE){
  
  taus <- seq(tau.abs.min, tau.abs.max, length = tau.length)
  ci95 <- ci90 <- ci80 <- array(NA)

  for(i in 1:length(taus)){
    
    tau <- taus[i]
    y <- dat[, outcome.var] 
    t <- dat[, tr.var]	
    
    ## given dat X/Y/Tr, calc W_{\tau_0}^D(T) [see Ho and Imai 2006: 892-893].
    wdt <- sum(t*y)/sum(t) - sum((1-t)*y)/sum(1-t)
    
    wdT.dist <- array(NA)
    
    for(j in 1:n.sb.p){
      
      ## given dat X, do an SB.
      sb1 <- seqblock(id.vars = id.vars, id.vals = dat[1, id.vals], exact.vars = exact.vars, exact.vals = dat[1, exact.vals], exact.restr = exact.restr, exact.alg = exact.alg, covar.vars = covar.vars, covar.vals = dat[1, covar.vals], covar.restr = covar.restr, covars.ord = covars.ord, n.tr = n.tr, tr.names = tr.names, assg.prob = assg.prob, seed = seed, seed.dist = seed.dist, assg.prob.stat = assg.prob.stat, trim = trim, assg.prob.method = assg.prob.method, assg.prob.kfac = assg.prob.kfac, file.name = file.name, verbose = FALSE)
  
      for(n.idx in 2:nrow(dat)){
        sb2k <- seqblock(object = file.name, id.vals = dat[n.idx, id.vals], exact.vals = dat[n.idx, exact.vals], exact.restr = exact.restr, exact.alg = exact.alg, covar.vals = dat[n.idx, covar.vals], covar.restr = covar.restr, covars.ord = covars.ord, n.tr = n.tr, tr.names = tr.names, assg.prob = assg.prob, seed = seed, seed.dist = seed.dist, assg.prob.stat = assg.prob.stat, trim = trim, assg.prob.method = assg.prob.method, assg.prob.kfac = assg.prob.kfac, file.name = file.name, verbose = FALSE)
      }
        
      t.tmp <- sb2k$orig[, "Tr"]
      t.tmp <- as.numeric(as.factor(t.tmp))-1
      
      ## t.tmp is "this rerandomized t".  t is the original randomization.
      wdT.dist[j] <- sum(t.tmp*(y + (1-t)*tau))/sum(t.tmp) - sum((1-t.tmp)*(y - t*tau))/sum(1-t.tmp)
    }
    
    if(.025 <= mean(wdT.dist >= wdt) & mean(wdT.dist >= wdt) <= .975){
      ci95 <- append(ci95, tau)
    }
    if(.05 <= mean(wdT.dist >= wdt) & mean(wdT.dist >= wdt) <= .95){
      ci90 <- append(ci90, tau)
    }
    if(.1 <= mean(wdT.dist >= wdt) & mean(wdT.dist >= wdt) <= .9){
      ci80 <- append(ci80, tau)
    }
    
  }
  
  ci95 <- ci95[2:length(ci95)]
  ci90 <- ci90[2:length(ci90)]
  ci80 <- ci80[2:length(ci80)]
  
  return(list(ci95=ci95, ci90=ci90, ci80=ci80))
  
}