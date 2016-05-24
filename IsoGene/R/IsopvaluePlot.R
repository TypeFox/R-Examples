IsopvaluePlot <- function(x, y, niter, stat=c("E2","Williams","Marcus","M","ModifM") ){

  # TODO: add check on niter
  Probe.ID <- row.names(y)
  obs <- IsoGene1(x,y)

  stat <- match.arg(stat)
  obs.up <- switch(stat,
      E2 = obs[[1]],
      Williams = obs[[2]],
      Marcus = obs[[3]],
      M = obs[[4]],
      ModifM = obs[[5]])
  
  obs.dn <- switch(stat,
      E2 = obs[[6]],
      Williams = obs[[7]],
      Marcus = obs[[8]],
      M = obs[[9]],
      ModifM = obs[[10]])
  
  exp.up <- exp.dn <- 1:niter

  ### replicate ?
  
  x.niter <- t(sapply(1:niter, function(i) sample(x)))
    
  for (j in 1:niter){
    exps <- IsoGene1(x.niter[j,], y)

    exp.up[j] <- switch(stat,
        E2 = exps[[1]],
        Williams = exps[[2]],
        Marcus = exps[[3]],
        M = exps[[4]],
        ModifM = exps[[5]])
    
    exp.dn[j] <- switch(stat,
        E2 = exps[[6]],
        Williams = exps[[7]],
        Marcus = exps[[8]],
        M = exps[[9]],
        ModifM = exps[[10]])
    
     cat(paste(j, ". "))
   }
   rawp.up <- sum(obs.up < exp.up) / niter
   rawp.dn <- sum(obs.dn > exp.dn) / niter
   
   if (stat == "E2") rawp.dn <- sum(obs.dn < exp.dn) / niter
   
   par(mfrow = c(2,1))
   hist(exp.up, main = "", nclass = 1000, col = 0, probability = TRUE,
        xlim = c(min(exp.up, obs.up), max(exp.up, obs.up)), xlab = paste(stat))
   dx <- density(exp.up, from = min(exp.up), to = max(exp.up))
   abline(v = obs.up, col = 7, lwd = 3)
   lines(dx$x,dx$y, lwd = 3, col = 5)
   title(paste("Gene: ", Probe.ID, ":p-value^{up}=", rawp.up, sep="")) 
  
   hist(exp.dn, main = "", nclass = 1000, col = 0, probability = TRUE,
        xlim = c(min(exp.dn,obs.dn), max(exp.dn,obs.dn)), xlab = paste(stat))
   dx <- density(exp.dn, from = min(exp.dn), to = max(exp.dn))
   abline(v = obs.dn, col = 7, lwd = 3)

   lines(dx$x, dx$y, lwd = 3, col = 5)
   title(paste("Gene: ", Probe.ID, ":p-value^{down}=", rawp.dn, sep = '')) 
}
