transroutes <- function(ID, GD, sample.times, inf.times, rec.times=NULL, mut.rate, eq.size, bottle.size=1, 
                        p.level=0.95, geninterval=NULL, summary=TRUE) {
  n.cases <- length(ID)
  maxpostsource <- numeric(n.cases)
  if (nrow(GD)!=n.cases || ncol(GD)!=n.cases) {
    stop("Distance matrix dimensions must match length of ID vector")
  }
  
  if (bottle.size>1) {
    if (is.null(geninterval)) {
      geninterval <- mean(inf.times[-1]-inf.times[-length(inf.times)])
    }
  }
  reject <- matrix(0,n.cases,n.cases)
  likelihood <- matrix(0,n.cases,n.cases)
  posterior <- matrix(0,n.cases,n.cases)
  closestsource <- list()
  for (i in 1:n.cases) { # recipient
    mindist <- Inf
    closestsource[[i]] <- 0
    for (j in 1:n.cases) { # source
      if (inf.times[j] < inf.times[i] && (is.null(rec.times[j]) || rec.times[j] > inf.times[i])) {
        distance <- GD[i,j]
        if (distance<mindist) {
          mindist <- distance
          closestsource[[i]] <- ID[j] # Closest source
        } else if (distance==mindist) {
          closestsource[[i]] <- c(closestsource[[i]], ID[j])
        }
        coaltime <- estcoaltime(bottlesize=bottle.size, popsize=eq.size, 
                                bottletimes=ifelse(bottle.size==1,0,seq(from=inf.times[j], to=0, by=-geninterval)), 
                                obstime=min(inf.times[i], sample.times[j]))
        lh <- expsnps(x=distance, m.rate=mut.rate, c.rate=1/coaltime, 
                      tau=abs(sample.times[i]-sample.times[j]))
        likelihood[i,j] <- lh
        cuml <- sum(expsnps(x=0:distance, m.rate=mut.rate, c.rate=1/coaltime, 
                            tau=abs(sample.times[i]-sample.times[j])))
        
        #if (cuml<(1-p.level) || (1-(cuml-lh) > p.level && distance!=0)) {
        if (cuml-lh > p.level && distance!=0) {
          reject[i,j] <- 1
        }
        
      }
    }
    if (sum(likelihood[i,])>0) {
      posterior[i,] <- likelihood[i,]/sum(likelihood[i,])
      maxpostsource[i] <- ID[which.max(posterior[i,])]
    } else {
      maxpostsource[i] <- 0
    }
    
    if (summary) {
      cat("ID ", ID[i], ":\nClosest source: ", paste(closestsource[[i]], collapse=", "), 
          " (", ifelse(closestsource[[i]][1]==0,"-",GD[i,which(ID==closestsource[[i]][1])]), " SNPs", ")",
          "\nMax. posterior source: ", maxpostsource[i], " (PP=", round(max(posterior[i,]),3), ")",
          "\nRejected sources: ", ifelse(sum(reject[i,]==1)==0,"None", paste(ID[which(reject[i,]==1)], collapse=", ")), 
          " (P<", 1-p.level, ")\n\n", sep="")
    }
  }
  return(invisible(list(maxpostsource=maxpostsource, likelihood=likelihood, posterior=posterior, 
                        closestsource=closestsource, reject=reject)))
}