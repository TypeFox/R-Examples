
# randcmatrix=random partition matrix for n observations to cln clusters
randcmatrix <- function (n,cln,p){
  ct <- 0
  while(ct<p+3){
    m <- rep(0, times=n*cln)
    summ <- rep(0, times=cln)
    dim(m) <- c(n,cln)
    for (i in 1:n){
      nummer <- round(0.5+cln*runif(1))
#      print(m[i])
      m[i,nummer] <- 1
      summ[nummer] <- summ[nummer] + 1
    }
    ct <- min(summ)
  }
  m
}


# Regression EM iteration;
# m = random partition matrix, cln = number of clusters, icrit = iteration
# stopping criterion, minsig = minimal error variance
# output: coef=regression coefficients, var=residual variances, eps=cluster 
# proportions, z=posterior probabilities, loglik= loglikelihood, warn= T
# if too small or collinear cluster 
regem <- function (indep, dep, m, cln, icrit=1.e-5, minsig=1.e-6,
                         warnings=FALSE) {
  n <- length(dep)
  p <- ncol(as.matrix(indep))
  loglik <- (-1.e8)
  eps <- rep(0,cln)
  fv <- rep(0,n*cln)
  dim(fv) <- c(n,cln)
  rc <- rep(0,(p+1)*cln)
  dim(rc) <- c(p+1,cln)
  rv <- rep(0,cln)
  stm <- rep(0,n)  
  change <- TRUE
  smallcluster <- FALSE
  while (change & !smallcluster) {
#    plot(indep,dep)
    for(i in 1:cln){
      eps[i] <- sum(m[,i])/n
      if (sum(m[,i]>0.01) < p+2){
        if (warnings) warning("Too small cluster")
        smallcluster <- TRUE
      } # if too small cluster
      else{
        reg <- lm(dep~indep, weights=m[,i])
        fv[,i] <- fitted.values(reg)
        rc[,i] <- coefficients(reg)
#        abline(rc[,i],col=i)
        for (j in 2:(p+1))
          if (is.na(rc[j,i])){
            smallcluster <- TRUE
            if (warnings) warning("Collinear regressors")
          } # if collinearity
        res <- residuals(reg)
        rv[i] <- weighted.mean(res^2,m[,i])
        if (rv[i]<minsig){
          rv[i] <- minsig
          if (warnings) warning("Error variance smaller than minimum.")
        } # if error variance below minimum
      } # else (cluster large enough)
    } # for i
    if(!smallcluster){
        for (i in 1:cln)
          for (j in 1:n){
# cat("i= ",i," j= ",j, "dep[j]= ",dep[j]," mean= ",fv[j,i]," sd= ",sqrt(rv[i]),    "\n")
            m[j,i] <- eps[i]*dnorm(dep[j],mean=fv[j,i], sd=sqrt(rv[i]))
          }
        for (j in 1:n){
          stm[j] <- sum(m[j,])
          for (i in 1:cln)
            m[j,i] <- m[j,i]/stm[j]
        } # for j
        oldlog <- loglik        
        loglik <- sum(log(stm))
        change <- (loglik - oldlog > icrit)
    } # if no collinearity & clusters large enough
  } # while change
  g <- c()
  for (i in 1:n)
    g[i] <- which.max(m[i,])
  out <- list(coef=rc, vars=rv, z=m, g=g, eps=eps, loglik=loglik,
              warn=smallcluster)
  out
} # regem     



# Regression mixture analysis (DeSarbo and Cron), 
# ir=iteration runs, nclust= cluster numbers vector, icrit=iteration stopping 
# criterion, minsig = minimum error variance
regmix <- function (indep, dep,
                    ir=1, nclust=1:7, icrit=1.e-5, minsig=1.e-6,
                    warnings=FALSE){
  n <- length(dep)
  p <- ncol(as.matrix(indep))
  clnopt <- min(nclust)
  czmax <- max(nclust)  
  bic <- loglik <- (-1.e9)
  clbic <- rep((-1.e9), czmax)
  eps <- rep(0, czmax)
  rc <- rep(0,(p+1)*czmax)
  dim(rc) <- c(p+1,czmax)
  rv <- rep(0,czmax)
  z <- rep(0, n*czmax)
  dim(z) <- c(n,czmax)
  for (cln in nclust){
    for (i in 1:ir){
      cat("Iteration ",i," for ",cln," clusters.\n")
      emi <- regem(indep, dep, m=randcmatrix(n,cln,p), cln=cln,
                         icrit=icrit, minsig=minsig, warnings=warnings)
      if (emi$warn)   
        emi <- regem(indep, dep, m=randcmatrix(n,cln,p), cln=cln,
                           icrit=icrit, minsig=minsig, warnings=warnings)
      if (!emi$warn){
        bicval <- 2*emi$loglik - log(n)*((p+3)*cln-1)
        if (bicval > clbic[cln])
          clbic[cln] <- bicval
        if (bicval > bic){
          clnopt <- cln
          bic <- bicval
          loglik <- emi$loglik
          eps[1:cln] <- emi$eps
          rc[,1:cln] <- emi$coef
          rv[1:cln] <- emi$var
          z[,1:cln] <- emi$z
        } # if bicval>bic
      }   # if no warning
    }     # for i
  }       # for cln
  g <- c()
  for (i in 1:n)
    g[i] <- which.max(z[i,1:clnopt])
  out <- list(clnopt=clnopt, loglik=loglik, bic=clbic,
              coef=rc[,1:clnopt], var=rv[1:clnopt], eps=eps[1:clnopt], 
              z=z[,1:clnopt], g=g)
  out
# clnopt: Optimal number of clusters, loglik: Loglikelihood, bic: Vector of
# BIC values, coef: Regression coefficients, var: Error variances: 
# eps: cluster proportions, z:a posteriori probabilities, g:optimal
# classification
}          
        





























