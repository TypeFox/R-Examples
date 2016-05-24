reweight <- function(ori, mar, raw=NA, wgt=NA, unique=T, bound=c(0, 100),
                     trace=F, tolerance=0.1, penalty=0, ...) {
    
 # wgt is the original weights for the corresponding "ori" line
 # raw is the raw counts in survey for the corresponding "ori" line
 if (is.na(raw[1])) raw <- rep(1,dim(ori)[1])
 if (is.na(wgt[1])) wgt <- raw

 if (! unique) {
 # aggregate over "ori" to make unique combination of the categorical variable
   x <- aggregate(data.frame(wgt=wgt,raw=raw),as.list(ori),sum)
   ori <- sapply(x[,names(ori)],as.integer) #convert factor matrix to numeric matrix
   wgt <- x$wgt
   raw <- x$raw
 }
 
  lvls <- apply(ori,2,max)              #number of levels for each categorial variable
  mar <- as.numeric(as.vector(mar))

  # Verify the "mar" vector is consistent.                                        
  total <- sum(mar[1:lvls[1]])
  segment <- rep(1:length(lvls),lvls)
  if (length(mar) != length(segment)) stop("Length of 'mar' is NOT equal to ",length(segment))
  if (sum((tapply(mar,segment,sum)-total)^2)>0.01) {
    cat("\n",tapply(mar,segment,sum),"\n")
    stop("echo mar vector does not add up correctly across all categories!")
  }

  # Construct Z and Y matrix.
  wgttot <- sum(wgt)
  wgt <- wgt/wgttot       #convert weight to proporions so they add up to 1
  mar <- mar/total                  #convert marginal counts to proportions
  dinit <- weight.init(ori,wgt,mar,lvls)    

  # Solve the linear system with SVD + Tikhonov regularization
  dx <- decompose(dinit$X)
  sresult <- search.best.h(dx, dinit,
                           trace=trace, tolerance=tolerance, penalty=penalty)
 
  finalest <- getestimate.h(dx,dinit$Y, sresult$best.h)
  ratio <- finalest$beta +1
  ratio[ratio<bound[1]] <- bound[1]     # Cut 'ratio' to within the range set by 'bound'
  ratio[ratio>bound[2]] <- bound[2] 
  newwgt <- ratio*wgt

  # How original sample weights fits to the marginal.                                         
  wtmp <- c()
  for (i in 1:length(lvls)) {
    wtmp <- append(wtmp,
                   tapply(c(wgt,rep(0,lvls[i])),c(ori[,i],1:lvls[i]),sum))
  }
    #above usage of lvls ensures every level is included in summary, even if
    #not present in data.
  wfit <- unlist(abs(wtmp - mar)/mar)   #Mean absolute deviation from true marginal cnt
  
  # How new sample weights fits to the marginal.                                         
  wtmp <- c()
  for (i in 1:length(lvls)) {
    wtmp <- append(wtmp,
                   tapply(c(newwgt,rep(0,lvls[i])),c(ori[,i],1:lvls[i]),sum))
    }
  wfitnew <- unlist(abs(wtmp - mar)/mar)#Mean absolute deviation from true marginal cnt
  
  ans <- list(newwgt=newwgt, 
              ratio=ratio,ori=ori,wgt=wgt,raw=raw,mar=mar,wgttot=wgttot,total=total,
              lvls=lvls,keplist=dinit$keplist,sresult=sresult,finalest=finalest,
              wfit=wfit,wfitnew=wfitnew,sinval=dx$d)
  class(ans) <- "reweight"
  invisible(ans)
}

print.reweight <- function(x,...) {
  cat("Usage of print(r): Use the returned object (data frame) from print(r) to cross reference the survey respondents from the factor levels and multiply the corresponding original weights by the 'weight.ratio' to get new weights.\n")
  data.frame(x$ori, weight.ratio=x$ratio) }

summary.reweight <- function(object,...) {
  o <- object		 
  newwgt=o$newwgt
  ori=o$ori
  wgt=o$wgt
  mar=o$mar
  wgttot=o$wgttot
  total=o$total
  lvls=o$lvls
  sresult=o$sresult
  
  cat("\n Total weights: (sample) = ",wgttot,",  (margin) = ",total,"\n")

  cat("Number of levels in each of ",length(lvls),"categorical variables:\n")
  print(lvls)
  
  cat("Number of true marginal counts used = ", sum(o$keplist),
      " out of total = ",length(o$keplist),"\n") 

  cat("Summary of nonzero singular values of the design matrix X: (total=",
      sum(lvls)+1-length(lvls),")\n")
  print(summary(o$sinval,...))
  
  cat("Best regularization parameter = ", sresult$best.h," at iteration = ", sresult$iter,"\n")
  
  cat("GCV = ",sresult$min.gcv," with DF = ", length(wgt)-o$finalest$df,"\n")
  
  cat("Percent of weighting ratios needs to be modified to be no less than zero = ",
      round(o$finalest$pctbeta,4)," out of total = ",length(o$finalest$beta),"\n")

  # How original sample weights fits to the marginal.                                         
  cat("\nSummary of Absolute Relative Deviation of Sample Marginal Away From True Marginal:\n ")
  print(summary(o$wfit,...))
  cat("\n")

  # How new sample weights fits to the marginal.                                         
  cat("Summary of Absolute Relative Deviation of New Marginal Away From True Marginal:\n ")
  print(summary(o$wfitnew,...))
  cat("\n")

  cat("Summary of Weighting Ratios:\n ")
  print(summary(o$ratio,...))
  cat("\n")

  invisible()
}


plot.reweight <- function(x,...) {
  o <- x
  newwgt=o$newwgt
  ori=o$ori
  wgt=o$wgt
  mar=o$mar
  wgttot=o$wgttot
  total=o$total
  lvls=o$lvls
  dinit=o$dinit
  sresult=o$sresult
  
  par(mfrow=c(2,2))

  # plot distribution of weighting ratios
  fullratio <- rep(o$ratio,o$raw)
  plot(density(fullratio),xlim=c(0,max(o$ratio)+1),
       main=paste("Weight Ratio Dist (",length(fullratio)," Points)"),
       xlab="Weight Ratio",...)
  rug(fullratio)

  # plot distribution of GCV versus r
  sord=order(sresult$plotx)
  plot(sresult$plotx[sord],sresult$ploty[sord],type='b',main="Search for r that minimizes GCV",
       xlab="Regularization Parameter r",ylab="GCV",...)

  # plot distribution of relative deviation of estimated vs. true marginals.
  xlimup <- round(max(o$wfit)+1)        #uppder boundary for MRD
  plot(density(o$wfit),main="Fit to Marginals (Original Weights)",
       xlim=c(0,xlimup),xlab="Relative Deviation",...)
  rug(o$wfit)
  plot(density(o$wfitnew),main="Fit to Marginals (New Weights)",
       xlim=c(0,xlimup),xlab="Relative Deviation",...)
  rug(o$wfitnew)

  invisible()
}

       
  
weight.init <- function(ori,wgt,mar,lvls) {
  nvar <- length(lvls)                  #number of variables
  totlvls <- sum(lvls)                  #total number of variable levels
  a <- sapply(1:nvar,function(x){tapply(c(wgt,rep(0,lvls[x])),
                                        c(ori[,x],1:lvls[x]),sum)})
                                        # enumerate full list of levels by
                                        # adding pseudo records with
                                        # zero weights
  act.num <- unlist(a)                  #actual number of panelists whose
  act.num[is.na(act.num)] <- 0          #variable Cj takes k-th level
  #Discard equations where the actual number of panelists are zero,
  #meaning there is no need to fit for that marginal component since there
  #is no such variable level present in the current sample.
  keplist <- (act.num!=0)
  Y <- unlist(c(0,mar - act.num))         # the vector Y
  loc <- cumsum(c(0,lvls[1:(nvar-1)])) #initial location of each factor
                                        #in a vector
  Z <- apply(sweep(ori,2,loc,"+"),1,
                 function(x){b <- rep(0,totlvls);b[x] <- 1;c(1,b)})
  X <- sweep(Z,2,wgt,"*")
  list(X=X[keplist,],Y=Y[keplist],keplist=keplist)
}


decompose <- function(x,truncate=1.0e-5)  {
    svdx <- La.svd(x)
    d <- svdx$d
    nonzero <- (d>truncate)
    d <- d[nonzero]    #only nonzero eigenvalues are kept
    #only singular vectors which corresponds to nonzero eigenvalues are kept    
    u <- svdx$u[,nonzero]
    v <- t(svdx$vt[nonzero,])
    d2 <- d^2    
    list(d=d,u=u,v=v,d2=d2)
  }

getestimate.h <- function (dx,Y,h) {
  d <- dx$d
  u <- dx$u
  v <- dx$v
  d2 <- dx$d2    
  h2 <- h^2
  f <- d/(d2+h2)
  betacp <- sapply(1:length(f),function(i){f[i]*sum(u[,i]*Y)*v[,i]})
  beta <- apply(betacp,1,sum)
  #Pct of beta needs modified to be greater than -1 is
  pctbeta <- sum(1+beta<=0)/length(beta)
  beta <- ifelse(1+beta>0,beta,-0.99)
  list(beta=beta,df=sum(f*d),pctbeta=pctbeta)
}

search.best.h <- function(dx,dinit,trace=F,tolerance=0.1,penalty=0)  {
    f <- function(h)  {
        est <- getestimate.h(dx,dinit$Y,h)
        # get GCV function value
    # the term (1+..)^penalty penalize large number of zero weighting ratios
      (  sum((dinit$Y - dinit$X %*% est$beta) ^ 2) * (1+est$pctbeta)^penalty /
       (length(est$beta)- est$df)^2 )
      }

    #implement Golden Selection Search
    #set smallest h to tenth of the smallest nonzero absolute eigenvalues of Z
    a <- delta <- min(0.1,tolerance * min(dx$d))
    #set largest h to the upper 0.8 quantile of the eigenvalues of Z
    b <- quantile(dx$d,0.8)
    t <- (sqrt(5)-1)/2 #golden ratio
    x1 <- a+(1-t)*(b-a)
    f1 <- f(x1)
    x2 <- a+t*(b-a)
    f2 <- f(x2)
    if (trace) cat("tolerance = ", tolerance,"\n")
    iter <- 0
    plotx <- c(x1,x2) #for plotting GCV vs. r
    ploty <- c(f1,f2)
    while ((b-a)> delta ) {
      iter <- iter+1
      if (trace) cat("Object Function on (",a,",",b,") = ",f1,"\n") 
        if (f1>=f2) {
          a <- x1
          x1 <- x2
          f1 <- f2
          x2 <- a+t*(b-a)
          f2 <- f(x2)
          plotx <- append(plotx, x2)
          ploty <- append(ploty, f2)
        }
        else {
          b <- x2
          x2 <- x1
          f2 <- f1
          x1 <- a+(1-t)*(b-a)
          f1 <- f(x1)
          plotx <- append(plotx, x1)
          ploty <- append(ploty, f1)
        }
      }
    best.h <- (a+b)/2
    list(best.h=best.h, min.gcv=f(best.h),iter=iter,plotx=plotx,ploty=ploty)
  }


