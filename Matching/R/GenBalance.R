MATCHpt <- function(q, df, ...)
  {
    #work around for strange windows-R bug found by example(GenMatch) when MatchBalance is hit.
    #don't know how general it is so let's try to work around it
    ret = suppressWarnings(pt(q,df, ...))
    
    if (is.na(ret)) {
      ret   <- pt(q, df, ...)
      if(is.na(ret))
        warning("pt() generated NaN. q:",q," df:",df,"\n")
    }
    
    return(ret)
  } #end of MATCHpt

Mt.test.pvalue  <- function(Tr, Co, weights)
  {
    v1  <- Tr-Co
    estimate  <- sum(v1*weights)/sum(weights)
    var1  <- sum( ((v1-estimate)^2)*weights )/( sum(weights)*sum(weights) )

    if (estimate==0 & var1==0)
      {
        return(1)
      }
    
    statistic  <- estimate/sqrt(var1)
#    p.value    <- (1-pnorm(abs(statistic)))*2
    p.value    <- (1-MATCHpt(abs(statistic), df=sum(weights)-1))*2

    return(p.value)
  } #end of Mt.test.pvalue

Mt.test.unpaired.pvalue  <- function(Tr, Co, weights)
  {
    obs <- sum(weights)
    
    mean.Tr <- sum(Tr*weights)/obs
    mean.Co <- sum(Co*weights)/obs
    estimate <- mean.Tr-mean.Co
    var.Tr  <- sum( ( (Tr - mean.Tr)^2 )*weights)/(obs-1)
    var.Co  <- sum( ( (Co - mean.Co)^2 )*weights)/(obs-1)
    dim <- sqrt(var.Tr/obs + var.Co/obs)

    if (estimate==0 & dim==0)
      {
        return(1)
      }
    
    statistic  <- estimate/dim

    a1 <- var.Tr/obs
    a2 <- var.Co/obs
    dof <- ((a1 + a2)^2)/( (a1^2)/(obs - 1) + (a2^2)/(obs - 1) )    
    p.value    <- (1-MATCHpt(abs(statistic), df=dof))*2    

    return(p.value)
  } #end of Mt.test.unpaired.pvalue

ks.fast <- function(x, y, n.x, n.y, n)
  {
    w <- c(x, y)
    z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
    z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]        
    
    return( max(abs(z)) )
  } #ks.fast

qqstats <- function (x, y, standardize = TRUE, summary.func=NULL)
{
  if (standardize) {
    vals <- sort(unique(c(x, y)))
    x <- ecdf(x)
    sx <- x(vals)
    y <- ecdf(y)
    sy <- y(vals)
  }  else {
    sx <- sort(x)
    sy <- sort(y)
    lenx <- length(sx)
    leny <- length(sy)
    
    if (leny < lenx) 
      sx <- approx(1:lenx, sx, n = leny, method = "constant")$y
    if (leny > lenx) 
      sy <- approx(1:leny, sy, n = lenx, method = "constant")$y
  }
  diffxy <- abs(sx-sy)
  meandiff <- mean(diffxy)
  maxdiff <- max(diffxy)  
  mediandiff <- median(diffxy)
  if(!is.null(summary.func))
    {
      summarydiff <- summary.func(diffxy)
      
      return(list(meandiff=meandiff, mediandiff=mediandiff, maxdiff=maxdiff,
                  summarydiff=summarydiff, summary.func=summary.func))      
    } else {
      return(list(meandiff=meandiff, mediandiff=mediandiff, maxdiff=maxdiff))
    }
}#qqstats

GenBalanceQQ <- function(rr, X, summarystat="mean", summaryfunc="mean")
  {
    index.treated <- rr[,1]
    index.control <- rr[,2]
    
    nvars <- ncol(X)
    qqsummary   <- c(rep(NA,nvars))
    
    for (i in 1:nvars)
      {    
        
        qqfoo <- qqstats(X[,i][index.treated], X[,i][index.control], standardize=TRUE)
        
        if (summarystat=="median")
          {
            qqsummary[i] <- qqfoo$mediandiff
          } else if (summarystat=="max") {
            qqsummary[i] <- qqfoo$maxdiff    
          } else {
            qqsummary[i] <- qqfoo$meandiff
          }
        
      } #end of for loop
    
    
    if (summaryfunc=="median")
      {
        return(median(qqsummary))
      } else if (summaryfunc=="max")  {
        return(sort(qqsummary, decreasing=TRUE))
      } else if (summaryfunc=="sort")  {
        return(sort(qqsummary, decreasing=TRUE))
      } else {
        return(mean(qqsummary))
      }    
  } #end of GenBalanceQQ

GenBalance <-
  function(rr, X, nvars=ncol(X), nboots = 0, ks=TRUE, verbose = FALSE, paired=TRUE)
  {
    index.treated <- rr[,1]
    index.control <- rr[,2]
    weights <- rr[,3]

    tol  <- sqrt(.Machine$double.eps)
    storage.t <- c(rep(9,nvars))
    storage.k <- c(rep(9,nvars))
    fs.ks     <- matrix(nrow=nvars, ncol=1)
    s.ks      <- matrix(nrow=nvars, ncol=1)  
    bbcount   <- matrix(0, nrow=nvars, ncol=1)
    dummy.indx  <- matrix(0, nrow=nvars, ncol=1)
    
    w  <- c(X[,1][index.treated], X[,1][index.control])
    obs <- length(w)
    n.x  <- length(X[,1][index.treated])
    n.y  <- length(X[,1][index.control])
    cutp <- n.x
    w  <- matrix(nrow=obs, ncol=nvars)

    for (i in 1:nvars)
      {
        w[,i] <- c(X[,i][index.treated], X[,i][index.control])

        if(paired)
          {
            t.out <- Mt.test.pvalue(X[,i][index.treated],
                                    X[,i][index.control],
                                    weights = weights)
          } else {
            t.out <- Mt.test.unpaired.pvalue(X[,i][index.treated],
                                             X[,i][index.control],
                                             weights = weights)
          }
        
        storage.t[i] <- t.out            
        
                                        #      print(length(unique(X[,i])) < 3)

        
        dummy.indx[i]  <- length(unique(X[,i])) < 3

        if (!dummy.indx[i] & ks & nboots > 9)
          {
            fs.ks[i]  <- ks.fast(X[,i][index.treated], X[,i][index.control],
                                 n.x=n.x, n.y=n.y, n=obs)
          } else if(!dummy.indx[i] & ks)
            {
                                        #storage.k[i] <- ks.test(X[,i][index.treated], X[,i][index.control])$p.value
              
              storage.k[i] <- Mks.test(X[,i][index.treated], X[,i][index.control])$p.value

            }
      }#end of i loop

    
    if (ks & nboots > 9)
      {
        for (b in 1:nboots)
          {
            sindx  <- sample(1:obs, obs, replace = TRUE)
            
            for (i in 1:nvars)
              {
                
                if (dummy.indx[i])
                  next;
                
                X1tmp <- w[sindx[1:cutp],i ]
                X2tmp <- w[sindx[(cutp + 1):obs], i]
                s.ks[i] <- ks.fast(X1tmp, X2tmp, n.x=n.x, n.y=n.y, n=obs)
                if (s.ks[i] >= (fs.ks[i] - tol) )
                  bbcount[i]  <-  bbcount[i] + 1
              }#end of i loop
          } #end of b loop
        
        for (i in 1:nvars)
          {
            
            if (dummy.indx[i])
              {
                storage.k[i]  <- 9
                next;
              }
            
            storage.k[i]  <- bbcount[i]/nboots
            
          }
        storage.k[storage.k==9]=storage.t[storage.k==9]
        output <- c(storage.t, storage.k)
      } else if(ks){
        storage.k[storage.k==9]=storage.t[storage.k==9]                
        output <- c(storage.t, storage.k)
      } else {
        output <- storage.t
      }

    if(sum(is.na(output)) > 0) {
      output[is.na(output)] = 2
      warning("output has NaNs")
    }
    
    if (verbose == TRUE)
      {
        cat("\n")
        for (i in 1:nvars)
          {
            cat("\n", i, " t-test p-val  =", storage.t[i], "\n" )
            if(ks)
              cat(" ", i, "  ks-test p-val = ", storage.k[i], " \n",sep="")
          }
        cat("\nsorted return vector:\n", sort(output), "\n")
        cat("number of return values:", length(output), "\n")
      }

    return(output)
  } #end of GenBalance


#
# writing fast KS test
#

KSbootBalanceSummary <- function(index.treated, index.control, X, 
                                 nboots = 1000)
{
  X <- as.matrix(X)
  nvars <- ncol(X)

  tol  <- sqrt(.Machine$double.eps)
  storage.k <- c(rep(NA,nvars))
  storage.k.naive <- c(rep(NA,nvars))  
  fs.ks     <- matrix(nrow=nvars, ncol=1)
  s.ks      <- matrix(nrow=nvars, ncol=1)  
  bbcount   <- matrix(0, nrow=nvars, ncol=1)
  dummy.indx  <- matrix(0, nrow=nvars, ncol=1)
  
  w  <- c(X[,1][index.treated], X[,1][index.control])
  obs <- length(w)
  n.x  <- length(X[,1][index.treated])
  n.y  <- length(X[,1][index.control])
  cutp <- n.x
  w  <- matrix(nrow=obs, ncol=nvars)

  for (i in 1:nvars)
    {
      w[,i] <- c(X[,i][index.treated], X[,i][index.control])
      
      dummy.indx[i]  <- length(unique(X[,i])) < 3

      if (!dummy.indx[i])
        {
          foo <- Mks.test(X[,i][index.treated], X[,i][index.control])
          fs.ks[i] <- foo$statistic
          storage.k.naive[i] <- foo$p.value
        } 
    }#end of i loop

  for (b in 1:nboots)
    {
      sindx  <- sample(1:obs, obs, replace = TRUE)
      
      for (i in 1:nvars)
        {
          
          if (dummy.indx[i])
            next;
          
          X1tmp <- w[sindx[1:cutp],i ]
          X2tmp <- w[sindx[(cutp + 1):obs], i]
          s.ks[i] <- ks.fast(X1tmp, X2tmp, n.x=n.x, n.y=n.y, n=obs)
          if (s.ks[i] >= (fs.ks[i] - tol) )
            bbcount[i]  <-  bbcount[i] + 1
        }#end of i loop
    } #end of b loop

  for (i in 1:nvars)
    {
      if (dummy.indx[i])
        {
          storage.k[i]  <- NA
          next;
        }
      
      storage.k[i]  <- bbcount[i]/nboots
    }

  ret = list(ks.boot.pval=storage.k, ks.naive.pval=storage.k.naive, ks.stat=fs.ks)

  return(ret)
} #end of KSbootBalanceSummary
