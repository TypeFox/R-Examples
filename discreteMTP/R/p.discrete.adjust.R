p.discrete.adjust.methods <- c("BH","BL","BHmidp","BLmidp","DBH","DBL","none")

p.discrete.adjust = function(p, pCDF, method = p.discrete.adjust.methods, cutoff = 1, n = length(p))
  {
    method <- match.arg(method)
    nm <- names(p)
    p <- as.numeric(p)
    names(p) <- nm
    p0 <- p
    if (all(nna <- !is.na(p))) 
        nna <- TRUE
    p <- p[nna]
    lp <- length(p)
    stopifnot(n >= lp)
    if (n <= 1) 
      return(p0)
    
    if (method %in% c("BHmidp","BLmidp","DBH","DBL"))
    {
      midp <- numeric(lp)
      in.CDF <- numeric(lp)
      for (i in 1L:lp)
      {
        in.CDF[i] <- match(p[i],pCDF[[i]])
        if (is.na(in.CDF[i]))
        {
          in.CDF[i] <- which.min(abs(pCDF[[i]]-p[i]))
          p[i] <- pCDF[[i]][in.CDF[i]]
          warning("Since the p-value ",p[[i]]," is not a value the CDF of the p-value,
                  the p-value is rounded to be ",pCDF[[i]][in.CDF[i]],call.=F)
        }  
        
        if (method %in% c("BHmidp","BLmidp"))
          midp[i] <- sum(pCDF[[i]][(in.CDF[i]-1):in.CDF[i]])/2
      }
    }

    if (method %in% c("DBL","DBH"))
    {
       o <- order(p)
       ro <- order(o)
       lc  <- sum(p[o] <= cutoff)
       stepf <- lapply(pCDF,function(x) stepfun(x,c(0,x)))
       p.mat <- matrix(NA,lp,lc)
      
       for (i in 1L:lp)
         p.mat[i,] <- stepf[[i]](p[o][1L:lc])
    }
      
    p0[nna] <- switch(method, DBH = {

        p.sums <- apply(p.mat,2,sum)
        p.sums <-c(p.sums/(1:lc),rep(1,lp-lc))  
        p1 <- numeric(lp)
        p1[lp:1L] = cummin(p.sums[lp:1L])   
        pmin(1,p1[ro])
        
    }, DBL = {
      
        p.prod <- numeric(lc)  
        for (i in 1L:lc)
        {
          p.prod[i] <- prod(1-p.mat[o[i:lp],i])
        }
        p.prod <- c((lp-(1:lc)+1)/lp*(1-p.prod),rep(1,lp-lc))
        cummax(p.prod)[ro]
       
    }, BLmidp = {
        
        i <- 1L:lp
        o <- order(midp, decreasing = FALSE)
        ro <- order(o)
        p1 <- ((lp-i+1)/lp*(1-(1-midp[o])^(lp-i+1)))
        cummax(p1)[ro]
    
    }, BHmidp = {
        
        i <- lp:1L
        o <- order(midp, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin(lp/i * midp[o]))[ro]

    }, BL = {
        
        i <- 1L:lp
        o <- order(p, decreasing = FALSE)
        ro <- order(o)  
        p1 <- ((lp-i+1)/lp*(1-(1-p[o])^(lp-i+1)))
        cummax(p1)[ro]
    
    }, BH = {
      
        i <- lp:1L
        o <- order(p, decreasing = TRUE)
        ro <- order(o)
        pmin(1, cummin(lp/i * p[o]))[ro]
        
    },  none = p)

    p0
}