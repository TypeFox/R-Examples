#############################################################
#                                                           #
#	wle.binomial function                               #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: February, 22, 2010                            #
#	Version: 0.2-1                                      #
#                                                           #
#	Copyright (C) 2010 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.binomial <- function(x, size, boot=30, group, num.sol=1, raf="HD", tol=10^(-6), equal=10^(-3), max.iter=500, verbose=FALSE) {

result <- list()

if (raf!="HD" & raf!="NED" & raf!="SCHI2") stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

if (missing(group)) {
group <- 0
}

x <- as.vector(x)
nsize <- length(x)
result <- list()

if (nsize<1) {
stop("Number of observation must be at least equal to 1")
}

if (group<1) {
    group <- max(round(nsize/4),1)
    if (verbose) cat("wle.binomial: dimension of the subsample set to default value: ",group,"\n")
}

maxboot <- sum(log(1:nsize))-(sum(log(1:group))+sum(log(1:(nsize-group))))

if (boot<1 | log(boot) > maxboot) {
    stop("Bootstrap replication not in the range")
}

if (!(num.sol>=1)) {
    if (verbose) cat("wle.binomial: number of solution to report set to 1 \n")
    num.sol <- 1
}

if (max.iter<1) {
    if (verbose) cat("wle.binomial: max number of iteration set to 500 \n")
    max.iter <- 500
}

if (tol<=0) {
    if (verbose) cat("wle.binomial: the accuracy must be positive, using default value: 10^(-6) \n")
    tol <- 10^(-6)
}

if (equal<=tol) {
    if (verbose) cat("wle.binomial: the equal parameter must be greater than tol, using default value: tol+10^(-3) \n")
    equal <- tol+10^(-3)
}

tot.sol <- 0
not.conv <- 0
iboot <- 0

while (tot.sol < num.sol & iboot < boot) {
   iboot <- iboot + 1
   x.boot <- x[round(runif(group,0.501,nsize+0.499))]
   p <- sum(x.boot)/(size*group)

   ff <- rep(0,nsize)
   x.diff <- tol + 1
   iter <- 0
   while (x.diff > tol & iter < max.iter) {
   iter <- iter + 1
   p.old <- p 
       tff <- table(x)/nsize
       nff <- as.numeric(names(tff))
       for (i in 1:nsize) {
           ff[i] <- tff[nff==x[i]] 
       }
       mm <- dbinom(x,size=size,prob=p)
       dd <- ff/mm - 1
       
       ww <- switch(raf,
                 HD =  2*(sqrt(dd + 1) - 1) ,
	         NED =  2 - (2 + dd)*exp(-dd) ,
	         SCHI2 =  1-(dd^2/(dd^2 +2)) )       

       if (raf=="HD" | raf=="NED") {
            ww <- (ww + 1)/(dd + 1)
       }
       ww[is.infinite(dd)] <- 0
       ww[ww > 1] <- 1
       ww[ww < 0] <- 0

       p <- ww%*%x/(sum(ww)*size)

       x.diff <- abs(p - p.old)
   }
#### end of while (x.diff > tol & iter < max.iter)

   if (iter < max.iter) {

   if (tot.sol==0) {
      p.store <- p
      w.store <- ww
      m.store <- mm
      f.store <- ff
      d.store <- dd
      tot.sol <- 1
   } else {
      if (min(abs(p.store-p))>equal) {
          p.store <- c(p.store,p)
          w.store <- rbind(w.store,ww)
          m.store <- rbind(m.store,mm)
          f.store <- rbind(f.store,ff)
          d.store <- rbind(d.store,dd)
          tot.sol <- tot.sol + 1
      }
   }

   } else not.conv <- not.conv + 1
   

}
##### end of while (tot.sol < num.sol & iboot < boot)

if (tot.sol) {
    result$p <- p.store
    result$tot.weights <- sum(ww)/nsize
    result$weights <- w.store
    result$delta <- d.store
    result$f.density <- f.store
    result$m.density <- m.store
    result$tot.sol <- tot.sol
    result$not.conv <- not.conv
    result$call <- match.call()
} else {
    if (verbose) cat("wle.binomial: No solutions are fuond, checks the parameters\n")
    result$p <- NA
    result$tot.weights <- NA
    result$weights <- rep(NA,nsize)
    result$delta <- rep(NA,nsize)
    result$f.density <- rep(NA,nsize)
    result$m.density <- rep(NA,nsize)
    result$tot.sol <- 0
    result$not.conv <- boot
    result$call <- match.call()
}

class(result) <- "wle.binomial"

return(result)
}

#############################################################
#                                                           #
#	print.wle.binomial function                         #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: August, 2, 2001                               #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

print.wle.binomial <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("p:\n")
    print.default(format(x$p, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("\nNumber of solutions ",x$tot.sol,"\n")
    cat("\n")
    invisible(x)
}



