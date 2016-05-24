#############################################################
#                                                           #
#	wmme.gamma function                                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: March 25, 2011                                #
#	Version: 0.3                                        #
#                                                           #
#	Copyright (C) 2011 Claudio Agostinelli              #
#                                                           #
#############################################################
## WARNING: MAX shape allowed is 70! This was true with 0.2-3 version. We did not check for the new version.
## SCHI2 work much better than HD, this is way it is the default!!!!!!
#######


wmme.gamma <- function(x, boot=30, group, num.sol=1, raf="SCHI2", smooth=0.008, tol=10^(-6), equal=10^(-3), max.iter=500, use.smooth=TRUE, tol.int, verbose=FALSE) {

max.aa <- 70
x <- as.vector(x)
size <- length(x)
result <- list()
  
#### weighted central moments
m1 <- function(x, w=rep(1, length(x))) c(w%*%x/sum(w))
m2 <- function(x, m1, w=rep(1, length(x))) c(w%*%(x-m1)^2/sum(w))
m3 <- function(x, m1, w=rep(1, length(x))) c(w%*%(x-m1)^3/sum(w))

#### wmme solutions
alphap <- function(m2, m3) 4*m2^3/m3^2          #shape
betap <- function(m2, m3) m3/(2*m2)             #scale
gammap <- function(m1, m2, m3) (m1 - 2*m2^2/m3) #location

raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

if (missing(group)) {
group <- 0
}

if (size<2) {
stop("Number of observation must be at least equal to 2")
}

if (group<2) {
    group <- max(round(size/4),2)
    if (verbose) cat("wmme.gamma: dimension of the subsample set to default value: ",group,"\n")
}

maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

if (boot<1 | log(boot) > maxboot) {
    stop("Bootstrap replication not in the range")
}

if (!(num.sol>=1)) {
    if (verbose) cat("wmme.gamma: number of solution to report set to 1 \n")
    num.sol <- 1
}

if (max.iter<1) {
    if (verbose) cat("wmme.gamma: max number of iterations set to 500 \n")
    max.iter <- 500
}

if (smooth<10^(-5)) {
    if (verbose) cat("wmme.gamma: the smooth parameter seems too small \n")
}

if (tol<=0) {
    if (verbose) cat("wmme.gamma: the accuracy must be positive, using default value: 10^(-6) \n")
    tol <- 10^(-6)
}

if (equal<=tol) {
    if (verbose) cat("wmme.gamma: the equal parameter must be positive, using default value: tol+10^(-3) \n")
    equal <- tol+10^(-3)
}

if (!is.logical(use.smooth)) {
    if (verbose) cat("wmme.gamma: the use.smooth must be a logical value, using default value \n")
    use.smooth <- TRUE
}

if (missing(tol.int)) {
   tol.int <- tol*10^(-4)
} else {
   if (tol.int <=0) {
       if (verbose) cat("wmme.gamma: tol.int must be positive, using default value \n")
       tol.int <- tol*10^(-4) 
   } 
}

tot.sol <- 0
not.conv <- 0
iboot <- iter <- 0
while (tot.sol < num.sol & iboot <= boot) {

  if (verbose) cat('Solution ', tot.sol+1, '\n')
  if (verbose) cat('Boot ', iboot, '\n')

  cont <- TRUE
  jj <- 0
  while (cont & jj <= 100) {
    jj <- jj + 1
    x.boot <- x[sample(1:size, group, replace = FALSE)]
    m1.boot <- m1(x.boot)
    m2.boot <- m2(x.boot, m1=m1.boot)
    m3.boot <- m3(x.boot, m1=m1.boot)
    aa <- alphap(m2.boot, m3.boot)
    bb <- betap(m2.boot, m3.boot)
    cc <- gammap(m1.boot, m2.boot, m3.boot)
    if (all(!is.na(c(aa, bb, cc))) && aa < max.aa) cont <- FALSE
  }
  iboot <- iboot + 1
  if (verbose) {
    cat('a', aa, '\n')
    cat('b', bb, '\n')
    cat('c', cc, '\n')
  }
  
  if (all(!is.na(c(aa, bb, cc))) && aa < max.aa) {
    xdiff <- tol + 1
    iter <- 0
    while (xdiff > tol & iter < max.iter+2) {
       iter <- iter + 1
       temp <- aa*bb^2
       xx <- x-cc+10^(-2)
       dsup <- max(xx)+ 3*smooth*temp
       aaold <- aa
       bbold <- bb
       ccold <- cc

#       if (verbose) {
#          cat('a', aa, '\n')
#          cat('b', bb, '\n')
#          cat('c', cc, '\n')
#       }

       z <- .Fortran("wlegamma",
	    as.double(xx),
	    as.double(xx),
	    as.integer(size),
	    as.integer(size),                     
	    as.integer(raf),
            as.double(1),
	    as.double(smooth*temp),
            as.integer(1*use.smooth),
            as.double(dsup),
	    as.double(tol),
            as.double(tol.int),
	    as.double(1/bb),
	    as.double(aa),
	    weights=double(size),
	    density=double(size),
	    model=double(size),
        PACKAGE = "wle")

       ww <- z$weights
####       if (verbose) print(summary(ww))
       m1n <- m1(x, w=ww)
       m2n <- m2(x, m1=m1n, w=ww)
       m3n <- m3(x, m1=m1n, w=ww)
       aa <- alphap(m2n, m3n)
       bb <- betap(m2n, m3n)
       cc <- gammap(m1n, m2n, m3n)
       xdiff <- max(abs(c(aa-aaold,bb-bbold,cc-ccold)))
       if (sum(ww)/length(ww) < 0.1) iter <- max.iter+2
       if (is.na(aa) || aa > max.aa) iter <- max.iter+2       
   }

  if (verbose) {
    cat('a', aa, '\n')
    cat('b', bb, '\n')
    cat('c', cc, '\n')
  }
    
   if (iter < max.iter) {

   if (tot.sol==0) {
      a.store <- aa
      b.store <- bb
      c.store <- cc
      w.store <- ww
      f.store <- z$density
      m.store <- z$model
      d.store <- f.store/m.store - 1
      tot.sol <- 1
   } else {
      if (min(abs(a.store-aa))>equal & min(abs(b.store-bb))>equal & min(abs(c.store-cc))>equal) {
          a.store <- c(a.store,aa)
          b.store <- c(b.store,bb)
          c.store <- c(c.store,cc)
          w.store <- rbind(w.store,ww)
          f.store <- rbind(f.store,z$density)
          m.store <- rbind(m.store,z$model)
          d.store <- rbind(d.store,z$density/z$model - 1)
          tot.sol <- tot.sol + 1
      }
   }
   
   } else not.conv <- not.conv + 1
 } else not.conv <- not.conv + 1

}
##### end of while (tot.sol < num.sol & iboot < boot)

if (tot.sol) {
   result$scale <- c(b.store)
   result$rate <- c(1/b.store)  
   result$shape <- c(a.store)
   result$location <- c(c.store)
   if (tot.sol>1) {
       tot.w <- apply(w.store,1,sum)/size
   } else tot.w <- sum(w.store)/size
  
   result$tot.weights <- tot.w
   result$weights <- w.store
   result$delta <- d.store
   result$f.density <- f.store
   result$m.density <- m.store
   result$tot.sol <- tot.sol
   result$not.conv <- not.conv
   result$call <- match.call()
} else{
   if (verbose) cat("wmme.gamma: No solutions are fuond, checks the parameters\n")
   result$scale <- NA
   result$rate <- NA
   result$shape <- NA
   result$location <- NA
   result$tot.weights <- NA
   result$weights <- rep(NA,size)
   result$delta <- rep(NA,size)
   result$f.density <- rep(NA,size)
   result$m.density <- rep(NA,size)
   result$tot.sol <- 0
   result$not.conv <- boot
   result$call <- match.call()
}

   class(result) <- "wmme.gamma"
   return(result)
}

#############################################################
#                                                           #
#	print.wmme.gamma function                           #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: April, 18, 2007                               #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2007 Claudio Agostinelli              #
#                                                           #
#############################################################

print.wmme.gamma <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Scale:\n")
    print.default(format(x$scale, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Rate:\n")
    print.default(format(x$rate, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")  
    cat("Shape:\n")
    print.default(format(x$shape, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")  
    cat("Location:\n")
    print.default(format(x$location, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("\nNumber of solutions ",x$tot.sol,"\n")
    cat("\n")
    invisible(x)
}






