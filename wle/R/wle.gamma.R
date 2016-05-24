#############################################################
#                                                           #
#	wle.gamma function                                  #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: February, 21, 2011                            #
#	Version: 0.4-3                                      #
#                                                           #
#	Copyright (C) 2011 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.gamma <- function(x, boot=30, group, num.sol=1, raf="HD", smooth=0.008, tol=10^(-6), equal=10^(-3), max.iter=500, shape.int=c(0.01, 100), use.smooth=TRUE, tol.int, verbose=FALSE, maxiter=1000) {

sem <- options()$show.error.messages

wsolve <- function (o, media, medialog) {
   medialog + log(o/media) - digamma(o)
}


## the wlegamma fortran function implements others type of RAF too.
raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

if (missing(group)) {
group <- 0
}

x <- as.vector(x)
size <- length(x)
result <- list()

if (size<2) {
stop("Number of observation must be at least equal to 2")
}

if (group<2) {
    group <- max(round(size/4),2)
    if (verbose) cat("wle.gamma: dimension of the subsample set to default value: ",group,"\n")
}

maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

if (boot<1 | log(boot) > maxboot) {
    stop("Bootstrap replication not in the range")
}

if (!(num.sol>=1)) {
    if (verbose) cat("wle.gamma: number of solution to report set to 1 \n")
    num.sol <- 1
}

if (max.iter<1) {
    if (verbose) cat("wle.gamma: max number of iterations set to 500 \n")
    max.iter <- 500
}

if (maxiter<1) {
    if (verbose) cat("wle.gamma: max number of iterations for the uniroot function set to 1000 \n")
    maxiter <- 1000
}

if (smooth<10^(-5)) {
    if (verbose) cat("wle.gamma: the smooth parameter seems too small \n")
}

if (tol<=0) {
    if (verbose) cat("wle.gamma: the accuracy must be positive, using default value: 10^(-6) \n")
    tol <- 10^(-6)
}

if (equal<=tol) {
    if (verbose) cat("wle.gamma: the equal parameter must be positive, using default value: tol+10^(-3) \n")
    equal <- tol+10^(-3)
}

if (!is.logical(use.smooth)) {
    if (verbose) cat("wle.gamma: the use.smooth must be a logical value, using default value \n")
    use.smooth <- TRUE
}

if (length(shape.int)!=2) stop("shape.int must be a vector of length 2 \n")

shape.int <- sort(shape.int, decreasing = FALSE)

if (shape.int[2] <= 0) {
    stop("the elements of shape.int must be positive \n")
}

if (shape.int[1] <= 0) {
    if (verbose) cat("wle.gamma: the elements of shape.int must be positive, using default value \n")
    shape.int[1] <- tol
}

if (missing(tol.int)) {
   tol.int <- tol*10^(-4)
} else {
   if (tol.int <=0) {
       if (verbose) cat("wle.gamma: tol.int must be positive, using default value \n")
       tol.int <- tol*10^(-4) 
   } 
}

tot.sol <- 0
not.conv <- 0
iboot <- 0

xlog <- log(x)

while (tot.sol < num.sol & iboot <= boot) {
   cont <- TRUE
   i <- 0
   while (cont & iboot <= boot) {
          i <- i + 1
          iboot <- iboot + 1
          x.boot <- x[sample(1:size, group, replace = FALSE)]
          xlog.boot <- log(x.boot)

          media <- sum(x.boot)/group
          medialog<- sum(xlog.boot)/group

          options(show.error.messages=FALSE)
          o <- try(uniroot(wsolve, interval=shape.int, media=media, medialog=medialog, tol=tol, maxiter=maxiter)$root)
          options(show.error.messages=sem)
          
          if (!is.character(o)) {
              cont <- FALSE
          }
   }

   if (!is.character(o)) {
   
       if (o < tol) o <- 2*tol

       l <- media/o

       xdiff <- tol + 1
       iter <- 0
       while (xdiff > tol & iter < max.iter) {

       iter <- iter + 1
       ### shape
       shape <- o
       ### rate
       lambda <- 1/l
       temp <- shape/lambda^2
       dsup <- max(x)+ 3*smooth*temp

       z <- .Fortran("wlegamma",
	    as.double(x),
            as.double(x),
	    as.integer(size),
	    as.integer(size),                     
	    as.integer(raf),
            as.double(1),
	    as.double(smooth*temp),
            as.integer(1*use.smooth),
            as.double(dsup),
	    as.double(tol),
            as.double(tol.int),
	    as.double(lambda),
	    as.double(shape),
	    weights=double(size),
	    density=double(size),
	    model=double(size),
        PACKAGE = "wle")

       ww <- z$weights
       wsum <- sum(ww)
       wmedia <- ww%*%x/wsum
       wmedialog <- ww%*%xlog/wsum

       options(show.error.messages=FALSE)
       o <- try(uniroot(wsolve, interval=shape.int, media=wmedia, medialog=wmedialog, tol=tol, maxiter=maxiter)$root)
       options(show.error.messages=sem)
       
       if (!is.character(o)) {
           if (o < tol) o <- 2*tol

           l <- wmedia/o

           xdiff <- max(abs(c(o-shape,l-1/lambda)))
       } else {
           xdiff <- 0
           iter <- max.iter+1
       }
   }

   if (iter <= max.iter) {

   if (tot.sol==0) {
      o.store <- o
      l.store <- l
      w.store <- ww
      f.store <- z$density #### /sqrt(2*pi*smooth*temp)
      m.store <- z$model #### /sqrt(2*pi*smooth*temp)
      d.store <- f.store/m.store - 1
      bw.store <- c(smooth*shape/lambda^2)
      tot.sol <- 1
   } else {
      if (min(abs(o.store-o))>equal & min(abs(l.store-l))>equal) {
          o.store <- c(o.store,o)
          l.store <- c(l.store,l)
          w.store <- rbind(w.store,ww)
          f.store <- rbind(f.store,z$density) ##### /sqrt(2*pi*smooth*temp))
          m.store <- rbind(m.store,z$model) ##### /sqrt(2*pi*smooth*temp))
          d.store <- rbind(d.store,z$density/z$model - 1)
          bw.store <- c(bw.store, smooth*shape/lambda^2)
          tot.sol <- tot.sol + 1
      }
   }

   } else not.conv <- not.conv + 1
   } else not.conv <- not.conv + i

}
##### end of while (tot.sol < num.sol & iboot < boot)

if (tot.sol) {
   result$scale <- c(l.store)
   result$rate <- c(1/l.store)  
   result$shape <- o.store
   
   if (tot.sol>1) {
       tot.w <- apply(w.store,1,sum)/size
   } else tot.w <- sum(w.store)/size
  
   result$tot.weights <- tot.w
   result$weights <- w.store
   result$delta <- d.store
   result$f.density <- f.store
   result$m.density <- m.store
   result$bw <- c(bw.store)
   result$tot.sol <- tot.sol
   result$not.conv <- not.conv
   result$call <- match.call()
} else{
   if (verbose) cat("wle.gamma: No solutions are fuond, checks the parameters\n")
   result$scale <- NA
   result$rate <- NA
   result$shape <- NA
   result$tot.weights <- NA
   result$weights <- rep(NA,size)
   result$delta <- rep(NA,size)
   result$f.density <- rep(NA,size)
   result$m.density <- rep(NA,size)
   result$tot.sol <- 0
   result$not.conv <- boot
   result$call <- match.call()
}

   class(result) <- "wle.gamma"
   return(result)
}

#############################################################
#                                                           #
#	print.wle.gamma function                                #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: August, 28, 2003                                  #
#	Version: 0.3                                            #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.wle.gamma <- function(x, digits = max(3, getOption("digits") - 3), ...) {
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
    cat("\nNumber of solutions ",x$tot.sol,"\n")
    cat("\n")
    invisible(x)
}






