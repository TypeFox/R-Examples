#########################################################
#                                                       #
#	wle.fracdiff function                           #
#	Author: Claudio Agostinelli                     #
#	E-mail: claudio@unive.it                        #
#	Date: December, 6, 2005                         #
#	Version: 0.2                                    #
#                                                       #
#	Copyright (C) 2002 Claudio Agostinelli          #
#                                                       #
#       Plus all the functions needed to                #
#       run a genetic algorithms                        #
#                                                       #
#########################################################

wle.fracdiff <- function(x, lower, upper, M, group, na.action=na.fail, tol=10^(-6), equal=10^(-3), raf="HD", smooth=0.0031, smooth.ao=smooth, boot=10, num.sol=1, x.init=rep(0,M), use.uniroot=FALSE, max.iter.out=20, max.iter.in=100, max.iter.step=5000, max.iter.start=max.iter.step,  verbose=FALSE, w.level=0.4, min.weights=0.5, init.values=NULL, num.max=length(x), include.mean=FALSE, ao.list=NULL, elitist=5, size.generation=5, size.population=10, type.selection="roulette", prob.crossover=0.8, prob.mutation=0.02, type.scale="none", scale.c=2) {

#    if (use.init) {
#        MM <- 0
#    } else {
#        MM <- M
#    }

    raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

    if (raf==-1) stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

    result <- list()
    series <- deparse(substitute(x))
    if(NCOL(x) > 1) stop("only implemented for univariate time series")

    x <- na.action(as.ts(x))
    nused <- length(x)

    if (length(x.init)!=M) stop("x.init must have M elements\n")

# start bootstrap iteration
    first.time <- TRUE
    iboot <- 1
    tot.sol <- 0
    not.conv <- 0

    while (iboot<=boot & tot.sol<num.sol) {
           pos.iboot <- round(runif(1,(group+1),nused))
           x.boot <- x[(pos.iboot-group+1):pos.iboot]

           if (!is.null(init.values)) {
               temp <- list(conv=TRUE)
               temp$d <- init.values[1]
               temp$sigma2 <- init.values[2]
               temp$x.mean <- init.values[3]
               temp$resid <- wle.fracdiff.residuals(d=temp$d, M=M, x=x, x.ao=x, x.init=x.init, x.mean=temp$x.mean)
           } else {
               temp <- wle.fracdiff.solve(x=x.boot, x.init=x.init, max.iter=max.iter.start, verbose=verbose, M=M,  lower, upper, tol=tol, use.uniroot=use.uniroot, include.mean=include.mean)
               temp$resid <- wle.fracdiff.residuals(temp$d, M=M, x=x, x.ao=x, x.init=x.init, x.mean=temp$x.mean)
               temp$sigma2 <- wle.fracdiff.sigma2(resid=temp$resid)
           }
    
           if (temp$conv) {
               d <- temp$d
               sigma2 <- temp$sigma2
               x.mean <- temp$x.mean
               resid <- temp$resid
               nresid <- length(resid)

               if (verbose) {
	           cat("Initial values from the subsample ",iboot,": \n parameters, d: ", d,"\n sigma2: ",sigma2," \n x.mean: ",x.mean, " \n") 
               }

               weights <- .Fortran("wlew",
    	            as.double(resid),
    	            as.integer(nresid),
	            as.double(resid),
	            as.integer(nresid),
	            as.integer(raf),
	            as.double(smooth.ao),
	            as.double(sigma2),
	            totweights=double(1),
	            weights=double(nresid),
				PACKAGE="wle")$weights

               if (sum(weights)/nresid >= min.weights) {
                   wres <- wle.fracdiff.ao(d=d, sigma2=sigma2, x=x, M=M, x.init=x.init, x.mean=x.mean, raf=raf, smooth=smooth.ao, w.level=w.level, verbose=verbose, ao.list=ao.list, num.max=num.max, elitist=elitist, size.generation=size.generation, size.population=size.population, type.selection=type.selection, prob.crossover=prob.crossover, prob.mutation=prob.mutation, type.scale=type.scale, scale.c=scale.c)                  

                   x.ao <- wres$x.ao
                   ao.position <- wres$ao.position
                   resid <- wle.fracdiff.residuals(d, M=M, x=x, x.ao=x.ao, x.init=x.init, x.mean=x.mean)
                    weights <- wle.weights(x=resid, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)$weights

                   if (!is.null(ao.position)) {
                       ao.list <- c(ao.list,list(ao.position))
                   }
                   ao.position.old <- c(ao.position,0)
                   conv <- TRUE
                   iter.out <- 0

                   while (!setequal(ao.position,ao.position.old) & conv) {
    	                  iter.out <- iter.out + 1
    	                  ao.position.old <- ao.position	
                          maxtol <- tol + 1
                          iter.in <- 0
                          while (maxtol > tol & conv) {
                                 iter.in <- iter.in + 1
                                 d.old <- d
                                 sigma2.old <- sigma2
                                 x.mean.old <- x.mean
	                         res <- wle.fracdiff.solve(x=x.ao, M=M, x.init=x.init, lower=lower, upper=upper, w=weights, tol=tol, max.iter=max.iter.step, verbose=verbose, use.uniroot=use.uniroot, include.mean=include.mean)
	                         d <- res$d
                                 x.mean <- res$x.mean
                                 resid <- wle.fracdiff.residuals(d=d, M=M, x=x, x.ao=x.ao, x.init=x.init, x.mean=x.mean)
	                         sigma2 <- wle.fracdiff.sigma2(resid=resid, w=weights)
                                 weights <- wle.weights(x=resid, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)$weights
                                 conv <- res$conv            

                                 if (iter.in > max.iter.in) {
                                     if (verbose) cat("Convergence problem: maximum iteration number reached in the outer loop\n")
                                     conv <- FALSE
                                 }
                                 maxtol <- max(abs(d-d.old),abs(sigma2-sigma2.old), abs(x.mean-x.mean.old))
                                  if (verbose)  {
	    	                      cat("inner loop, iteration: ",iter.in," \n parameters, d: ",d," \n sigma2: ",sigma2," \n x.mean: ",x.mean," \n")
	                          }
                          }

                          if (conv) {      
	                          if (verbose) {
	    	                      cat("outer loop, iteration: ",iter.out," convergence achieved for the inner loop \n")
	                          }

                              wres <- wle.fracdiff.ao(d=d, sigma2=sigma2, x=x, M=M, x.init=x.init, x.mean=x.mean, raf=raf, smooth=smooth.ao, w.level=w.level, verbose=verbose, ao.list=ao.list, num.max=num.max, elitist=elitist, size.generation=size.generation, size.population=size.population, type.selection=type.selection, prob.crossover=prob.crossover, prob.mutation=prob.mutation, type.scale=type.scale, scale.c=scale.c)         


                              x.ao <- wres$x.ao
                              ao.position <- wres$ao.position
                              resid <- wle.fracdiff.residuals(d, M=M, x=x, x.ao=x.ao, x.init=x.init, x.mean=x.mean)
                              sigma2 <- wle.fracdiff.sigma2(resid=resid, w=weights)
                               weights <- wle.weights(x=resid, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)$weights

                              if (!is.null(ao.position)) {
                                  ao.list <- c(ao.list,list(ao.position))
                              }

                              if ((lao <- length(ao.list))>2) {
                                  for (i.ao in 1:(lao-2)) {
                                       if (setequal(ao.list[[lao]],ao.list[[i.ao]])) {
                                           conv <- FALSE
                                       }
                                   } 
                               }

                               if (iter.out > max.iter.out) {
                                   if (verbose) cat("Convergence problem: maximum iteration number reached in the outer loop\n")
                                   conv <- FALSE
                               }

                          }

                   }
# end while (!setequal(ao.position,ao.position.old) & conv)

                   if (conv) {
                       resid.with.ao <- wle.fracdiff.residuals(d, M=M, x=x, x.ao=x, x.init=x.init, x.mean=x.mean) 
                       resid.with.ao <- ts(resid.with.ao, start=(start(x)), end=end(x), frequency=frequency(x))
                       class(resid.with.ao) <- "ts"
                       weights.with.ao <- wle.weights(x=resid.with.ao, smooth=smooth, sigma2=sigma2, raf=raf, tol=tol, location=TRUE)$weights

                       resid <- ts(resid, start=(start(x)), end=end(x), frequency=frequency(x))        
                       class(resid) <- "ts" 

                       resid.without.ao <- wle.fracdiff.residuals(d, M=M, x=x.ao, x.ao=x.ao, x.init=x.init, x.mean=x.mean) 
                       resid.without.ao <- ts(resid.without.ao, start=(start(x)), end=end(x), frequency=frequency(x))
                       class(resid.without.ao) <- "ts"
                       weights.without.ao <- wle.weights(x=resid.without.ao, smooth=smooth, sigma2=sigma2, raf=raf, tol=tol, location=TRUE)$weights
 
                       x.ao <- ts(x.ao, start=start(x), end=end(x), frequency=frequency(x))
                       class(x.ao) <- "ts" 

    if (first.time) {
        d.final <- d
        sigma2.final <- sigma2
        x.mean.final <- c(x.mean)
        weights.final <- weights
        weights.with.ao.final <- weights.with.ao
        weights.without.ao.final <- weights.without.ao   
        resid.final <- resid
        resid.with.ao.final <- resid.with.ao
        resid.without.ao.final <- resid.without.ao
        x.ao.final <- x.ao
        ao.position.final <- list(ao.position)
        first.time <- FALSE	
        tot.sol <- 1
    } else {
        if (min(abs(d.final-d))>equal) {
	    tot.sol <- tot.sol+1
	    d.final <- c(d.final,d)
            sigma2.final <- c(sigma2.final,sigma2)
            x.mean.final <- c(x.mean.final,c(x.mean))
	    weights.final <- rbind(weights.final,weights)
            weights.with.ao.final <- rbind(weights.with.ao.final, weights.with.ao)
            weights.without.ao.final <- rbind(weights.without.ao.final, weights.without.ao)
	    resid.final <- rbind(resid.final,resid)
            resid.with.ao.final <- rbind(resid.with.ao.final,resid.with.ao)
            resid.without.ao.final <- rbind(resid.without.ao.final,resid.without.ao)
            ao.position.final <- c(ao.position.final,list(ao.position))
            x.ao.final <- rbind(x.ao.final,x.ao)
        }
    }

                   } else {
    	               not.conv <- not.conv+1
                   }
               } else {
    	           not.conv <- not.conv+1
               }
           } else {
    	       not.conv <- not.conv+1
           }    

        iboot <- iboot+1
    }
# end bootstrap iteration

if(tot.sol==0) {
   
    result$d <- NULL
    result$sigma2 <- NULL
    result$x.mean <- NULL
    result$resid <- NULL
    result$resid.with.ao <- NULL
    result$resid.without.ao <- NULL
    result$x.ao <- NULL
    result$call <- match.call()
    result$weights <- NULL
    result$weights.with.ao <- NULL
    result$weights.without.ao <- NULL   
    result$tot.sol <- 0
    result$not.conv <- not.conv
    result$ao.position <- NULL
} else { 
    result$d <- d.final
    result$sigma2 <- sigma2.final
    result$x.mean <- x.mean.final
    result$resid <- resid.final
    result$resid.without.ao <- resid.without.ao.final
    result$resid.with.ao <- resid.with.ao.final
    result$x.ao <- x.ao.final
    result$call <- match.call()
    result$weights <- weights.final
    result$weights.with.ao <- weights.with.ao.final
    result$weights.without.ao <- weights.without.ao.final
    result$tot.sol <- tot.sol
    result$not.conv <- not.conv
    result$ao.position <- ao.position.final
    }

    class(result) <- "wle.farima"
    
return(result)
}

#############################################################
#                                                           #
#	wle.fracdiff.solve function                         #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@stat.unipd.it                       #
#	Date: December, 11, 2001                            #
#	Version: 0.1-1                                      #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.solve <- function(x, M=100, x.init=rep(0,M), lower, upper, w=rep(1,length(x)), tol=.Machine$double.eps^0.25, max.iter=1000, verbose=FALSE, use.uniroot=FALSE, include.mean=FALSE) {

    result <- list()
    if (include.mean) {
        x.mean <- w%*%x/sum(w)
        x <- x - x.mean
    } else {
        x.mean <- 0
    }

    result$x.mean <- x.mean
    result$call <- match.call()

    if (use.uniroot) {
        temp <- uniroot(wle.fracdiff.equation, x=x, M=M, w=w, x.init=x.init, lower=lower, upper=upper, use.uniroot=use.uniroot, tol=tol, maxiter=max.iter, verbose=verbose)
        if (temp$iter < max.iter) {
            result$d <- temp$root
            result$conv <- TRUE
        } else {
            result$d <- NA
            result$conv <- FALSE
        }
    } else {
        temp <- optimize(wle.fracdiff.equation, x=x, M=M, w=w, x.init=x.init, lower=lower, upper=upper, use.uniroot=use.uniroot, tol=tol, verbose=verbose)
        result$d <- temp$minimum
        result$conv <- TRUE
    }
    return(result)
}

#############################################################
#                                                           #
#	wle.fracdiff.equation function                      #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 5, 2005                             #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################


wle.fracdiff.equation <- function(d, M, x, x.init=rep(0,M), w=rep(1,length(x)), use.uniroot=FALSE, verbose=FALSE) {
 
    nused <- length(x)
    pi.coef <- wle.fracdiff.pi.coef(d,M)
    if (use.uniroot) {
        xi.coef <- wle.fracdiff.xi.coef(d,M)
    }
    y <- c(x.init,x)

    somma <- 0
    if (use.uniroot) {
        for (k in 1:nused) {
             somma <- somma + w[k]*(x[k]+pi.coef%*%y[(k-1+M):k])*(xi.coef%*%y[(k-1+M):k])
        }
    } else {
        for (k in 1:nused) {
             somma <- somma + w[k]*(x[k]+pi.coef%*%y[(k-1+M):k])^2
        }
    }

    if (verbose) cat("value of d: ",d," value of the function: ",somma,"\n")

    return(as.vector(somma))
}


#############################################################
#                                                           #
#	wle.fracdiff.pi.coef function                       #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: November, 30, 2001                            #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.pi.coef <- function(d,M) {
     pi.coef <- rep(0,M)
     pi.coef[1] <- -d
     for (j in 2:M) {
          pi.coef[j] <- pi.coef[j-1]*(j-1-d)/j
     }
return(pi.coef)
}

#############################################################
#                                                           #
#	wle.fracdiff.xi.coef function                       #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: March, 4, 2009                                #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 2009 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.xi.coef <- function(d,M) {
     xi.coef <- rep(0,M)
     for (j in 1:M) {
          primo.termine <- gamma(j-d)/gamma(j+1)
          if (is.nan(primo.termine)) {
              primo.termine <- j^(-(1+d))*exp(d)
          }
          xi.coef[j] <- primo.termine*(digamma(j-d) - digamma(-d))
     }
     xi.coef <- xi.coef*d/gamma(1-d)
return(xi.coef)
}

#############################################################
#                                                           #
#	wle.fracdiff.residuals function                     #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 11, 2001                            #
#	Version: 0.1-1                                      #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.residuals <- function(d, M, x, x.ao, x.init=rep(0,M), x.mean=0) {

    x <- x - x.mean
    x.ao <- x.ao - x.mean

    nused <- length(x) 
    pi.coef <- wle.fracdiff.pi.coef(d,M)
    y <- c(x.init,x.ao)
    resid <- rep(0,nused)

    for (t in 1:nused) {
         resid[t] <- x[t]+pi.coef%*%y[(t-1+M):t]
    }
    resid <- resid[1:nused]

    return(resid)
}

#############################################################
#                                                           #
#	wle.fracdiff.sigma2 function                        #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: February, 23, 2009                            #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.sigma2 <- function(resid, w=rep(1,length(resid))) {
    sigma2 <- sum(w*resid^2)/(sum(w) - 2)
    return(sigma2)
}

#############################################################
#                                                           #
#	wle.fracdiff.fitted function                        #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 5, 2001                             #
#	Version: 0.1-1                                      #
#                                                           #
#	Copyright (C) 2001 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.fitted <- function(t, d, M, x, x.init=rep(0,M), x.mean=0) {
 
    x <- x - x.mean
    nused <- length(x)
    pi.coef <- wle.fracdiff.pi.coef(d,M)
    y <- c(x.init,x)
    return((-pi.coef%*%y[(t-1+M):t]+x.mean))
}

#############################################################
#                                                           #
#	wle.crossover.ao function                           #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: October, 13, 2012                             #
#	Version: 0.1-1                                      #
#                                                           #
#	Copyright (C) 2012 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.crossover.ao <- function(x, y, prob.crossover) {
    size <- length(x)
    if (rbinom(n=1, size=1, prob=prob.crossover)) {
		split <- sample(x=1:(size-1), size=1, replace=FALSE)
		x.temp <- c(x[1:split],y[(split+1):size])
        y.temp <- c(y[1:split],x[(split+1):size])
    } else {
        x.temp <- x
        y.temp <- y
    }

    result <- list(x=x.temp, y=y.temp)
    return(result)
}

#############################################################
#                                                           #
#	wle.mutation.ao function                            #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: October 13, 2012                              #
#	Version: 0.1-1                                      #
#                                                           #
#	Copyright (C) 2012 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.mutation.ao <- function(x, prob.mutation) {
   mutation <- as.logical(rbinom(n=length(x), size=1, prob=prob.mutation))
   replace <- sample(x=c(0, 1), size=sum(mutation), replace=TRUE)   
   x[mutation] <- replace
   return(x)
}


#############################################################
#                                                           #
#	wle.selection.ao function                           #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 5, 2005                             #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.selection.ao <- function(population, type.selection, fitness) {
   n <- nrow(population)
   if (missing(fitness)) fitness <- rep(1,n)
   if (type.selection=="uniform") { 
       pos <- sample(x=1:n, size=2, replace=TRUE)   
   } else {
       pos <- sample(x=(1:n), size=2, replace=TRUE, prob=(fitness/sum(fitness)))
   }
   x <- population[pos[1],] 
   y <- population[pos[2],]
   result <- list(x=x, y=y)
   return(result)
}


#############################################################
#                                                           #
#	wle.fitness.population.ao function                  #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 5, 2005                             #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fitness.population.ao <- function(population, x, ao, d, M, x.init, x.mean, nresid, smooth, sigma2, raf) {
   size.population <- nrow(population)
   fitness.value <- rep(0, size.population)
   for (i in 1:size.population)  {
        decode <- wle.decode.ao(x=population[i,], ao=ao)
        fitness.value[i] <- wle.fitness.ao(x=decode, serie=x, d=d, M=M, x.init=x.init, x.mean=x.mean, nresid=nresid, smooth=smooth, sigma2=sigma2, raf=raf) 
   }
   return(fitness.value)
}

#############################################################
#                                                           #
#	wle.generate.population.ao function                 #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 5, 2005                             #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.generate.population.ao <- function(length.string, size.population) {
   population <- matrix(sample(x=c(0, 1), size=(length.string*size.population), replace=TRUE), ncol=length.string)
   return(population)
}

#############################################################
#                                                           #
#	wle.ga.ao function                                  #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 5, 2005                             #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.ga.ao <- function(x, ao, d, elitist=5, size.generation=100, size.population=50, type.selection="roulette", prob.crossover=0.8, prob.mutation=0.02, type.scale="none", scale.c=2, M, x.init, x.mean, nresid, smooth, sigma2, raf, ao.list=NULL) {
   population <- wle.generate.population.ao(length.string=length(ao), size.population=size.population)
   if (!is.null(ao.list)) {
       population <- rbind(ao.list, population)
       size.population <- nrow(population)
   }
   fit.of.best <- 0
   
   if (elitist > 0) {    
       if (elitist>=1) {
            elitist <- min(floor(elitist), size.population)
       } else {
            elitist <- floor(elitist*size.population)
       }
   } else {
       elitist <- 0
   }

   for (j in 1:(size.generation+1)) {
        new.population <- vector(length=0)
        fit <- wle.fitness.population.ao(population=population, x=x, ao=ao, d=d, M=M, x.init=x.init, x.mean=x.mean, nresid=nresid, smooth=smooth, sigma2=sigma2, raf=raf)
        if (min(fit) < 0) stop("fitness function can not be negative")

        if (elitist) {
            telitist <- rev(order(fit))[1:elitist]
            new.population <- population[telitist,]
        }
        
        mfit <- max(fit)
        best.pos <- which(fit==mfit)[1]

        if (mfit > fit.of.best) {
            best.of.best <- population[best.pos,]
            fit.of.best <- mfit
        }

        fit.scale <- fit
        if (type.scale!="none") fit.scale <- wle.fitness.scale.ao(x=fit, type=type.scale, scale.c=scale.c) 
        i <- elitist
        while (i < size.population) {
             sel <- wle.selection.ao(population, type.selection=type.selection, fitness=fit.scale)
             if (length(ao) > 1) {
                 cross <- wle.crossover.ao(x=sel$x, y=sel$y, prob.crossover=prob.crossover)
             } else {
                 cross <- list(x=sel$x, y=sel$y) 
             }    
             mut.x <- wle.mutation.ao(cross$x, prob.mutation=prob.mutation)
	     mut.y <- wle.mutation.ao(cross$y, prob.mutation=prob.mutation)
             mut.x <- as.vector(mut.x)
             mut.y <- as.vector(mut.y)
             new.population <- rbind(new.population, mut.x, mut.y)
             i <- nrow(new.population)
        }
        population <- new.population[1:size.population,]
        if (length(ao)==1) population <- matrix(population, ncol=1)
   }

   return(list(ao=wle.decode.ao(x=best.of.best, ao=ao), fit=fit.of.best))
}

#############################################################
#                                                           #
#	wle.fitness.scale.ao function                       #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 5, 2005                             #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fitness.scale.ao <- function(x, type="linear", scale.c=2, tol=0.001) {
   if (type=="linear") {	
       media <- mean(x)
       if ((media-min(x)<tol) | ((max(x)-media < tol))) return(x) 
       scale.min <- media/(media-min(x))
       scale <- (scale.c-1)*media/(max(x)-media)
       if (scale > scale.min) scale <- scale.min
       x <- scale*x+media*(1-scale)
   } else {
       if (type=="sigma.truncation") {
           x <- x - (mean(x) - scale.c*sqrt(var(x)))

       }
   }
   x[x <0] <- 0
   return(x)
}

#############################################################
#                                                           #
#	wle.encode.ao function                              #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 5, 2005                             #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.encode.ao <- function(x, ao) {
   y <- rep(0, length(ao))
   for (i in 1:length(x)) {
        y[ao==x[i]] <- 1
   }
   return(y)
}

#############################################################
#                                                           #
#	wle.decode.ao function                              #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 5, 2005                             #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.decode.ao <- function(x, ao) {
   return(ao[x==1])
}

#############################################################
#                                                           #
#	wle.fitness.ao function                             #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 5, 2005                             #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fitness.ao <- function(x, serie, d, M, x.init, x.mean, nresid, smooth, sigma2, raf) {
              x.ao <- serie
              x <- sort(x)              
              for (t in x) {
                   x.ao[t] <- wle.fracdiff.fitted(t=t, d=d, M=M, x=x.ao, x.init=x.init, x.mean=x.mean)
              }
              
              resid.ao <- wle.fracdiff.residuals(d=d, M=M, x=serie, x.ao=x.ao, x.init=x.init, x.mean=x.mean)
              if (length(x)) {
                  resid.ao <- resid.ao[-x]
              }
              ww <- wle.weights(x=resid.ao, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)$weights
              return(sum(ww)/nresid)
}

#############################################################
#                                                           #
#	wle.fracdiff.ao function                            #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: December, 5, 2005                             #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2005 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.fracdiff.ao <- function(d, sigma2, x, M=100, x.init=rep(0,M), x.mean=0, raf=1, smooth=0.0031, w.level=0.5, verbose=FALSE, ao.list=NULL, num.max=length(x), elitist=5, size.generation=100, size.population=50, type.selection="roulette", prob.crossover=0.8, prob.mutation=0.02, type.scale="none", scale.c=2) {

    nused <- length(x)
    resid <- wle.fracdiff.residuals(d=d, M=M, x=x, x.ao=x, x.init=x.init, x.mean=x.mean)  
    nresid <- length(resid)

    weights <- .Fortran("wlew",
	as.double(resid), 
	as.integer(nresid),
	as.double(resid), 
	as.integer(nresid), 
	as.integer(raf),
	as.double(smooth),
	as.double(sigma2),
	totweight=double(1),
	weights=double(nresid),
	PACKAGE="wle")$weights

    ao.position <- NULL
    pos.temp <- 1:nresid
    pos.temp <- pos.temp[rev(order(weights))]
    weights.sort <- rev(sort(weights))
    ao.temp <- weights.sort <= w.level
    pos.temp <- pos.temp[ao.temp]

    ao <- rep(FALSE,nused)
    if (length(pos.temp)) {
        pos.temp <- pos.temp[1:min(length(pos.temp),num.max)]
        ao[pos.temp] <- TRUE
    }

    pos <- which(ao)

    if (verbose) {
        cat("We have the following observations under the w.level=",w.level,":\n",pos,"\n")
    }

    if (length(pos)) {

      if (!is.null(ao.list)) {
          temp <- matrix(0, ncol=length(pos), nrow=length(ao.list))
          for (ilist in 1:length(ao.list)) {
               temp[ilist,] <- wle.encode.ao(intersect(ao.list[[ilist]], pos), pos)
          }
          ao.list <- temp
      }
      
      ga.result <- wle.ga.ao(x=x, ao=pos, d=d, elitist=5, size.generation=size.generation, size.population=size.population, type.selection=type.selection, prob.crossover=prob.crossover, prob.mutation=prob.mutation, type.scale=type.scale, scale.c=scale.c, M=M, x.init=x.init, x.mean=x.mean, nresid=nresid, smooth=smooth, sigma2=sigma2, raf=raf, ao.list=ao.list)        
         if (ga.result$fit<(sum(weights)/nresid)) {
             ao.position <- NULL
         } else {
             ao.position <- ga.result$ao 
         }

    } else {
        ao.position <- NULL
    }

    x.ao <- x
    for (t in ao.position) {
         x.ao[t] <- wle.fracdiff.fitted(t=t, d=d, M=M, x=x.ao, x.init=x.init, x.mean=x.mean)
    }

    resid.ao <- wle.fracdiff.residuals(d=d, M=M, x=x, x.ao=x.ao, x.init=x.init, x.mean=x.mean)
    w.temp <- wle.weights(x=resid.ao, smooth=smooth, sigma2=sigma2, raf=raf, location=TRUE)
    resid.ao <- resid.ao - w.temp$location

    if (verbose) {
        cat("Additive outliers: \n", ao.position, "\n")
    }

    return(list(x.ao=x.ao, resid.ao=resid.ao, ao.position=ao.position))
}

