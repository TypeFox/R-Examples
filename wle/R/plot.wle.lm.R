#############################################################
#                                                           #
#	plot.wle.lm function                                    #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: April, 17, 2003                                   #
#	Version: 0.5-1                                          #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################
#### fix a problem with iid when id.n=0!!!

plot.wle.lm <- function(x, roots, which=1:4, which.main, level.weight=0.5, ask=dev.interactive(), col=c(2, 1, 3), id.n=3, labels.id, cex.id = 0.75, verbose=FALSE, ...) {

   old.par <- par(no.readonly=TRUE)
   on.exit(par(old.par))

   if (!inherits(x, "wle.lm")) stop("Invalid 'wle.lm' object")
   if (ask) par(ask = TRUE)
   if (!is.numeric(which) || any(which < 0) || any(which > 4))
       stop("`which' must be in 0:4")

   if (level.weight < 0 | level.weight > 1) {
       if (verbose) cat("plot.wle.lm: level.weight should be between zero and one, set to 0.5 \n")
           level.weight <- 0.5
   }

   param <- x$coefficients
   res <- x$residuals
   y.fit <- x$fitted.values
   weight <- x$weights
   tot.weight <- x$tot.weights
   tot.sol <- x$tot.sol
   if (tot.sol > 1) {
       size <- ncol(y.fit)
       nomi <- dimnames(res)[[2]]
   } else {
       size <- length(y.fit)
       nomi <- names(res)
   }

   if (is.null(id.n)) {
       id.n <- 0
   } else {
       id.n <- as.integer(id.n)
       if (id.n < 0 || id.n > size)
	   stop(paste("`id.n' must be in { 1,..,", size,"}"))
   }
   if(id.n > 0) {
        if (missing(labels.id)) {
            if (is.null(nomi)) {
                labels.id <- paste(1:size)
            } else {
                labels.id <- nomi
            }
        } else {
            if (length(labels.id)!=size) {
                stop("the length of 'labels.id' must be equal to the number of observations")
            }
        }
        iid <- 1:id.n
        text.id <- function (x, y, ind, labels.id, adj.x = FALSE)
            text(x - if(adj.x) strwidth(" ")*cex.id else 0, y, labels.id[ind],
                 cex = cex.id, xpd = TRUE, adj = if(adj.x) 1)
   }

   if (missing(roots)) {
       plot.tot.sol <- x$tot.sol
       roots <- 1:plot.tot.sol
   } else {
       if (is.numeric(roots)) {
           plot.tot.sol <- length(roots)
           if (max(roots > tot.sol)) {
               stop("'roots' values must be not greater than the 'x$tot.sol'")
           }
       } else {
           stop("'roots' must be numeric")
       }
   }

   if (missing(which.main)) {
       which.main <- 1:(plot.tot.sol^2)
   } 
   if (!is.numeric(which.main) || any(which.main < 0) || any(which.main > plot.tot.sol^2))
       stop(paste("`which.main' must be in 0:", plot.tot.sol^2, sep=""))
   
   if (plot.tot.sol>1) {
       if (prod(old.par$mfcol)==(length(which.main)+length(roots)*length(which))) {
           par(mfcol=old.par$mfcol)           
       } else if (plot.tot.sol^2==length(which.main)) {
                  par(mfcol=c(plot.tot.sol, plot.tot.sol))
       }
       for (isol in roots) {
            for (jsol in roots) {
                 y.fit.i <- y.fit[isol,]
                 res.i <- res[isol,]
                 weight.i <- weight[isol,]

                 y.fit.j <- y.fit[jsol,]
                 res.j <- res[jsol,]
                 weight.j <- weight[jsol,]
 
                 level.i <- weight.i>=level.weight
                 level.j <- weight.j>=level.weight
 
                 color <- color.res <- color.w <- rep(col[1], size)
                 color[level.i] <- col[2]
                 color[level.j] <- col[3]

                 color.res[level.i & (res.i>res.j)] <- col[2]
                 color.res[level.j & (res.i<res.j)] <- col[3]

                 color.w[level.i & (weight.i>weight.j)] <- col[2]
                 color.w[level.j & (weight.i<weight.j)] <- col[3]


            if (any(which.main==((isol-1)*plot.tot.sol+jsol))) {
                 
                 if (isol==jsol) {
                     ylim <- range(weight.i, na.rm=TRUE)
	             if (id.n > 0)
	                 ylim <- ylim + c(-1,1)* 0.08 * diff(ylim)
                     plot(weight.i, col=color, xlab="Observations", ylab="Weights",main=paste("Weights of the root: ",isol), ylim=ylim)
      	             if (id.n > 0) {
                         show.w <- order(weight.i)[iid] 
	                 x.id <- (1:size)[show.w]
                         y.id <- weight.i[show.w]
#	                 y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
	                 text.id(x=x.id, y=y.id, ind=show.w, labels.id=labels.id,  adj.x = TRUE)
	             }
                 } else {
                     if (isol>jsol) {
                         plot(res.i, res.j, col=color.res, xlab=paste("Residuals of the root: ",isol), ylab=paste("Residuals of the root: ", jsol), main="Residuals")
                         abline(0,1)
                     } else {
                         plot(weight.i, weight.j, col=color.w, xlab=paste("Weights of the root: ", isol), ylab=paste("Weights of the root: ", jsol), main="Weights")
                         abline(0,1)
                     }
                 }
               }        # fine which.main
            }
       }

       for (isol in roots) {
            if (prod(old.par$mfcol)==length(which)) {     
                par(mfcol=old.par$mfcol)
            }
            
            y.fit.temp <- y.fit[isol,]
            res.temp <- res[isol,]
            weight.temp <- weight[isol,]
            level <- weight.temp>=level.weight
            color <- rep(col[1], size)
            color[level] <- col[2]

            if (any(which==1)) {
                ylim <- range(res.temp, na.rm=TRUE)
	        if (id.n > 0)
	            ylim <- ylim + c(-1,1)* 0.08 * diff(ylim)
                plot(y.fit.temp, res.temp, col=color, xlab="Fitted values", ylab="Residuals", ylim=ylim)
      	        if (id.n > 0) {
                    show.w <- order(weight.temp)[iid] 
	            x.id <- y.fit.temp[show.w]
                    y.id <- res.temp[show.w]
                    y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
	            text.id(x=x.id, y=y.id, ind=show.w, labels.id=labels.id, adj.x = TRUE)
	        }
            }

            if (any(which==2)) {
                res.temp.weight <- res.temp*weight.temp
                ylim <- range(res.temp.weight, na.rm=TRUE)
	        if (id.n > 0)
	            ylim <- ylim + c(-1,1)* 0.08 * diff(ylim)            
                plot(y.fit.temp, res.temp.weight, col=color, xlab="Fitted values", ylab="Weighted residuals", ylim=ylim)
      	        if (id.n > 0) {
                    show.w <- order(weight.temp)[iid] 
                    x.id <- y.fit.temp[show.w]
                    y.id <- res.temp.weight[show.w]
                    y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
	            text.id(x=x.id, y=y.id, ind=show.w, labels.id=labels.id, adj.x = TRUE)
	        }
            }

            if (any(which==3)) {
                qqnorm(res.temp, col=color)
                qqline(res.temp)
            }

            if (any(which==4)) {
                qqnorm(res.temp*weight.temp, col=color)
                qqline(res.temp*weight.temp)
            }
       }
   } else {

       if (tot.sol > 1) {
           param <- x$coefficients[roots,]
           res <- x$residuals[roots,]
           y.fit <- x$fitted.values[roots,]
           weight <- x$weights[roots,]
           tot.weight <- x$tot.weights[roots]
       }

       if (prod(old.par$mfcol)!=(length(which)+1)) {
           par(mfcol=c(1,1))
       }

       level.i <- weight>=level.weight
       color <- rep(col[1], size)
       color[level.i] <- col[3]
       show.w <- order(weight)[iid] 

       if (any(which.main!=0)) {            
           ylim <- range(weight, na.rm=TRUE)
           if (id.n > 0)
               ylim <- ylim + c(-1,1)* 0.08 * diff(ylim)
           plot(weight, col=color, xlab="Observations", ylab="Weights", main="Weights of the root", ylim=ylim)
           if (id.n > 0) {
               x.id <- (1:size)[show.w]
               y.id <- weight[show.w]
               y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
	       text.id(x=x.id, y=y.id, ind=show.w, labels.id=labels.id, adj.x = TRUE)
           }
       }
       if (prod(old.par$mfcol)==length(which)) {     
           par(mfcol=old.par$mfcol)
       }
       
       level <- weight>=level.weight
       color <- rep(col[1], size)
       color[level] <- col[2]

       if (any(which==1)) {
           ylim <- range(res, na.rm=TRUE)
           if (id.n > 0)
	       ylim <- ylim + c(-1,1)* 0.08 * diff(ylim)
               plot(y.fit, res, col=color, xlab="Fitted values", ylab="Residuals", ylim=ylim)
           if (id.n > 0) {
               x.id <- y.fit[show.w]
               y.id <- res[show.w]
               y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
	       text.id(x=x.id, y=y.id, ind=show.w, labels.id=labels.id, adj.x = TRUE)
	   }
        }

        if (any(which==2)) {
           res.weight <- res*weight
           ylim <- range(res.weight, na.rm=TRUE)
	   if (id.n > 0)
	       ylim <- ylim + c(-1,1)* 0.08 * diff(ylim)            
           plot(y.fit, res*weight, col=color, xlab="Fitted values", ylab="Weighted residuals", ylim=ylim)
     	   if (id.n > 0) {
               x.id <- y.fit[show.w]
               y.id <- res.weight[show.w]
               y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
	       text.id(x=x.id, y=y.id, ind=show.w, labels.id=labels.id, adj.x = TRUE)
	   }
        }

        if (any(which==3)) {
            qqnorm(res, col=color)
            qqline(res)
        }

       if (any(which==4)) {
            qqnorm(res*weight, col=color)
            qqline(res*weight)
        }
    }
}



