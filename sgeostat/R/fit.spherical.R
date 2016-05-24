"fit.spherical" <-
function (v.object, c0 = 0, cs = 1000, as = 1000, type = "c", 
        iterations = 10, tolerance = 1e-06, echo = FALSE, plot.it = FALSE, 
        weighted = TRUE, delta = 0.1, verbose=TRUE) 
{
        # This program fits a univariate spherical model to an empirical variogram
        # estimate.  The SEMI variogram model is the model fit...
        if (!inherits(v.object, "variogram")) 
                stop("v.object must be of class, \"variogram\".\n")
        if (is.null(c0) & is.null(cs) & is.null(as)) 
                estimate.initial <- TRUE
        else estimate.initial <- FALSE
        # Interpret the "type" argument...
        if (!estimate.initial & (is.null(c0) | is.null(cs) | 
                is.null(as))) 
                stop("c0, cs, and as must all be entered.\n")
        if (type == "c") 
                empgamma <- v.object$classic/2
        else if (type == "r") 
                empgamma <- v.object$robust/2
        else if (type == "m") 
                empgamma <- v.object$med/2
        else stop("type must be \'c\', \'r\', or \'m\'.\n")
        # Set up the variogram function...
        spherical.v <- function(h, parameters) ifelse(h == 0, 
                0, ifelse(h <= parameters[3], parameters[1] + 
                        parameters[2] * (3/2 * h/parameters[3] - 
                                1/2 * (h/parameters[3])^3), parameters[1] + 
                        parameters[2]))
        # Set up the first derivative functions for each of the parameters...
        # dc0 = 1
        # differences, not derivatives
        dcs <- function(h, p, p1) {
                return((spherical.v(h, p) - spherical.v(h, c(p[1], 
                        p1[2], p[3])))/(p[2] - p1[2]))
        }
        das <- function(h, p, p1) {
                return((spherical.v(h, p) - spherical.v(h, c(p[1], 
                        p[2], p1[3])))/(p[3] - p1[3]))
        }
        # Get the number of observations and the bins (h's)...
        numobs <- v.object$n
        h <- v.object$bins
        # If any numobs is 0, get rid of that lag...
        empgamma <- empgamma[numobs > 0]
        h <- h[numobs > 0]
        numobs <- numobs[numobs > 0]
        # Start the sums of squares at 0...
        rse <- 0
        # Begin iterations...
	# initial estimates:
        parameters <- c(c0, cs, as)
        # need two initial parameter estimates to calculate differences
	# second = first + first * delta %)
        if(iterations>0){
          parameters1 <- c(c0, cs, as) * (1 + delta/100)
          cat("Initial parameter estimates: \n first:", parameters, "\n second:", 
              parameters1,"\n")
          loop <- TRUE
          converge <- FALSE
          i <- 1
        } else {
          loop <- FALSE
          converge <- FALSE
          parameters1 <- parameters
        }
        # Plot it before we start if requested...
        if (plot.it) {
                v.m.object <- list(parameters = parameters, model = spherical.v)
                attr(v.m.object, "class") <- "variogram.model"
                attr(v.m.object, "type") <- "spherical"
                plot(v.object, var.mod.obj=v.m.object, type = type)
        }
        while (loop) {
                if(verbose)
                  cat("Iteration:", i, "\n")
                # establish the Y vector...
                y <- (empgamma - spherical.v(h, parameters))
                # establish the x matrix...
                xmat <- cbind(rep(1, length(h)), dcs(h, parameters, 
                        parameters1), das(h, parameters, parameters1))
                # establish the weights (Cressie, p. 99)...
                if (weighted) {
                        w <- numobs/(spherical.v(h, parameters))^2
                }
                else {
                #        w <- (1:length(numobs))
                  w <- rep(1,length(numobs))
                }
                if (echo) 
                        cat("  X matrix:\n")
                if (echo) 
                        print(cbind(y, xmat, w))
                if (echo) 
                        cat("\n\n")
                fit <- lsfit(xmat, y, wt = w, intercept = FALSE)
                # calculate the new parameter estimates...
                parameters.old <- parameters
                parameters <- fit$coef + parameters
                parameters <- ifelse(parameters > 0, parameters, 
                        1e-06)
                if(verbose){
                  cat("Gradient vector: ", fit$coef, "\n")
                  cat("New parameter estimates: ", parameters, 
                      "\n\n")
                }
                # Check for convergence, see if the sum of squares has converged
                rse.old <- rse
                rse <- sum(fit$residuals^2)
                rse.dif <- rse - rse.old
                # Check for convergence of parmeters...
                parm.dist <- sqrt(sum((parameters - parameters.old)^2))
                if(verbose)
                  cat("rse.dif = ", rse.dif, "(rse =", rse, ")  ;  parm.dist = ", 
                      parm.dist, "\n\n")
                #    cat('rse.dif = ',rse.dif,'(rse =',rse,')\n\n')
                #    if(rse.dif < tolerance & parm.dist < tolerance) {
                if (abs(rse.dif) < tolerance) {
                        loop <- FALSE
                        converge <- TRUE
                        cat("Convergence achieved by sums of squares.\n")
                }
                i <- i + 1
                if (i > iterations) {
                        loop <- FALSE
                }
                v.m.object <- list(parameters = parameters, model = spherical.v)
                attr(v.m.object, "class") <- "variogram.model"
                attr(v.m.object, "type") <- "spherical"
                if (plot.it) 
                        plot(v.object, var.mod.obj=v.m.object, type = type)
                parameters1 <- parameters
                parameters <- parameters.old
        }
        if (converge) 
                cat("Final parameter estimates: ", parameters1, 
                        "\n\n")
        else cat("Convergence not achieved!\n")


        v.m.object <- list(parameters = parameters1, model = spherical.v)
        names(v.m.object$parameters)<-c("nugget","sill","range")
        attr(v.m.object, "class") <- "variogram.model"
        attr(v.m.object, "type") <- "spherical"
        return(v.m.object)
}
