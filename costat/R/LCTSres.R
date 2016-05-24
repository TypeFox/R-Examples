LCTSres <-
function (res, tsx, tsy, inc = 0, solno = 1:nrow(res$endpar), 
    filter.number = 1, family = c("DaubExPhase", "DaubLeAsymm"),
    plot.it = FALSE, spec.filter.number = 1,
    spec.family = c("DaubExPhase", "DaubLeAsymm"), plotcoef = FALSE, 
    sameplot = TRUE, norm = FALSE, plotstystat = FALSE, plotsolinfo = TRUE, 
    onlyacfs = FALSE, acfdatatrans = I, xlab = "Time", ...) 
{
family <- match.arg(family)
spec.family <- match.arg(spec.family)
#
# Produces plots of solutions identified in csFSS class object.
# This function is interactive and can look at more than one solution per call 
#
#
# Check object is of required class
#
if (class(res) != "csFSS") 
        stop("res object must be of class `csFSS`")
#
# Work out some dimensions of the object
#
nopt <- nrow(res$endpar)
lcfs <- ncol(res$endpar)
lts <- length(tsx)
#
# Add an increment onto the time axis
#
xx <- inc + 1:lts
#
# Communicate with the user
#
cat("The solutions you wish to look at are:")
print(solno)
#
# Produce plots for each requested solution
#
for (i in 1:nopt) {
	if (!is.na(match(i, solno))) {
		cat("Solution Number ", i, "\n")
		#
		# Only plot for converged solutions
		3
		if (res$convergence[i] == 0) {
			#
			# Get coefficients associated with this soln,
			# breakout into separate vectors and turn into functions
			#
			cfs <- res$endpar[i, ]
			alpha <- cfs[1:(lcfs/2)]
			beta <- cfs[(lcfs/2 + 1):lcfs]
			v <- coeftofn(alpha = alpha, beta = beta, n = lts, 
				filter.number = filter.number, family = family)
			#
			# Form time-varying combination Z_t, and compute its
			# spectrum
			#	
			lcts <- v$alpha * tsx + v$beta * tsy
			lctsspec <- ewspec(lcts,
				filter.number = spec.filter.number, 
				family = spec.family)$S

			#
			# Produce various plots depending on arguments
			#
			if (plotcoef == TRUE) {
			    if (sameplot == TRUE) {
				if (norm == TRUE) 
				    v$alpha <- v$alpha/sqrt(sum(v$alpha^2))
				    if (i == 1) {
					plot(xx, v$alpha, main = "alpha",
						type = "l", xlab = xlab,
						ylab = "alpha_t")
					}
				    else {
				        lines(xx, v$alpha)
				        }
				    }
			    else {
				plot(xx, v$alpha, main = "alpha", type = "l", 
				    xlab = xlab)
			    }
                	}
			else {
			    if (plotsolinfo == TRUE) {
                    		plot(xx, v$alpha, main = "alpha", type = "l", 
				    xlab = xlab, ylab = "alpha_t")
				plot(xx, v$beta, main = "beta", type = "l", 
				    xlab = xlab, ylab = "beta_t")
				plot(xx, lcts, xlab = xlab, type = "l",
				    main = paste("Combined. Minvar: ", 
				    signif(res$minvar[i], 3)), ylab = "Z_t")
			        plot(lctsspec, main =
				    paste("p-val is", signif(res$pvals[i], 3)),
				    xlab = xlab, sub = "", ylab = "Scale j")
				}
			    scan()
			    if (plotstystat == TRUE) {
				if (!onlyacfs) 
				    ts.plot(lcts, main = "Time series")
				acf(acfdatatrans(lcts), ...)
				acf(acfdatatrans(lcts), type = "partial", ...)
			        if (!onlyacfs) 
				    spectrum(lcts, span = c(5, 7))
				scan()
				}
			    }
			}
		}
	}
return(lcts)
}
