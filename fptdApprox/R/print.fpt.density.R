print.fpt.density <-
function (x, ...) 
{
    	if (!is.fpt.density(x)) 
        	stop(paste(sQuote("x"), "is not of class", shQuote("fpt.density")))

	cp <- is.null(x$y.x0)
    	A <- attr(x, "Call")
	f <- A[[1]]
	Args <- formals(eval(f))
	logic <- (as.character(f) == "Approx.fpt.density")
    	Args[names(A)[-1]] <- A[-1]

	t0 <- attr(attr(x, "summary.fptl"), "FPTLCall")[[1]]$t0
	T <- attr(attr(x, "summary.fptl"), "FPTLCall")[[1]]$T

    	cutoff <- options()$deparse.cutoff
    	options(deparse.cutoff = 175)
  
    	cat("\nAn object of class", shQuote("fpt.density"), "containing")
    	cat(paste("\n   $x: a sequence of", length(x$x), "time instants from", ifelse(t0 < x$x[1], format(x$x[1], ...), t0), "to", 
		ifelse(x$x[length(x$x)] < T, format(x$x[length(x$x)], ...), T), sep = " "))	
    	cat("\n   $y: the values of the approximate first-passage-time density function on sequence x")

    	cat("\n\nCall:")
    	cat("\n")
    	print(attr(x, "Call"))

    	cat("\nSpecified options to apply the numerical algorithm:")

    	label <- c("Variable integration step", "Calculate the approximation from the lower end of the interval considered",
				"Calculate the approximation to the upper end of the interval considered", 
				"Skip the subintervals at which the FPTL decreases and the density of f.p.t. is almost zero",
				"Stop the algorithm if the cumulative integral at potential final time instants is over")				
	if (!cp) label <- c(label, "Number of equally spaced values considered in the range of the initial distribution") 
	if (logic) label <- c(label, "Maximum slope required between two points to consider that a growing function is constant", 
					"Ratio of the global increase of the FPTL function in the growth subinterval [t[i], tmax[i]]", 
					"that should be reached to consider that it begins to increase significantly",
					"Maximum allowable distance between tmax[i] and tmax[i]^+")
	label <- c(label, "Number of points per amplitude unit used to determine optimal integration steps in [t[i]*, tmax[i]^+]",
			"Number of points per amplitude unit used to determine optimal integration step in [tmax[i]^+, t[i+1]*]", 
			"when (t[i+1]* - tmax[i]^+) <= (tmax[i]^- - t[i]*)",
			"Number of points per amplitude unit used to determine optimal integration step in [tmax[i]^+, t[i+1]*]", 
			"when (t[i+1]* - tmax[i]^+) > (tmax[i]^- - t[i]*)") 
    	label <- paste("   ", format(label), "   ")
	values <- as.character(c(Args[c("variableStep", "from.t0", "to.T", "skip")], 1 - Args$tol))
	if (!cp) values <- c(values, as.character(Args$m))
	if (logic) values <- c(values, as.character(Args$zeroSlope), paste("10^(-", Args$p0.tol, ")", sep = ""), "", paste(Args$k, "*(tmax[i] - t[i]*)(1-FPTL(tmax[i]))", sep = ""))
	values <- c(values, paste(Args$n, "/(tmax[i]^- - t[i]*)", sep = ""), paste(Args$p, "*", Args$n, "/(tmax[i]^- - t[i]*)", sep = ""), "",
			ifelse(Args$alpha == 1L, paste(Args$p, "*", Args$n, "/(tmax[i]^- - t[i]*)", sep = ""), paste(Args$p, "*", Args$n, 
			"*(t[i+1]* - tmax[i]^+)^(", Args$alpha - 1L, ")/(tmax[i]^- - t[i]*)^", Args$alpha, sep = "")), "")
	label <- paste(label, values)
    	cat("\n\n")
    	cat(label, sep= "\n")

    	nI <- length(attr(x, "cumIntegral"))
    	index <- 1:nI
    	y <- data.frame(matrix(, nrow = nI, ncol = 5))
    	names(y) <- c("Subinterval", "Integration step", "Cumulative integral", "Iterations", "User time")
    
    	lower <- attr(x, "Steps")[index, 1]
    	upper <- attr(x, "Steps")[index, 2]

	if (cp) m <- 1L
	else{
		m <- Args$m
	      if (is.null(m)) m <- 100L
	}
    	jumps <- which(sapply(attr(x, "skips"), identical, 1:m))

    	y[, 2] <- attr(x, "Steps")[index, 3]
    	y[, 3] <- attr(x, "cumIntegral")
    	y[, 4] <- (upper - lower)/y[, 2]
    	y[, 5] <- attr(x, "CPUTime")[,1]

	if (length(jumps) > 0L) y[jumps, 4] <- 0L

    	it <- format(sum(y[setdiff(index, jumps), 4]), ...) 
    	ut <- format(sum(y[, 5]), ...)

    	y <- format(y, ...)
    	y[,1] <- paste("(", format(lower, ...), ", ", format(upper, ...), "]", sep="")      					                    		  
    
    	cat("\nApplication summary of the numerical algorithm:\n\n")       	
    	y <- rbind(c("Subinterval", "Integration step", "Cumulative integral", "Iterations", "User time"), y)
	y <- apply(y, 2, format, justify = "right")
	cat(apply(y, 1, paste, collapse = "  "), sep = "\n")
    	label <- paste(format(c("Total number of iterations", "Total user time")), "   ", c(it, ut))
    	cat("\n")
    	cat(label, sep= "\n")

    	options(deparse.cutoff = cutoff)
}
