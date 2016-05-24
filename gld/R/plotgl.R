plotgld <- function(lambda1=0,lambda2=NULL,lambda3=NULL,lambda4=NULL,
  param="fmkl",lambda5=NULL, add=NULL, truncate = 0, bnw = FALSE,
  col.or.type = 1, granularity = 10000,xlab = "x", ylab=NULL,
  quant.probs = seq(0,1,.25), new.plot = NULL, ...)
{
# standard parameter fixin - copied directly from dgl, but we want the 
# warnings to happen in this function.
# Tidy the parameters so gl.check.lambda will work
lambdas <- .gl.parameter.tidy(lambda1,lambda2,lambda3,lambda4,param,lambda5)
# Check the parameters
if(!gl.check.lambda(lambdas,param=param,vect=TRUE)) {
        stop(paste("The parameter values", lambdas,
"\ndo not produce a proper distribution with the",param,
"parameterisation - see \ndocumentation for gl.check.lambda"))
        }
if (is.null(new.plot)) {
	if (is.null(add)) {new.plot <- TRUE } else {
	new.plot <- !add } } else { 
	if (!is.null(add)) {warning("new.plot and add both set.  Using add and ignoring new.plot") 
		new.plot <- !add
}} 
	# If neither set - then add should be FALSE, new.plot TRUE
	# If add is set, set new.plot to be used internally
	# If new.plot is set, use it internally
	# If both are set, give a warning and use add.
u <- seq(from = 0, to = 1, by = 1/granularity)
# Only difference across parameterisations is calculating the 
# quantiles and density
quantiles <- qgl(u,lambda1=lambdas,param=param)
density <- dqgl(u,lambda1=lambdas,param=param)
# Check for pathological densities - include check for all,
dots.required <- FALSE
nearzero <- .Machine$double.eps*c(0,1,3)
check.for.jump <- qgl(nearzero,lambda1=lambdas,param=param)
check.for.jump2 <- qgl(1-nearzero,lambda1=lambdas,param=param)
if (is.finite(check.for.jump[1])) # warning and different plot
	{
	first.space <- check.for.jump[2] - check.for.jump[1]
	from2to3 <- check.for.jump[3] - check.for.jump[2]
	if ((first.space > (from2to3 * 20)) & (from2to3>0)) {
		warning(paste("These parameter values give a pathological density.  \n",
		"The theoretical minimum: F^{-1}(0)= ",signif(check.for.jump[1],4),"\n is much less than F^{-1}("
		,signif(nearzero[2],4),")= ",signif(check.for.jump[2],4),".\n","The density is undefined between these points",sep=""))
		dots.required <- TRUE
		dots.start <- TRUE
		dots.end <- FALSE
		}
	}
if (is.finite(check.for.jump2[1])) # warning and different plot
	{
	last.space <- check.for.jump2[1] - check.for.jump2[2]
	from2to3 <- check.for.jump2[2] - check.for.jump2[3]
	if ((last.space > (from2to3 * 20)) & (from2to3>0)){
		warning(paste("These parameter values give a pathological density.  \n",
		"The theoretical maximum: F^{-1}(1)= ",signif(check.for.jump2[1],4),
		"\n is much more than F^{-1}(1-", signif(nearzero[2],4),")= ",signif(check.for.jump2[2],4),".\n","The density is undefined between these points",sep=""))
		dots.required <- TRUE
		dots.end <- TRUE
		if (!exists("dots.start")) { dots.start <- FALSE }
		}
	}
if(truncate > 0) { 
	# If truncated, not RS pathological problem, because we are not plotting density at the endpoint
	if(new.plot) {
		if (is.null(ylab)){
			ylab <- paste( "probability density (values below", deparse(substitute(truncate)), "not shown)")
			}
		if(bnw) {
			plot(quantiles[density > truncate], 
			density[density > truncate], type = "l", 
			xlab = xlab, ylab = ylab, lty = col.or.type,...)
			}
		else {
			plot(quantiles[density > truncate], 
			density[density > truncate], type = "l", 
			xlab = xlab, ylab = ylab, col = col.or.type,...)
			}
		}
	else {
		if(bnw) {
			lines(quantiles[density > truncate], density[
			  density > truncate], lty = col.or.type)
			}
		else {
			lines(quantiles[density > truncate], density[
			  density > truncate], col = col.or.type)
			}
		}
	}
else {
	if(new.plot) {
		if (is.null(ylab)){ ylab <- "probability density" }
		# Set up graphics system for plot with dots
		if (dots.required) {
			plot(quantiles, density, type = "n", xlab = xlab,ylab = ylab, lty=col.or.type, ...)
			if (dots.start) {
				if(bnw) {points(quantiles[1],density[1]) }
				else {points(quantiles[1],density[1],col=col.or.type) } 
				start.line <- 2
				}
			else {start.line <- 1}
			if (dots.end) {
				if (bnw) {points(quantiles[granularity+1],density[granularity+1]) }
				else {points(quantiles[granularity+1],density[granularity+1],col=col.or.type) }
				end.line <- granularity
				}
			else {end.line <- granularity+1}
			}
		else {
			start.line <- 1
			end.line <- granularity+1
			}
		if(bnw) {
			if (dots.required) { lines(quantiles[start.line:end.line], density[start.line:end.line], lty=col.or.type, ...) }
			else { plot(quantiles[start.line:end.line], density[start.line:end.line], type = "l", xlab = xlab,ylab = ylab, lty=col.or.type, ...) }
			}
		else {	
			if (dots.required) { lines(quantiles[start.line:end.line], density[start.line:end.line], col=col.or.type, ...) }
			else { plot(quantiles[start.line:end.line], density[start.line:end.line], type = "l", xlab = xlab,ylab = ylab, col=col.or.type, ...) }
			}
		}
	else {
		# Not a new plot - so the initial plots with type "n" arent needed in the
		# dots required case
		if (dots.required) {
			if (dots.start) {
	                        if(bnw) { points(quantiles[1],density[1]) } else {
					points(quantiles[1],density[1],col=col.or.type) } 
                                start.line <- 2 } else { start.line <- 1 
				}
			if (dots.end) {
                                if (bnw) {points(quantiles[granularity+1],density[granularity+1]) }
                                else {points(quantiles[granularity+1],density[granularity+1],col=col.or.type) }
                                end.line <- granularity
                                } else {end.line <- granularity+1}
			} else {
			start.line <- 1
                        end.line <- granularity+1
		}
		if (bnw) {lines(quantiles[start.line:end.line], density[start.line:end.line], lty = col.or.type)}
		else {lines(quantiles[start.line:end.line], density[start.line:end.line], type="l",xlab = xlab,ylab = ylab, col=col.or.type, ...) }
		}
	}
if (!is.null(quant.probs)){quantile(quantiles,quant.probs) } 
}

plotglc <- function(lambda1=0,lambda2=NULL,lambda3=NULL,lambda4=NULL,
  param="fmkl",lambda5=NULL, granularity=10000, xlab="x",
  ylab="cumulative probability",add=FALSE,...)
{
# standard parameter fixin' - copied directly from dgl, but we want the 
# warnings to happen in this function.
# Tidy the parameters so gl.check.lambda will work
lambdas <- .gl.parameter.tidy(lambda1,lambda2,lambda3,lambda4,param,lambda5)
# Check the parameters
if(!gl.check.lambda(lambdas,param=param,vect=TRUE)) {
        stop(paste("The parameter values", lambdas,
"\ndo not produce a proper distribution with the",param,
"parameterisation - see \ndocumentation for gl.check.lambda"))
        }
	u <- seq(from = 1/granularity, to = 1 - 1/granularity, length = 
		granularity - 1)
	x <- qgl(u,lambda1=lambdas,param=param)
	if (add) {lines(x, u, pch = ".",xlab=xlab,ylab=ylab,...)} else {
	plot(x, u, pch = ".",xlab=xlab,ylab=ylab,...) }
}
