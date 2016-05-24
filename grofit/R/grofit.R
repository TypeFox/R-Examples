grofit <-
function(time, data,  ec50=FALSE, control=grofit.control())
{ 
# /// check input
if (is.numeric(as.matrix(time))==FALSE) stop("Need a numeric matrix for 'time'")
if (is.numeric(as.matrix(data[-1:-3]))==FALSE) stop("Need a numeric matrix for 'data'")

if (is.logical(ec50)==FALSE) stop("Need a logical value for 'ec50'")
if (is(control)!="grofit.control") stop("control must be of class grofit.control!")

nboot.gc  <- control$nboot.gc
nboot.dr  <- control$nboot.dr
out.gcFit <- NA
out.drFit <- NA

# /// fit of growth curves -----------------------------------
out.gcFit <- gcFit(time, data, control)

# /// Estimate EC50 values -----------------------------------	
if (ec50==TRUE){	
	out.drFit <- drFit(summary(out.gcFit), control)	
	# /// rename output and write to harddisk
	EC50.table    <- out.drFit$drTable
	boot.ec       <- out.drFit$boot.ec
}
else{
	ec50.fit  <- NULL
}

grofit        <- list(time=time, data=data, gcFit=out.gcFit, drFit=out.drFit, control=control)
class(grofit) <- "grofit"
grofit
}

