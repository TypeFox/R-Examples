# Merged wrapper from Alex's dde and Ben's desolver (RH)
dde <- function( y, times, func, parms=NULL, switchfunc=NULL, 
                 mapfunc=NULL, tol=1e-8, dt=0.1, hbsize=10000 )
{

	if (is.null(hbsize)) hbsize <- 10000;
	if (hbsize<=0) hbsize <- 1;  #if 0 or less, C code crashes
	if (is.null(dt)) dt <- 0.1

	# Check times is valid, unique() and sort() it
	if (!is.null(times) && is.numeric(times)) { # sort out the times vector
		times <- sort(unique(times)); ltimes <- length(times);
		if (ltimes < 2) stop("times must be at least length 2"); }
	else stop("valid times vector must be given");

	# Create (from,to,by) for C-code
	from <- times[1];  to <- max(times);

	out <- .Call("startDDE", 
	             gradFunc=func,
	             switchFunc=switchfunc,
	             mapFunc=mapfunc,
	             env=new.env(), 
	             yinit=y,
	             parms=parms,
	             settings=c(tol, from, to, dt, hbsize),
	             outtimes=as.numeric(times), # bjc 2007-05-08
	             PACKAGE = "PBSddesolve")

	# If .Call("startDDE", ...) fails, return a dummy matrix w/ 2 columns
	tout <- try(data.frame(out));
	if (class(tout)=="try-error") 
		tout <- matrix(data=c(times, NA*times), nrow=ltimes, ncol=2) 
	#tout <- as.matrix(out);
	attr(tout, "istate") <- 2; # rh: what's this?
	return(tout) 
}

# t - at what time
# markno - used for optimization
pastvalue <- function(t)
{
	markno=0 #used for optimization when more than one lag time is used
	         #but will need code to update `data.nlag' in me95.c:setupglobaldata
    out <- .Call("getPastValue", 
                 t=t, 
                 markno=as.integer(markno),
                 PACKAGE = "PBSddesolve")
	return(out)
}

# t - at what time
# markno - used for optimization
pastgradient <- function(t)
{
	markno=0 #used for optimization when more than one lag time is used
	         #but will need code to update `data.nlag' in me95.c:setupglobaldata
    out <- .Call("getPastGradient", 
                 t=t, 
                 markno=as.integer(markno),
                 PACKAGE = "PBSddesolve")
	return(out)
}

.onClosePBSddeExamples = function(){
	setwd(.PBSddeEnv$oldwd) }


