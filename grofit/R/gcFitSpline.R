gcFitSpline <-
function(time,data, gcID ="undefined", control=grofit.control())
{

# /// check input parameters
if (is(control)!="grofit.control") stop("control must be of class grofit.control!")
if (!control$fit.opt%in%c("s","b")) stop("Fit option is not set for a spline fit. See grofit.control()")

# /// conversion to handle even data.frame inputs
time <- as.vector(as.numeric(as.matrix(time)))
data <- as.vector(as.numeric(as.matrix(data)))

# /// check length of input data
if (length(time)!=length(data)) stop("gcFitSpline: length of input vectors differ!")

# /// determine which values are not valid
bad.values <-  (is.na(time))|(time<0)|(is.na(data))|(data<0)|(!is.numeric(time))|(!is.numeric(data))

# /// remove bad values or stop program
if (TRUE%in%bad.values)
{
   if (control$neg.nan.act==FALSE)
   {
      time    <- time[!bad.values]
      data    <- data[!bad.values]
   }
   else{
   stop("Bad values in gcFitSpline")
   }
}

if (length(data)<5){
	cat("gcFitSpline: There is not enough valid data. Must have at least 5!")
        gcFitSpline <- list(raw.time = time, raw.data = data, gcID = gcID, fit.time = NA, fit.data = NA, parameters = list(A= NA, mu=NA, lambda=NA, integral=NA), parametersLowess=list(A=NA, mu=NA, lambda=NA), spline = NA, reliable=NULL, fitFlag=FALSE, control = control)
        class(gcFitSpline) <- "gcFitSpline"
        return(gcFitSpline)
}
else
{
	# /// apply transformation
	if (control$log.x.gc==TRUE){time <- log(1+time)}
	if (control$log.y.gc==TRUE){ data   <- log(1+data)}

	#will be used as start value in nls
	halftime <- (min(time)+max(time))/2

	# spline fit and computation of the maximum derivative
	try(y.spl <- smooth.spline(time,data,spar=control$smooth.gc))
	if (is.null(y.spl)==TRUE){
		warning("Spline could not be fitted to data!")
		if (is.null(control$smooth.gc)==TRUE){
			cat("This might be caused by usage of smoothing parameter NULL\n")
			fit.nonpara        <- list(raw.x = time, raw.y = data, fit.x = NA, fit.y = NA, parameters = list(A= NA, mu=NA, lambda=NA, integral=NA), spline = NA, parametersLowess=list(A= NA, mu=NA, lambda=NA), spline = NA, reliable=NULL, fitFlag=FALSE, control = control)
			class(gcFitSpline) <- "gcFitSpline"
			return(gcFitSpline)
		}		    	
	}
	# spline fit
	dydt.spl   <- predict(y.spl, time, deriv = 1)
	index      <- which.max(dydt.spl$y)          #index of maximum derivative
	t.max      <- dydt.spl$x[index]
	dydt.max   <- max(dydt.spl$y)
	y.max      <- y.spl$y[index]
	mu.spl     <- dydt.max;
	b.spl      <- y.max-dydt.max*t.max           #intercept
	lambda.spl <- -b.spl/mu.spl
	integral   <- low.integrate(y.spl$x,y.spl$y) #the integral under the curve
	
	# lowess fit
	low        <- lowess(time,data,f=0.25)
	y.low      <- low$y
	x.low      <- low$x
	dydt.low   <- diff(y.low)/diff(time)
	mu.low     <- max(dydt.low)
	index      <- which.max(dydt.low)            #index of maximum derivative
	t.max      <- x.low[index]
	y.max      <- y.low[index]
	b.low      <- y.max-mu.low*t.max             #intercept
	lambda.low <- (-1)*b.low/mu.low
	

}

gcFitSpline        <- list(raw.time = time, raw.data = data, gcID = gcID, fit.time = y.spl$x, fit.data = y.spl$y, parameters = list(A= max(y.spl$y), mu=mu.spl, lambda=lambda.spl, integral=integral), parametersLowess=list(A= max(y.low), mu=mu.low, lambda=lambda.low), spline = y.spl, reliable=NULL, fitFlag=TRUE, control = control)

class(gcFitSpline) <- "gcFitSpline"
gcFitSpline

}

