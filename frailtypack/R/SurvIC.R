"SurvIC" <- function(t0,lower,upper,event) {

# nombre d'arguments manquants
ng <- missing(t0) + missing(lower) + missing(upper) + missing(event)

# si un seul argument est fourni, c'est que c'est un lower
if (ng==3) lower <- t0
# si deux arguments fournis, c'est lower et upper
if (ng==2){
	upper <- lower
	lower <- t0
}
# si trois arguments fournis, lower upper event
if (ng==1){
	event <- upper
	upper <- lower
	lower <- t0
}
if (missing(lower)) stop("Must have a lower time argument")
if (!is.numeric(lower)) stop("Lower time variable is not numeric")
if (any(is.na(lower))) stop("There is some NA in your lower time")

if (missing(upper)) stop("Must have an upper time argument")
if (!is.numeric(upper)) stop("Upper time variable is not numeric")
if (any(is.na(upper))) stop("There is some NA in your upper time")

if (missing(event)) stop("Must have an event argument")
if (!is.numeric(event)) stop("Event variable is not numeric")
if (any(is.na(event))) stop("There is some NA in your event argument")

if (length(lower)!=length(upper)) stop("Lower and upper time are different lengths")
if (any(lower>upper)) stop("Lower time has to be less than upper time")

if (length(lower==upper) != length(event==0)) warning("There may be an error in the right censored times")

if ((ng==0) & (any(t0>lower))) warning("Be careful, some of your truncation times are higher than the lower")

# s'il manque le t0
if (ng==1){
	ss <- cbind(lower=lower,upper=upper,status=event)
	attr(ss,"type") <- "interval"
}else{
	ss <- cbind(t0=t0,lower=lower,upper=upper,status=event)
	attr(ss,"type") <- "intervaltronc"
}

if (is.R()) { class(ss) <- "SurvIC" }
else { oldClass(ss) <- "SurvIC" }
return(ss)
}