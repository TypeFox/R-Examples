#The version that gets stored is the one below
#it has runeff = e1; dayeff = e2; unitvar = e3
"strength.gen" <- function(nodays, noruns, temperature, ident)
{
#
# set up numeric vectors and generate strengths
#
	days <- rep(1:nodays, rep(noruns, nodays))
	runs <- rep(1:noruns, nodays)
	runeff <- (1:noruns) * 6 - 18
	strengths <- rep(runeff, nodays)
	dayeff <- sin(3.141593 * ((1:nodays) - 1)/8) * 15 - 5
	id <- ident %% 1024
	set.seed(id)
	nounits <- nodays * noruns
  strengths <- strengths + 25 + rep(dayeff, rep(noruns, nodays)) + rnorm(
		nounits, 0, 1)
	temp <- as.vector(as.character(temperature), mode = "numeric")
  strengths <- strengths + 0.25/(0.00584 * temp - 1.082 + 52.4/temp)	#
# add unit variability 
#
	set.seed(541)
  unitvar <- rnorm(nounits, 0, 3)
	day <- factor(days)
	run <- factor(runs)
	aov.unitvar <- aov(unitvar ~ day + run + temperature)
	summary(aov.unitvar)
	res <- residuals(aov.unitvar)
	strengths <- floor(strengths + res)  #/sqrt(var(res)) * 2)	#
# return the factors and generated strengths in a data frame
#
	exp.strength <- data.frame(day, run, temperature, strengths)
}

"strength" <- function(nodays, noruns, temperature, ident)
{
#test arguments
#
	if(!is.factor(temperature)) stop("Temperature must be a factor")
	nounits <- nodays * noruns
	if(length(temperature) != nounits)
		stop("Length of temperature must be nodays * noruns")
	temp <- as.vector(as.character(temperature), mode = "numeric")
	if(min(temp) < 75 || max(temp) > 125) stop(
			"Temperatures outside operating range of 75-125 C specified"
			)	
	if (noruns >5) stop("The number of runs per day must be 5 or less")
	if (nounits > 35) stop ("The maximum total number of runs is 35") 
#
# generate strengths and return the factors and generated strengths in a data frame
#
	exp.strength <- strength.gen(nodays, noruns, temperature, ident)
}
