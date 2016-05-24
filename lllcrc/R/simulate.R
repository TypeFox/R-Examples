#' Simulate a CRC experiment
#' 
#' Generate fake capture-recapture datasets to be used in simulation
#' experiments
#' 
#' 
#' @param k The number of lists to be generated
#' @param n The true size of the [fake] population
#' @param heters Optional: a string that is the name of some function for
#' defining a heterogeneity structure, possibly in terms of covariates.
#' Currently the only non-null option implemented is \code{\link{french.1}}, although
#' a user could easily create additional options
#' @param covs Optional: a string that is the name of some function for
#' defining covariates for each population unit.  Currently the only non-null
#' option implemented is \code{age.sex.zip}
#' @param remove.zeros Logical, TRUE by default.  If FALSE, keep the all-zero
#' capture patterns in the dataset to be returned
#' @return A data frame containing fake capture recapture data
#' @author Zach Kurtz
#' @export poptop
poptop = function(k=2, n=100, heters = NULL, covs = NULL, remove.zeros = TRUE)
{ #k = 3; n=100; heters = "french.1"; covs = "age.sex.zip"
	## Envision a population prior to any CRC business
	pop = init.pop(n, covs)
	## Perform a CRC experiment that may depend on covariates and heterogeneity
	pop = captures(n, pop, k, heters, covs)
	## Get the capture patterns
	pop = patterns(pop, k)
	if(remove.zeros) pop = pop[pop$y != paste(rep(0,k), collapse = ""),]
	for(i in 1:k) pop[, paste("c", i, sep = "")] = NULL
	return(pop)
}



#' Set up a fake population
#' 
#' Simply initialized a data frame to represent the CRC data in a simulation
#' experiment.
#' 
#' This is just a workhorse function, typically called by \code{poptop}.
#' 
#' @param n The population size.
#' @param covs A string specifying a function for generating covariate values.
#' By default, there are no covariates.  The only non-default method at present
#' is \code{french.1}.
#' @return A data.frame representing a population, possibly with covariates.
#' If there are no covariates, it is just a single column labeled `id'.
#' @author Zach Kurtz
#' @export init.pop
init.pop = function(n, covs = NULL)
{
	# With no covariates, the population is simply 1:n
	if(is.null(covs)){
		pop = data.frame(id = 1:n)   
	# With covariates, the population is a collection of covariate values
		}else if (covs == "age.sex.zip"){
		pop = age.sex.zip(n)
		}else{stop("covs method invalid")
	}
	return(pop)
}



#' Simulating captures
#' 
#' Simulating captures as part of a CRC experiment.
#' 
#' Given (possibly) some covariates specified by \code{covs} and a
#' heterogeneity function specified by \code{heters}, conduct a
#' capture-recapture experiment, and return the capture-recapture data.
#' 
#' @param n The population size.
#' @param pop The data.frame with a basic template for the population.
#' @param k The number of lists.
#' @param heters A character string to denote a specific heterogeneity
#' structure.  Currently only french.1 is implemented.
#' @param covs A character string to denote the covariates to be included.
#' Currently only NULL and age.sex.zip are implemented.
#' @return A data.frame with capture-recapture data
#' @author Zach Kurtz
#' @export captures
captures = function(n, pop, k, heters, covs)
{
	## Initialize capture column headings
	c.names = apply(cbind(rep("c", k), 1:k), MARGIN = 1, FUN = paste, collapse = "")
	for(i in 1:k) pop[,c.names[i]] = rep(NA, n)
	## Perform CRC experiment
	if(is.null(heters)) {  
		for(i in 1:k) pop[, paste("c", i, sep = "")] = rbinom(n, 1, 0.5)
	}else if (heters == "french.1"){pop = french.1(pop, n, k)
	}else{print("Error from fun=captures. Heters method invalid")}
	return(pop)
}



#' Collapse capture events into capture patterns (strings)
#' 
#' For any capture-recapture data matrix with columns such as \code{c1},
#' \code{c2} ... as binary indicators of capture or noncapture, running
#' \code{patterns} on the data appends a new column \code{y} that binds
#' \code{c1}, \code{c2}, ... together as a column of capture patterns.
#' 
#' 
#' @param pop Capture-recapture data, at an intermediate form in the inner
#' workings of \code{poptop}
#' @param k The number of lists
#' @author Zach Kurtz
#' @export patterns
patterns = function(pop, k)
{
	c.names = names(pop)[substring(names(pop),1,1) == "c"]
	pop.caps = pop[, c.names]
	pop$y = apply(pop.caps, MARGIN = 1, FUN = paste, collapse = "")
	return(pop)
}



#' Simulate CRC data with age, sex, and zip code
#' 
#' Generate a fake data set with covariates age, sex, and zip code.
#' 
#' This simulation is vaguely inspired by the French dataset on multiple
#' sclerosis prevalence (El Adssi et. al., 2012).  It is used to create a
#' three-list capture-recapture system with unit-level independence between
#' lists.  Capture probabilities are constant (for simplicity!) across lists,
#' but not across population units.  The heterogeneity is a function of the
#' covariates: age, sex, and zip code.  We code zip here as `a', `b', `c', `d'.
#' 
#' @param n The number of population units.
#' @return \item{pop}{ A data.frame containing the covariates and capture
#' histories associated with each observed unit. }
#' @author Zach Kurtz
#' @export age.sex.zip
age.sex.zip = function(n)
{
	#age
	age = round(rnorm(n, mean = 50, sd = 15),2)
	extreme.ages = which((age<15)|(age>90))
	age[extreme.ages] = round(runif(length(extreme.ages), min = 15, max = 90),2)	
	#sex	
	sex = sample(size = n, c("M", "F"), prob = c(1/3, 2/3), replace = TRUE)
	#zip
	zip = sample(size = n, c("a", "b", "c", "d"), prob = c(0.1, 0.2, 0.3, 0.4), replace = TRUE)
	out = data.frame(age=age, sex=sex, zip=zip)
	colnames(out) = c("age", "sex", "zip")	
	return(out)
}



#' A fake dataset, french.1
#' 
#' Generate a fake dataset, for simulation purposes, that is somewhat similar
#' to the data on multiple sclerosis in the Lorraine Region of France as in El
#' Adssi et. al. (2012).
#' 
#' This function is the guts for a capture-recapture experiment that is
#' specified by the \code{heters} argument in the function \code{poptop}; with
#' a minor revision ot the function \code{captures}, you could create your own
#' guts function for a very different capture-recapture experiment that sits
#' within the same coding framework.
#' 
#' @param pop A data frame containing a template for the population.
#' @param n The size of the target population.
#' @param k The number of lists.
#' @return A data frame, which gets passed back up to the main function
#' \code{poptop} for a bit of cleaning.
#' @author Zach Kurtz
#' @references
#' El Adssi H, Debouverie M and Guillemin F (2012). "Estimating the
#' prevalence and incidence of multiple sclerosis in the Lorraine region,
#' France, by the capture-Recapture method." \emph{Multiple Sclerosis Journal},
#' \bold{18}(9), pp. 1244-1250.
#' @export french.1
french.1 = function(pop, n, k)
{
	# Define heterogeneous capture probabilities in ways
	#	that depend on the covariates
	p = rep(0.21,n)
	m.sex.M = 1.1; ind.M = which(pop$sex == "M")
	m.zip.b = 1.1; ind.b = which(pop$zip == "b")
	m.zip.c = 1.2; ind.c = which(pop$zip == "c")
	m.zip.d = 1.3; ind.d = which(pop$zip == "d")
	m.age = 3 - 3*(abs(pop$age - 50)/45)
	p[ind.M] = p[ind.M]*m.sex.M
	p[ind.b] = p[ind.b]*m.zip.b
	p[ind.c] = p[ind.c]*m.zip.c
	p[ind.d] = p[ind.d]*m.zip.d
	p = p*m.age
	for(i in 1:k) pop[,paste("c",i, sep = "")] = rbinom(n,1,p)
	return(pop)
}
