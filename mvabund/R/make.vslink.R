#####################################################################
## make.vslink. a function to create a glm-link object             ##
## for the variance stabilising link                               ##
## this function works like make link, it is not expected 	       ##
## to be used directly by the user	                               ##
#####################################################################
## 
## Arguments:
## family	 the exponential family used, allowed are gaussian,
##			   poisson, quasipoisson, binomial # quasibinomial
## link 	 character "varstab",
##			   other links are not recognised by this function,
##			   use 'make.link' instead
## 
## returned Value: 
## An object of class "link-glm", a list with components 
## 
## linkfun 	Link function function(mu) 
## linkinv 	Inverse link function function(eta) # eta = x beta
## mu.eta 	Derivative function(eta) dmu/deta 
## valideta function(eta){ TRUE if eta is in the domain of linkinv }. 
## name 	  the name "varstab" that is used for the link,
##          if applicable the name has as attribute the name of the
##          identical link for this family
#####################################################################
## Examples for usage:
## poisson(make.vslink(poisson))
## or
## linkfct <- make.vslink(poisson)
## # note that the link name must be different from "link"
## poisson(linkfct)
####################################################################

make.vslink <- function(family, link="varstab", theta)           # , theta=NULL
{

	if(link!="varstab") stop(sQuote(link), " link not recognised")
	fam <- substitute(family)
	if(is.character(fam)) family <- get(fam)

	if(is.function(family))
			family <- family()

	if(inherits(family,"family")) {
			fam <- family$family
	} else if (is.character(fam)) stop("'family' is invalid")	

	if(is.null(fam)) stop("'family' is invalid")
		
	if(fam=="gaussian") {
    ml <- make.link("identity")
    attrml <- ml$name
    ml$name <- "varstab"
    attr(ml$name, "name") <- attrml
    attr(ml$name, "expression") <- expression(make.link("identity"))
		return(ml)
	} else if (fam=="poisson" | fam=="quasipoisson"){ 
	  ml <- make.link("sqrt")
    attrml <- ml$name
    ml$name <- "varstab"
    attr(ml$name, "name") <- attrml
    if (fam=="poisson"){
      attr(ml$name, "expression") <- 
		expression(make.vslink(family="poisson"))
		return(ml)
    } else {
      attr(ml$name, "expression") <- 
		expression(make.vslink(family="quasipoisson"))
		return(ml)
		}
	} else if (fam=="Gamma") {
    ml <- make.link("log")
    attrml <- ml$name
    ml$name <- "varstab"
    attr(ml$name, "name") <- attrml
    attr(ml$name, "expression") <- 
	expression(make.vslink(family="Gamma"))
		return(ml)
	} else if (fam=="inverse.gaussian") {
    ml <- make.link("1/mu^2")
    attrml <- ml$name
    ml$name <- "varstab"
    attr(ml$name, "name") <- attrml
    attr(ml$name, "expression") <- 
	expression(make.vslink(family="inverse.gaussian"))
		return(ml)
	} else if (fam=="binomial"|fam=="quasibinomial") { 
			linkfun <- function(mu)  asin(sqrt(mu/NROW(mu)))
# only works if nobs is globally available
# mu/n is always between 0 and 1, otherwise asin would return NaN
			linkinv <- function(eta) NROW(eta)*(sin(eta))^2
			mu.eta 	<- function(eta) NROW(eta)*sin(2*eta)
			valideta <- function(eta) TRUE
		  name     <- "varstab"
		  if (fam=="binomial") {
		    attr(name, "expression") <- 
		expression(make.vslink(family="binomial"))
      } else attr(name, "expression") <- 
		expression(make.vslink(family="quasibinomial"))
			structure(list(linkfun = linkfun, 
			linkinv = linkinv, mu.eta = mu.eta,
			valideta = valideta, name = name), 
			class = "link-glm")
			
	} else if (fam=="quasi") {
		var <- family$varfun	
		if(var=="constant") {
      ml <- make.link("identity")
      attrml <- ml$name
      ml$name <- "varstab"
      attr(ml$name, "name") <- attrml
      attr(ml$name, "expression") <-
        expression(make.vslink(family=quasi(var="constant")))
		  return(ml)
		} else if(var=="mu(1-mu)") {
			linkfun <- function(mu)  asin(sqrt(mu)) 
			# mu must be between -1 and +1
			linkinv <- function(eta) (sin(eta))^2
			mu.eta 	<- function(eta) sin(2*eta)
			valideta <- function(eta) all(abs(eta) <= pi/2)
			name     <- "varstab"
			attr(name, "expression") <-
			 expression(make.vslink(
				family=quasi(var="mu(1-mu)") ))
			structure(list(linkfun = linkfun, 
			linkinv = linkinv, mu.eta = mu.eta, 
			valideta = valideta, name = name),
			class = "link-glm")
		} else if(var=="mu") {
      ml <- make.link("sqrt")
      attrml <- ml$name
      ml$name <- "varstab"
      attr(ml$name, "name") <- attrml
      attr(ml$name, "expression") <- 
		expression(make.vslink(family=quasi(var="mu")))
		  return(ml)
		} else if(var=="mu^2") {
      ml <- make.link("log")
      attrml <- ml$name
      ml$name <- "varstab"
      attr(ml$name, "name") <- attrml
      attr(ml$name, "expression") <- 
		expression(make.vslink(family=quasi(var="mu^2")))
		  return(ml)
		} else if(var=="mu^3") {
      ml <- make.link("1/mu^2")
      attrml <- ml$name
      ml$name <- "varstab"
      attr(ml$name, "name") <- attrml
      attr(ml$name, "expression") <- 
		expression(make.vslink(family=quasi(var="mu^3")))
		  return(ml)
		} 
	} else if(substring(fam,1,17)=="Negative Binomial") {
		env <- new.env(parent = .GlobalEnv)
		# assign(".Theta", fam$theta, envir = env)

    if(!missing(theta))
    .Theta <<- theta      # and in the neg.bin function use the environment from there
    else .Theta <<- 1
    
		linkfun <- function(mu)  asinh(sqrt(mu/.Theta))
		linkinv <- function(eta) (sinh(eta))^(2)*.Theta
		mu.eta 	<- function(eta) (.Theta)*sinh(2*eta)
		valideta <- function(eta) TRUE
		
		environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <-
		 	environment(valideta) <- environment(.Theta) <- env
		name     <- "varstab"
		attr(name, "expression") <- 
		  expression( make.vslink(family="negative.binomial", theta=.Theta) )
		  
		structure( list(linkfun = linkfun,
			linkinv = linkinv, mu.eta = mu.eta, 
		  valideta = valideta, name = name),
			class = "link-glm" )
			
	} else  stop("'family' is invalid")

}

