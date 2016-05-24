`mafc.logit` <- 
function( .m = 2 )
{
	.m <- as.integer(.m)
	if (.m < 2) 
		stop(".m must be an integer > 1")
	linkfun <- function(mu) {
		mu <- pmax(mu, 1/.m +  .Machine$double.eps)
		qlogis((.m * mu - 1)/(.m - 1) ) }  
	linkinv <- function(eta) {
		1/.m  + (.m - 1)/.m * binomial()$linkinv(eta) 
#		.Call("logit_linkinv", eta, PACKAGE = "stats")
		}
	mu.eta <- function(eta) ((.m -1) / .m) *  binomial()$mu.eta(eta)
#		.Call("logit_mu_eta", eta, PACKAGE = "stats")
	valideta <- function(eta) TRUE
	link <- paste("mafc.logit(", .m, ")", sep = "")
	structure(list(linkfun = linkfun,
				  linkinv = linkinv,
				  mu.eta = mu.eta,
				  valideta = valideta, name = link),
				  class = "link-glm")
}

`mafc.probit` <- 
function( .m = 2 )
{
	.m <- as.integer(.m)
	if (.m < 2)
		stop("m must be an integer > 1")
	linkfun <- function(mu) {
		mu <- pmax(mu, 1/.m +  .Machine$double.eps)
		qnorm((.m * mu - 1)/(.m - 1) ) } 
	linkinv <- function(eta) {
		1/.m  + (.m - 1)/.m * pnorm(eta)
		}
	mu.eta <- function(eta) ((.m -1) / .m) * dnorm(eta)
	valideta <- function(eta) TRUE
	link <- paste("mafc.probit(", .m, ")", sep = "")
	structure(list(linkfun = linkfun,
				  linkinv = linkinv,
				  mu.eta = mu.eta,
				  valideta = valideta, name = link),
				  class = "link-glm")
}

`mafc.cloglog` <- 
function( .m = 2 ) 
{
	.m <- as.integer(.m)
    if (.m < 2) 
        stop(".m must be an integer > 1")
    linkfun <- function(mu) {
    	mu <-  pmax(pmin(mu, 1 - .Machine$double.eps), 
    		1/.m + .Machine$double.eps)
        	log(-log((.m - mu * .m)/(.m - 1))) 
     }
     linkinv <- function(eta) {
        	tmp <- 1/.m + ((.m - 1)/.m) * (-expm1(-exp(eta)))
  	     	pmax(pmin(tmp, 1 - .Machine$double.eps), 
        		1/.m + .Machine$double.eps) 
     }
     mu.eta <- function(eta) {
           eta <- pmin(eta, 700)
           pmax(((.m - 1)/.m) * exp(eta) * exp(-exp(eta)),
            .Machine$double.eps)
     }
     valideta <- function(eta) TRUE
     link <- paste("mafc.cloglog(", .m, ")", sep = "")
     structure(list(linkfun = linkfun, linkinv = linkinv, 
         	mu.eta = mu.eta, 
        valideta = valideta, name = link), class = "link-glm")
        }

`mafc.weib` <- function(...)
	mafc.cloglog(...)
	
`mafc.cauchit` <- function( .m = 2) {
	 .m <- as.integer(.m)
    if (.m < 2) 
        stop(".m must be an integer > 1")
	linkfun <- function(mu) {
		mu <- pmax(mu, 1/.m + .Machine$double.eps)
		qcauchy((.m * mu - 1)/(.m - 1))
	}
	linkinv <- function(eta) {
           thresh <- -qcauchy(.Machine$double.eps)
           eta <- pmin(pmax(eta, -thresh), thresh)
           1/.m + (.m - 1)/.m * pcauchy(eta)
    }
     mu.eta <- function(eta) {
     	pmax(((.m - 1)/.m)*  dcauchy(eta), .Machine$double.eps)	}
   	valideta <- function(eta) TRUE
   	link <- paste("mafc.cauchy(", .m, ")", sep = "")
    structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
        valideta = valideta, name = link), class = "link-glm")

}
