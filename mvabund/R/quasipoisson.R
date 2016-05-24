######################################################################
## quasipoisson a function to create a family object		            ##
######################################################################
## link			= character. possible values:			                      ##
##				  "log", "identity", and "sqrt"		                        ##
## phi			= the nuisance Parameter, var(Y) = phi*E(Y)             ##
##          (if phi = 1 we have a poisson model)	                  ##
######################################################################
## it is allowed not to specify phi                                 ##
## note, that variance, dev.resids, aic are not calculated in 	    ##
## this case.									                                      ##
## passing family=quasiopoisson to manyglm will calculate 	        ##
## the model and the family returned in the list of values		      ##
## will have phi, as well as variance, dev.resids and	aic specified	##
######################################################################

quasipoisson <- function(link = "log", phi=NULL)
{
    linktemp <- substitute(link)
    if (!is.character(linktemp)) {
        linktemp <- deparse(linktemp)
        if (linktemp == "link") {
            warning("use of quasipoisson(link=link) is deprecated\n", 
                domain = NA)
            linktemp <- eval(link)
            if (!is.character(linktemp) || length(linktemp) != 
                1) 
                stop("'link' is invalid", domain = NA)
        }
    }
    okLinks <- c("log", "identity", "sqrt")
    if(linktemp=="varstab"){
        stats     <- make.vslink("quasipoisson")
        linktemp  <- stats$name
        # no attribute, as there is no easy analogue for it.
    } else if (linktemp %in% okLinks) {
        stats <- make.link(linktemp)
    } else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    } else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name)) 
                linktemp <- stats$name
        }
        else {
            stop(gettextf("link \"%s\" not available for quasipoisson family; ",
             "available links are %s",
              linktemp, paste(sQuote(okLinks), collapse = ", ")),
                domain = NA)
        }
    }
	
    validmu <- function(mu) all(mu > 0)
	# check true, in new version with phi?
    env <- new.env(parent = .GlobalEnv)
	  environment(validmu) <- env
	
	if(!is.null(phi)){
			if(!is.numeric(phi)) stop("'phi' must be numeric")
			if(any(phi < 1) ) stop("'phi' must be chosen greater than one")
			
			# assign(".Phi", phi, envir = env)
			.Phi <<- phi
			 environment(.Phi) <- env
			
			variance <- function(mu) {
				mu <- as.matrix(mu)
				length.phi <- length(.Phi)  
				if(length.phi!=1 && length.phi != NCOL(mu)) 
					stop("the number of columns of 'mu' must be the length of '.Phi'")
				ph	<- matrix(rep.int(.Phi, times=nrow(mu)), 
					nrow=nrow(mu),ncol=NCOL(mu),byrow=TRUE)
				mu * ph
			}	
	
			aic <- function(y, n, mu, wt, dev) {
				mu <- as.matrix(mu)
				y <- as.matrix(y)
				nobs <- nrow(y)
				length.phi <- length(.Phi)  
					if(length.phi!=1 && length.phi != NCOL(mu)) 
						stop("the number of columns of 'mu' must be the length of '.Phi'")
				ph	<- matrix(rep.int(.Phi, times=nrow(mu)), 
					nrow=nrow(mu),ncol=NCOL(mu),byrow=TRUE)
        nopoi <- ph[1,] > 1
				qupoi <- which(nopoi) # others are poisson distr.
				aictmp <- numeric(NCOL(mu))
				
			  if(length(qupoi) >0 ){
          ph    <- ph[,qupoi]
				  term  <- lgamma(y[,qupoi] + mu[,qupoi]*(1/(ph-1)))
            - lgamma( mu[,qupoi]*(1/(ph-1)) )
            + ( mu[,qupoi]/(ph-1))* log(1/ph) + y[,qupoi]*log(1-1/ph)
					# term = log likelihood, not yet summed up		
					aictmp[qupoi] <- c( 2 * matrix(1,1,nobs) %*% as.matrix(term * wt) )
				}
				
        if(length(qupoi) < NCOL(mu)) {    # for poisson:
            poi <- which(! nopoi)
            aictmp[poi] <-  -2 * matrix(1,1,nrow(y)) %*% (dpois(y[,poi], mu[,poi], log = TRUE) * wt)
          }
					return(aictmp)
			}
			
			dev.resids <- function(y, mu, wt) { 

				length.phi <- length(.Phi)  
				mu <- as.matrix(mu)
				y <- as.matrix(y)
					if(length.phi!=1 && length.phi != NCOL(mu)) 
						stop("the number of columns of 'mu' must be the length of '.Phi'")
				ph	<- matrix(rep.int(.Phi, times=nrow(mu)), 
					nrow=nrow(mu),ncol=NCOL(mu),byrow=TRUE)
        nopoi <- ph[1,] > 1
				qupoi <- which(nopoi) # others are poisson distr.
        ph <- ph[,qupoi]
        drtmp <- matrix(ncol=NCOL(mu), nrow = nrow(mu))

        if(length(qupoi) >0 ){
      	   # drtmp[,qupoi]  <-  2 * wt * ( (mu-y)[,qupoi]*log(ph-1) + (y-mu)[,qupoi]*log(ph))
      	    isy0 <- y[,qupoi] == 0
      	    yno0 <- y[,qupoi]
      	    yno0[isy0] <- yno0[isy0] + 0.00000000001
      	    muno0      <- mu[,qupoi]

      	    # muno0[isy0] <- muno0[isy0] + 0.00000001
      	    drtmp[,qupoi]  <-  2  * wt *(  ( mu - y )[,qupoi]*( 1 - ph )*log(1 / ph) +
                    lgamma( muno0 / ( ph - 1 ) ) + lgamma( yno0 + yno0 / ( ph - 1 ))
                    -lgamma( yno0 + muno0 / ( ph - 1 ) ) -lgamma( yno0 / ( ph - 1 ) )  )
				}
				
				 if(length(qupoi) < NCOL(mu)) {    # for poisson:
            poi <- which(! nopoi)
            drtmp[,poi]  <-  2 * wt * (y[,poi] * log(ifelse(y[,poi] ==
                      0, 1, y[,poi]/mu[,poi])) - (y[,poi] - mu[,poi]))
        }
				return(drtmp)
			}
		
			 environment(dev.resids) <- environment(aic) <- environment(variance) <- env	
			 
		famname <- paste("quasipoisson(", format(round(phi, 4)), ")", sep = "")	 
	} else { 
		dev.resids <- variance <- NULL
		aic	<- function(y, n, mu, wt, dev) { rep(NA, times=NCOL(y)) }
		famname <-  "quasipoisson" 
	}
	# R's solution:
	# dev.resids <- function(y, mu, wt) 2*wt*(y*log(ifelse(y==0,1,y/mu))-(y-mu))
    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the quasiPoisson family")
        n <- rep.int(1, nobs)
        mustart <- y + 0.1
    })
	
	linkfun 	<- stats$linkfun
	linkinv 	<- stats$linkinv
	mu.eta  	<- stats$mu.eta
	valideta 	<- stats$valideta	
	if(linktemp=="varstab")
		environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- 
			environment(valideta) <- env

    structure(list(family = famname, link = linktemp, 
        linkfun = linkfun, linkinv = linkinv, phi=phi,
		variance = variance, 
        dev.resids = dev.resids, aic = aic, mu.eta = mu.eta, 
        initialize = initialize, validmu = validmu, valideta = valideta), 
        class = c("family.mvabund", "family"))
}

