locglmfit<-function( xfit, r, m, x, h, returnH = FALSE, link = c( "logit" ),
                     guessing = 0, lapsing = 0, K = 2, p = 1,
                     ker = c( "dnorm" ), maxiter = 50, tol = 1e-6 ) {
#
# Local polynomial estimator for the psychometric function and eta function (psychometric function 
# transformed by link) for binomial data; also returns the hat matrix H. Actual calculations are 
# done in LOCGLMFIT_PRIVATE or LOCGLMFIT_SPARSE_PRIVATE depending on the size of the data set. 
# Here, the data are split into several parts to speed up the calculations.
#
#
#INPUT
#
# xfit - points at which to calculate the estimate pfit
# r    - number of successes at points x
# m    - number of trials at points x 
# x    - stimulus levels 
# h    - bandwidth(s)
#
# OPTIONAL INPUT
#
# returnH  - logical, if TRUE then hat matrix is calculated; default is FALSE
# link     - name of the link function; default is 'logit'
# guessing - guessing rate; default is 0
# lapsing  - lapsing rate; default is 0
# K    - power parameter for Weibull and reverse Weibull link; default is 2
# p        - degree of the polynomial; default is 1
# ker      - kernel function for weights; default is 'dnorm'
# maxiter  - maximum number of iterations in Fisher scoring; default is 50
# tol      - tolerance level at which to stop Fisher scoring; default is 1e-6
#
# OUTPUT
#
# pfit    - value of the local polynomial estimate at points xfit
# etafit  - estimate of eta (link of pfit)
# H       - hat matrix (OPTIONAL)

#### 
# KZ 19-Mar-12
# changed so that in every call to locglmfit the warnings about zero determinant and exceeded number of 
# iterations are displayed only once; that is:
# added a variable Warn which is [0 0] if there are no warnings from private functions, 
# if the first entry is positive, then a warning about too small a bandwidth is displayed,
# if the second entry is positive, then a warning about exceeded number of iterations is displayed; 
# changed the private function acordingly, so that they return the required
# information
####

#### 
# KZ 24-Mar-12
# the threshold value for which sparse matrices are used in calculations 
# was chnaged from 15 to 20 as this improved speed
####

#### 
# KZ 28-Mar-12
# included on.exit function which restores warning settings to their
# original state
####

# First 5 arguments are mandatory
    if( missing("xfit") || missing("r") || missing("m") || missing("x") || missing("h") ) {
        stop("Check input. First 5 arguments are mandatory");
    }

# CHECK ROBUSTNESS OF INPUT PARAMETERS
    if ( !is.vector( xfit ) ) {
        stop("Vector xfit should be a column vector");
    }
    checkdata<-list();
    checkdata[[1]] <- x;
    checkdata[[2]] <- r;
    checkdata[[3]] <- m;
    checkinput( "psychometricdata", checkdata );
    rm( checkdata )
    if ( !is.vector( h ) ) {
        stop("Bandwidths h should be a vector");
    }
    if( !( length( h ) == 1 || length( h ) == length( xfit ) ) ) {
        stop( "Bandwidth h must be either a scalar or a vector with the same number of elements as xfit" );
    }
    checkinput( "linkfunction", link );
    if( length( guessing ) > 1 ) {
        stop( "Guessing rate must be a scalar" );
    }
    if( length( lapsing ) > 1 ) {
        stop( "Lapsing rate must be a scalar" );
    }

    checkinput( 'guessingandlapsing', c(guessing, lapsing) );
   if (link == "weibull" || link == "revweibull"){
	    checkinput( "exponentk", K );
	    }
    pn <- list()
	pn[[1]] <- p
	pn[[2]] <- x
    checkinput( "degreepolynomial", pn );
    checkinput( 'kernel', ker );
    checkinput( 'maxiter', maxiter );
    checkinput( 'tolerance', tol );

# INITIAL VALUES
    split <- 20;

    Lxfit <- length(xfit);
    Lx <- length(x);

	value <- NULL
    pfit <- NULL;
    etafit    <- NULL;
	Warn <- c(0,0)
    if( returnH  ) H <- NULL;

	if( link == "logit"      ||
        link == "probit"     ||
        link == "loglog"     ||
        link == "comploglog" ||
        link == "weibull"    ||
        link == "revweibull" ) {
            	
    			link <- paste( link, "_link_private", sep = "" );
	}

# KZ 28-03-2012 included on.exit routine so that the warning settings are
# restored when the function terminates even if interrupted by the user

warn.current <- getOption("warn")
on.exit(options(warn = warn.current));

options(warn=-1)

# KZ 24-03-12
# changed the threshold value for which sparse matrices 
# are used to 20 (previously 15) to speed up calculations

    if( Lx > 20 ) {
# big data
# First try to load package SparseM
	        options(warn = -1);
        existSparseM <- library( SparseM, logical.return = TRUE );
         if( existSparseM ) {
            fun_estim <- locglmfit_sparse_private;
        }
        else {
            fun_estim <- locglmfit_private;
            message("Package SparseM not installed. No sparse matrices used");
            message("SparseM can be found at CRAN web site http://cran.r-project.org/");
        }
    }
    else {
# small data
        fun_estim <- locglmfit_private;
    }

# SPLIT AND EVALUATION
################################################################ SCALAR h
    if( length( h ) == 1 ) {
# with Hat matrix
        if( returnH ) {
            if( Lxfit <= split ) {
# small x

                value <- fun_estim( xfit, r, m, x, h, returnH, link, guessing,
                       lapsing, K, p, ker, maxiter, tol);
			Warn <- Warn + value$warncount
            }
            else {
# large x
# number of parts into which the fitting is divided
                fLx = floor( Lxfit / split );
# initialise output
                    for( i in 0:(fLx-1) ) {
# part of the fit
                        value1 <- fun_estim( xfit[i*split+c(1:split)], r, m, x, h,
                                returnH, link, guessing, lapsing, K, p, ker,
                                maxiter, tol );
# put the fits together
                        pfit <- c( pfit, value1$pfit );
                        etafit    <- c( etafit,    value1$etafit );
                        H      <- rbind(H,   value1$H );
				Warn <- Warn + value1$warncount
                    }
# final part of the fit
                    if( ( split * fLx ) < Lxfit ) {
                        value1 <- fun_estim( xfit[c(1+split*fLx):Lxfit], r, m, x, h,
                                returnH, link, guessing, lapsing, K, p, ker,
                                maxiter, tol );

# put the fits together
                        pfit <- c( pfit, value1$pfit );
                        etafit    <- c( etafit,    value1$etafit );
                        H      <- rbind(H,   value1$H );
				Warn <- Warn + value1$warncount
                    }
#   values to return

value$pfit <- pfit
value$etafit <- etafit
value$H <- H

				}
        }
        else { # no Hat matrix

            if( Lxfit <= split ) {
# small x

                value <- fun_estim( xfit, r, m, x, h,returnH, link, guessing,
                       lapsing, K, p, ker, maxiter, tol );
			Warn <- Warn + value$warncount

            }
            else {
# large x
# number of parts into which the fitting is divided
                fLx = floor( Lxfit / split );

# initialise output
                for( i in 0:(fLx-1) ) {
# part of the fit

                    value1 <- fun_estim(  xfit[i*split + c(1:split)], r, m, x, h,
                                returnH, link, guessing, lapsing, K, p, ker,
                                maxiter, tol );
# put the fits together
                    pfit <- c( pfit, value1$pfit );
                    etafit    <- c( etafit,    value1$etafit );
			Warn <- Warn + value1$warncount
                }
# final part of the fit
                if( ( split * fLx ) < Lxfit ) {
                    value1 <- fun_estim( xfit[c(1+split*fLx):Lxfit], r, m, x, h,
                            returnH, link, guessing, lapsing, K, p, ker,
                            maxiter, tol );
# put the fits together
                    pfit <- c( pfit, value1$pfit );
                    etafit    <- c( etafit,    value1$etafit );
				Warn <- Warn + value1$warncount
                }
           
#   values to return

value$pfit <- pfit
value$etafit <- etafit

			}
        } # if( returnH )
    }
    
 ################################################################ VECTOR h   
    else { # if( length( h ) == 1 )
    	
# with Hat matrix
        if( returnH ) {
            if ( Lxfit <= split ) {
# small x
                value <- fun_estim( xfit, r, m, x, h, returnH, link, guessing,
                       lapsing, K, p, ker, maxiter, tol );
                  Warn <- Warn + value$warncount
            }
            else {
# large x
# number of parts into which the fitting is divided
                fLx = floor( Lxfit / split );
                    for( i in 0:(fLx-1) ) {
# part of the fit
                        value1 <- fun_estim( xfit[i*split + c(1:split)], r, m, x,
                                h[i*split + c(1:split)], returnH, link,
                                guessing, lapsing, K, p, ker, maxiter,
                                tol );
# put the fits together
                        pfit <- c( pfit, value1$pfit );
                        etafit    <- c( etafit,    value1$etafit );
                        H      <- rbind(H,   value1$H );
				Warn <- Warn + value1$warncount
                    }
# final part of the fit
                    if( ( split * fLx ) < Lxfit ) {
                        value1 <- fun_estim( xfit[c(1+split*fLx):Lxfit], r, m, x,
                                h[c(1+split*fLx):Lxfit], returnH, link,
                                guessing, lapsing, K, p, ker, maxiter,
                                tol );
# put the fits together
                        pfit <- c( pfit, value1$pfit );
                        etafit    <- c( etafit,    value1$etafit );
                        H      <- rbind(H,   value1$H );
				Warn <- Warn + value1$warncount
                    }
#   values to return

value$pfit <- pfit
value$etafit <- etafit
value$H <- H

				}        
        }# end if(retunH)
        else { # no Hat matrix
            if( Lxfit <= split ) {
# small x
                value <- fun_estim( xfit, r, m, x, h, returnH, link, guessing,
                       lapsing, K, p, ker, maxiter, tol );
			Warn <- Warn + value$warncount
            }
            else {
# large x
# number of parts into which the fitting is divided
            fLx = floor( Lxfit / split );

# initialise output
                for( i in 0:(fLx-1) ) {
# part of the fit
                    value1 <- fun_estim( xfit[i*split + c(1:split)], r, m, x,
                            h[i*split + c(1:split)], returnH, link,
                            guessing, lapsing, K, p, ker, maxiter, tol );
# put the fits together
                    pfit <- c( pfit, value1$pfit );
                    etafit    <- c( etafit,    value1$etafit );
				Warn <- Warn + value1$warncount
                }
# final part of the fit
                if( ( split * fLx ) < Lxfit ) {
                    value1 <- fun_estim( xfit[c(1+split*fLx):Lxfit], r, m, x,
                            h[c(1+split*fLx):Lxfit], returnH, link,
                            guessing, lapsing, K, p, ker, maxiter, tol );
# put the fits together
                    pfit <- c( pfit, value1$pfit );
                    etafit    <- c( etafit,    value1$etafit );
				Warn <- Warn + value1$warncount
                }

#   values to return

value$pfit <- pfit
value$etafit <- etafit

			}
        } # if( returnH )
        
    } # if( length( h ) == 1 )

# KZ 19-03-12
# warn user once if there were any warnings generated in the private
# functions

options(warn = warn.current)

if (Warn[1]) warning('Determinant close to 0: bandwidth is too small')
if (Warn[2]) warning("iteration limit reached")

   return( value )
}
