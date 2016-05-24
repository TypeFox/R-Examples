exponential <- function( sigma, theta, N ) {

	x <- 0:round( N / 2 )
	return( sigma^2 * exp( -3 * x / theta ) )

} # end of 'exponential' function.

WRSS.exp <- function( params, N, dcov, dat ) {

	sigma <- params[1]
	theta <- params[2]
	acf.fit <- exponential( sigma, theta, N = N)
	num.pairs <- N - dcov$lag
	sum( num.pairs * ( acf.fit - dat$y )^2 )

} # end of 'WRSS.exp' function.

ORSS.exp <- function( params, N, dcov, dat ) {

	sigma <- params[1]
	theta <- params[2]
	acf.fit <- exponential( sigma, theta, N = N) 
	sum(( acf.fit - dat$y )^2 )

} # end of 'ORSS.exp' function.

hg.test <- function( loss1, loss2, plot = FALSE , type) {

    # Arguments (input):
    #
    # 'loss1', 'loss2', numeric vectors of equal length giving the two loss series
    #			(e.g., loss1 = abs( M1 - O ) and loss2 = abs( M2 - O ) ).
    #
    # type says whether the optimization uses WLS or OLS

    # Value (output):
    #
    # numeric vector of length four giving the statistics and p-values calculated by:
    #    (1) fitting model using up to half the maximum lag and (2) setting the 
    #    autocovariances to zero that correspond to small empirical ones.

    # The loss differential series.
    d <- loss1 - loss2

    # length of the series.
    N <- length( d )

    # Calculate the autocovariances up to half the maximum lag.
    dcov <- acf( d, type = "covariance", lag.max = round( N / 2 ), na.action = na.omit, plot = FALSE )

    # Put the lags and estimated autocovariances into a data frame object.
    dat <- data.frame( x = dcov$lag, y = dcov$acf[,, 1] )

    # Get some kind of reasonable starting values for the exponential covariance parameters.
    #theta.start <- ifelse( any( dat$y < 0.5 ), min( dat$x[ dat$y < 0.5 ] ), length( dat$x ) )
	theta.start <- ifelse(dat$y[2]>0,-(dat$y[2] - dat$y[1]), -(0-dat$y[1]))
#ASH:  Changed b/c the theta parameter should correspond to the speed at which the covariance decreases.
# If the second lag is negative, then just compute the slope from the first lag to 0.
    sigma.start <-sqrt( dat$y[ 1 ] )

    # Use 'nls' to try to fit the model to the data.  Use try because often this doesn't work.
	
	if(type=="WLS"){

	    f <- try(nlminb(c( sigma.start, theta.start ),
		    WRSS.exp, N = N, dcov = dcov, dat = dat,lower = c(0,0), upper = c(Inf,Inf)) )

	} else if(type=="OLS"){

	    f <- try(nlminb(c( sigma.start, theta.start ),
		    ORSS.exp, N = N, dcov = dcov, dat = dat,lower = c(0,0), upper = c(Inf,Inf)) )

	}

    # If it worked, find the statistics and p-values.  Otherwise, return NA's.
    if( class( f ) != "try-error" ) {

	xseq <- seq(0, 100, len = 1000)
        co <- f$par

        if( plot ) {

            par( mfrow = c(3, 3) )

            acf( loss1, main = "loss1", xlab = "" )
            acf( loss2, main = "loss2", xlab = "", ylab = "" )
            acf( d, main = "loss1 - loss2", xlab = "", ylab = "" )

            pacf( loss1, main = "" )
            pacf( loss2, main = "", ylab = "" )
            pacf( d, main = "", ylab = "" )

            plot(dat$x, dat$y, xlim = c(0, 100))
            lines(xseq, (co[1]^2) * exp(-3 * xseq / co[2]))
            abline(h = 0, col = 2)

        } # end of if 'plot' stmt.


	# Estimate the mean loss differential.
	m <- mean( d, na.rm = TRUE )

	# Use all lags to estimate the variance from the fitted model.
	d1.var <- (co[1]^2) * exp( -3 * (0:N) / co[2] )
	nd1var <- length( d1.var )

	# Find lags that are small and set to zero (for d2.var).
	id <- c(FALSE,abs(dat$y[2:length(dat$y)]) < 2 / sqrt(N))
#ASH: the first lag, which is the variance should always be included
	nid <- length( id )
	if( nid < nd1var ) id <- c( id, rep( TRUE, nd1var - nid ) )

	d2.var <- d1.var
	d2.var[ id ] <- 0
	d2.var[ length(id):length(d1.var) ] <- 0
	
	# Summing var/covs over the lags
	var1 <- ( d1.var[1] + 2 * sum( d1.var[ 2:length(d1.var) ], na.rm = TRUE) ) / N
	var2 <- ( d2.var[1] + 2 * sum( d2.var[ 2:length(d2.var) ], na.rm = TRUE) ) / N

	# Estimate the two statistics.
	S1 <- m / sqrt( var1 )

	# if all autocovariances are zero, return NA for S2 and p-value 2.
	# Actually, the way I'm doing it now, they will never all be zero.
	if( !all( id ) ) S2 <- m / sqrt( var2 )
	else S2 <- NA

	# Find the associated p-values.  Use a t distribution if N < 30.
	if( N < 30 ) pval1 <- 2 * pt( -abs( S1 ), df = N - 1 )
	else pval1 <- 2 * pnorm( -abs( S1 ) )

	if( !all( id ) ) {

	    if( N < 30 ) pval2 <- 2 * pt( -abs( S2 ), df = N - 1 )
	    else pval2 <- 2 * pnorm( -abs( S2 ) )

	} else S2 <- pval2 <- NA

	out <- c( S1, pval1, S2, pval2 )

    } else out <- rep( NA, 4 )

    names( out ) <- c( "S 1", "pval 1", "S 2", "pval 2" )

    return( out )

} # end of 'hg.test' function.
