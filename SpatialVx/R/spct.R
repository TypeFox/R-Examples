spct <- function( d, loc, trend = 0, lon.lat = TRUE,
    dmax = NULL, vgmodel = "expvgram", vgmodel.args = NULL, init,
    alpha = 0.05, alternative = c("two.sided", "less", "greater"), mu = 0, 
    verbose = FALSE, ... ) {

    alternative <- tolower( alternative )
    alternative <- match.arg( alternative )

    if( missing( loc ) ) stop( "spct: must supply loc argument." )

    nloc <- dim( loc )

    if( length( nloc ) != 2 || nloc[ 2 ] != 2 ) stop( "spct: invalid loc argument.  Must have two dimensions." )

    nloc <- nloc[ 1 ]

    out <- list()

    out$data.name <- deparse( substitute( d ) )

    out$loss.differential <- d
    out$nloc <- nloc
    out$trend <- trend

    out$optional.arguments <- list(...)

#     if( trend == "ols" ) {
# 
# 	if( verbose ) cat( "Attempting to find and remove any trend.\n" )
# 
# 	dat <- data.frame( lon = loc[,1], lat = loc[,2], loss.differential = d )
# 	fit <- try( lm( loss.differential ~ lon + lat, data = dat ), silent = !verbose )
# 
# 	if( class( fit ) != "try-error" ) {
# 
# 	    tr <- predict( fit )
# 
# 	    d <- d - tr
# 
# 	    out$loss.differential.detrended <- d
# 
# 	} else warning( "spct: unable to find a trend.  Continuing without de-trending." )
# 
# 	out$trend.fit <- fit
# 
#     } else if( is.numeric( trend ) ) {

    d <- d - trend
    out$loss.differential.detrended <- d


    if( is.null( dmax ) ) {

	if( verbose ) cat( "Attempting to find half of the maximum lags.\n" )

	if( lon.lat ) r <- c( rdist.earth( loc, loc ) )
	else r <- c( rdist( loc, loc ) )

	if( any( is.infinite( r ) ) ) r[ is.infinite( r ) ] <- NA

	if( all( is.na( r ) ) ) stop( "spct: unable to determine half the maximum lags.  Maybe supply a value to dmax?" )

	dmax <- 0.5 * max( r, na.rm = TRUE )

    } # end of if 'dmax' is null stmt.

    if( verbose ) cat( "Finding the empirical variogram.\n" )
    vg <- vgram( loc = loc, y = d, lon.lat = lon.lat, dmax = dmax, ... )

    out$empirical.variogram <- vg

    gamma <- vg$vgram
    h <- vg$d

    # Objective function to optimize.
    # Default uses the exponential variogram, but user can supply their own model
    # so long as it has arguments p (the parameters), h (distances) and '...'

    obfun <- function( p, tau, gamma, model, model.args, ... ) {

	fit <- do.call( model, c( list( p = p, h = tau ), model.args ) )

	out <- sum( (fit - gamma)^2, na.rm = TRUE )

	return( out )

    } # end of 'obfun' function.

    if( missing( init ) ) p <- c( sqrt( gamma[ 1 ] ), ifelse( gamma[ 2 ] > 0, -( gamma[ 2 ] - gamma[ 1 ] ), -( 0 - gamma[ 1 ] ) ) )
    else p <- init

    if( verbose ) cat( "Fitting exponential variogram to empirical variogram.\n" )

    fit <- try( nlminb( p, obfun, tau = vg$d, gamma = vg$vgram, model = vgmodel, model.args = vgmodel.args,
	 lower = c(0, 0), upper = c(Inf, Inf) ), silent = !verbose )

    out$parametric.vgram.fit <- fit

    if( verbose ) cat( "Calculating the mean loss differential value.\n" )

    dbar <- mean( d, na.rm = TRUE )

    out$estimate <- dbar
    out$null.value <- mu

    if( class( fit ) != "try-error" ) {

	if( verbose ) cat( "Performing the test.\n" )

	tmp <- fit$par
	if( vgmodel == "expvgram" ) names( tmp ) <- c( "nugget", "range" )
	out$parameter <- tmp

	y <- do.call( vgmodel, c( list( p = fit$par, h = vg$d ), vgmodel.args ) )
        out$fitted.values <- y

	v <- sqrt( mean( y, na.rm = TRUE ) )
	out$se <- v

	Sv <- ( dbar - mu ) / v

	out$statistic <- Sv

	if( nloc >= 30 ) {

	    dn <- "(Normal Approx.)"

	    if( alternative == "two.sided" ) pval <- 2 * pnorm( -abs( Sv ) )
	    else if( alternative == "less" ) pval <- pnorm( Sv )
	    else if( alternative == "greater" ) pval <- pnorm( Sv, lower.tail = FALSE )
 
	    out$alternative <- alternative
	    out$p.value <- pval

	    z.alpha <- qnorm( 1 - alpha / 2 )
            out$conf.int <- c( dbar - z.alpha * v, dbar + z.alpha * v )

	} else {

	    dn <- paste( "(Student-t with ", nloc - 1, " df)", sep = "" )

	    if( alternative == "two.sided" ) pval <- 2 * pt( -abs( Sv ), df = nloc - 1 )
            else if( alternative == "less" ) pval <- pt( Sv, df = nloc - 1 )
            else if( alternative == "greater" ) pval <- pt( Sv, df = nloc - 1, lower.tail = FALSE )

            out$alternative <- alternative
            out$p.value <- pval

            t.alpha <- qt( 1 - alpha / 2, df = nloc - 1 )
            out$conf.int <- c( dbar - t.alpha * v, dbar + t.alpha * v )

	}

	out$method <- paste( "Hering-Genton Spatial Prediction Comparison Test ", dn, sep = "" )

	class( out ) <- "htest"
	print( out )

    }

    invisible( out )    

} # end of 'spct' function.

expvgram <- function( p, h, ... ) {

    return( p[ 1 ] * ( 1 - exp( - h * p[ 2 ] ) ) )

} # end of 'expvgram' function.
