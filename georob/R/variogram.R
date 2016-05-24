##  ############################################################################

sample.variogram <- function( object, ... ) UseMethod( "sample.variogram" )


## ##############################################################################

sample.variogram.georob <- function(
  object,
  lag.dist.def,
  xy.angle.def = c( 0., 180. ),
  xz.angle.def = c( 0., 180. ),
  max.lag = Inf,
  estimator = c( "qn", "mad", "matheron", "ch" ),
  mean.angle = TRUE, ...
)
{
  
  ## purpose:          georob method for generic sample.variogram
  ##                   
  ## author:           A. Papritz
  ## date:             2015-11-27
  
  ## checking mandatory arguments
  
  if( missing( lag.dist.def ) ) stop(
    "some mandatory arguments are missing"
  )
  
  ## preparing response vector and matrix of coordinates
  
  x <- residuals( object, level = 0 )
  locations <- object[["locations.objects"]][["coordinates"]]
  
  ##computing sample.variogram
  
  sample.variogram.default(
    x, locations, 
    lag.dist.def = lag.dist.def,
    xy.angle.def = xy.angle.def, xz.angle.def = xz.angle.def,
    max.lag = max.lag, estimator = estimator, mean.angle = mean.angle
  )
  
}


## ##############################################################################

sample.variogram.formula <- function(
  object,
  data, subset, na.action, 
  locations,
  lag.dist.def,
  xy.angle.def = c( 0., 180. ),
  xz.angle.def = c( 0., 180. ),
  max.lag = Inf,
  estimator = c( "qn", "mad", "matheron", "ch" ),
  mean.angle = TRUE, ...
)
{
  
  ## purpose:          formula method for generic sample.variogram
  ##                   
  ## author:           A. Papritz
  ## date:             2015-11-27
  
  ## checking mandatory arguments
  
  if( missing( locations ) || missing( lag.dist.def ) ) stop(
    "some mandatory arguments are missing"
  )
  
  # get model frame, response vector, matrix of coordinates
  
  ### build combined formula for response and locations
  
  extended.formula <- update( 
    object,
    paste( as.character( object )[2], as.character( locations )[2], sep = " ~ " )
  )
  
  ## setting-up model frame
  
  cl <- match.call()
  mf <- match.call( expand.dots = FALSE )
  m <- match( 
    c( "data", "subset", "na.action" ),
    names(mf), 0L 
  )
  mf <- mf[c(1L, m)]
  mf[["formula"]] <- extended.formula
  mf[["drop.unused.levels"]] <- TRUE
  mf[[1L]] <- as.name( "model.frame" )
  
  mf <- eval( mf, parent.frame() )
  
  attr( attr( mf, "terms" ), "intercept" ) <- 0
  
  ## preparing response vector and matrix of coordinates
  
  x <- model.response( mf )
  locations <- model.matrix( terms( mf ), mf )
  
  ## computing sample.variogram
  
  sample.variogram.default(
    x, locations, 
    lag.dist.def = lag.dist.def,
    xy.angle.def = xy.angle.def, xz.angle.def = xz.angle.def,
    max.lag = max.lag, estimator = estimator, mean.angle = mean.angle
  )
  
}


## ##############################################################################

sample.variogram.default <- 
  function(
    object,
    locations,
    lag.dist.def,
    xy.angle.def = c( 0., 180. ),
    xz.angle.def = c( 0., 180. ),
    max.lag = Inf,
    estimator = c( "qn", "mad", "matheron", "ch" ),
    mean.angle = TRUE, ...
  )
{
  
  # purpose:          function computes the sample variogram of response
  #                   by various (non-)robust estimators
  #                   
  # author:           A. Papritz
  # date:             2012-04-13
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-04-07 AP correcting error when computing directional variograms
  ## 2015-11-27 default method for generic sample.variogram, checking mandatory arguments
  
  # checking mandatory arguments
  
  if( missing( object ) || missing( locations ) || missing( lag.dist.def ) ) stop(
	"some mandatory arguments are missing"
  )
    
  d2r <- pi / 180.

  estimator <- match.arg( estimator )
  
  # pad missing coordinates
  
  if( ( ndim <- NCOL( locations ) ) < 3 ){
    locations <- cbind( 
      locations, 
      matrix( 0., nrow = NROW( locations ), ncol = 3 - NCOL( locations ) )
    )
  }
  colnames( locations ) <- c( "x", "y", "z" )
  
  # compute lag vectors for all pairs of coordinates
  
  indices.pairs <- combn( NROW( locations ), 2 )
  lag.vectors <- locations[ indices.pairs[2,], ] - locations[ indices.pairs[1,], ]
  
  # reflect lag vectors onto half circle
  
  neg.x <- lag.vectors[, 1] < 0.
  
  lag.vectors[neg.x, 1] <- -lag.vectors[neg.x, 1]
  lag.vectors[neg.x, 2] <- -lag.vectors[neg.x, 2]
  lag.vectors[neg.x, 3] <- -lag.vectors[neg.x, 3]
  
  # compute pairwise differences of responses
  
  t.diff <- object[indices.pairs[2,]] - object[indices.pairs[1,]]
  
  # compute Euclidean distances
  
  t.dist <- sqrt( rowSums( lag.vectors^2 ) )
  
  # compute angles in xy- and xz-planes
  
  xy.angle <- -( atan2( lag.vectors[,2], lag.vectors[,1] ) - pi/2. ) / d2r
  xz.angle <- -( atan2( lag.vectors[,3], lag.vectors[,1] ) - pi/2. ) / d2r

  # define lag class upper limits
  
  if( length( lag.dist.def ) == 1 ) {
    t.lag.limits <- seq( 0, max( c( t.dist ) ) + lag.dist.def, by = lag.dist.def )
  } else {
    t.lag.limits <- lag.dist.def
  }
  
  # group the lag vectors into classes
  
  t.lag.class <- cut( t.dist, breaks = t.lag.limits, include.lowest = TRUE )
  xy.angle.class <- cut( xy.angle, breaks = xy.angle.def, include.lowest = TRUE )
  xz.angle.class <- cut( xz.angle, breaks = xz.angle.def, include.lowest = TRUE )
  
  # join first and last angle classes if end points match and there is
  # more than one angle
  
  # xy-plane
  
  n <- length( xy.angle.def )
  d <- diff( xy.angle.def )
  xy.angle.mid.class <- 0.5 * ( xy.angle.def[-1] + xy.angle.def[-n] )
  if( 
    n > 2 &&
    identical( xy.angle.def[1], 0. ) && 
    identical( xy.angle.def[n], 180. ) &&
    !all( d[1] == d[-1] )
  ){
    
    nb <- 2
    right <- TRUE
    include.lowest <- TRUE
    dig.lab <- 3
    breaks <- c( -diff( xy.angle.def[(n-1):n] ), xy.angle.def[2] )
    
    for (dig in dig.lab:max(12, dig.lab)) {
      ch.br <- formatC(breaks, digits = dig, width = 1)
      if (ok <- all(ch.br[-1L] != ch.br[-nb])) 
      break
    }
    labels <- if(ok){
      paste( 
        if(right){
          "("  
        } else {
          "["
        }, ch.br[-nb], ",", ch.br[-1L], 
        if (right){ 
          "]"
        } else {
          ")"
        }, sep = ""
      )
    } else {
      paste( "Range", seq_len(nb - 1L), sep = "_")
      if( ok && include.lowest ){
        if (right)  substr(labels[1L], 1L, 1L) <- "["
        else substring( labels[nb - 1L], nchar(labels[nb - 1L], "c")) <- "]"
      }
      
    }
    sel <- as.integer( xy.angle.class ) == nlevels( xy.angle.class )
    xy.angle[sel] <- xy.angle[sel] - 180.
    levels( xy.angle.class )[c( 1, nlevels( xy.angle.class )) ] <- labels
    xy.angle.mid.class[1] <- xy.angle.mid.class[1] - (180. - xy.angle.mid.class[n-1])
    xy.angle.mid.class <- xy.angle.mid.class[-(n-1)]
  }
  
  # xz-plane
  
  n <- length(xz.angle.def)
  d <- diff( xz.angle.def )
  xz.angle.mid.class <- 0.5 * ( xz.angle.def[-1] + xz.angle.def[-n] )
  if( 
    n > 2 &&
    identical( xz.angle.def[1], 0. ) && 
    identical( xz.angle.def[n], 180. ) &&
    !all( d[1] == d[-1] )
  ){
    
    nb <- 2
    right <- TRUE
    include.lowest <- TRUE
    dig.lab <- 3
    breaks <- c( -diff( xz.angle.def[(n-1):n] ), xz.angle.def[2] )
    for (dig in dig.lab:max(12, dig.lab)) {
      ch.br <- formatC(breaks, digits = dig, width = 1)
      if (ok <- all(ch.br[-1L] != ch.br[-nb])) 
      break
    }
    labels <- if(ok){
      paste( 
        if(right){
          "("  
        } else {
          "["
        }, ch.br[-nb], ",", ch.br[-1L], 
        if (right){ 
          "]"
        } else {
          ")"
        }, sep = ""
      )
    } else {
      paste( "Range", seq_len(nb - 1L), sep = "_")
      if( ok && include.lowest ){
        if (right)  substr(labels[1L], 1L, 1L) <- "["
        else substring( labels[nb - 1L], nchar(labels[nb - 1L], "c")) <- "]"
      }
      
    }
    sel <- as.integer( xz.angle.class ) == nlevels( xz.angle.class )
    xz.angle[sel] <- xz.angle[sel] - 180.
    levels( xz.angle.class )[c( 1, nlevels( xz.angle.class )) ] <- labels
    xz.angle.mid.class[1] <- xz.angle.mid.class[1] - (180. - xz.angle.mid.class[n-1])
    xz.angle.mid.class <- xz.angle.mid.class[-(n-1)]
  }
  
  t.classes <- list( t.lag.class, xy.angle.class, xz.angle.class )
  
  # compute mean distance, mean lag vector, angle classes and and number of
  # pairs per class
  
  lag.mean <- as.vector( tapply(
      t.dist,
      t.classes,
      mean
    ))
  xy.angle.mean <- as.vector( tapply(
      xy.angle,
      t.classes,
      mean
    ))
  xz.angle.mean <- as.vector( tapply(
      xz.angle,
      t.classes,
      mean
    ))
  xy.angle.centre <- as.vector( tapply(
      xy.angle.class,
      t.classes,
      function( x ) unique( x )
    ))
  xz.angle.centre <- as.vector( tapply(
      xz.angle.class,
      t.classes,
      function( x ) unique( x )
    ))
  
  xy.angle.centre[!is.na(xy.angle.centre)] <- 
    levels( xy.angle.class )[xy.angle.centre[!is.na(xy.angle.centre)]]
  xz.angle.centre[!is.na(xz.angle.centre)] <- 
    levels( xz.angle.class )[xz.angle.centre[!is.na(xz.angle.centre)]]
  
  t.lag.npairs <- as.vector( table( t.classes ) )
  t.lag.select <- lag.mean <= max.lag
  
  # compute semivariance per class
  
  if( estimator %in% c( "mad", "qn" ) ) {
    
    if( estimator == "mad" ) { 
      
      t.gamma <- as.vector( tapply( 
          t.diff, 
          t.classes, 
          mad, center = 0 
        ))
      
    } else {
      
      t.gamma <- as.vector( tapply( 
          t.diff, 
          t.classes, 
          Qn, finite.corr = TRUE
        ))
    } 
    
    t.gamma <- 0.5 * c( t.gamma )^2
    
  } else if ( estimator == "ch" ) {
    
    t.gamma <- as.vector( tapply( 
        t.diff, 
        t.classes, 
        function( x ) 0.5 * mean( sqrt(abs(x)) )^4 / ( 0.457+0.494/length(x) )
      ))
    
  } else if( estimator == "matheron" ) {
    
    t.gamma <- as.vector( tapply( 
        t.diff, 
        t.classes, 
        function( x ) 0.5 * mean( x^2 )
      ))
    
  }
  
  # collect results
  
  t.sel <- t.lag.npairs > 0 & lag.mean <= max.lag
  r.result <- data.frame(
    lag.dist = lag.mean[ t.sel],
    xy.angle = factor( xy.angle.centre[ t.sel ], levels = levels( xy.angle.class ) ),
    xz.angle = factor( xz.angle.centre[ t.sel ], levels = levels( xz.angle.class ) ),
    gamma = t.gamma[ t.sel ],
    npairs = t.lag.npairs[ t.sel ]
  )
    
  # compute lag vectors
  
  # center of circle sectors
  
  if( mean.angle ){
    
    # mean angle of lag pairs
    
    t.aux <- r.result[["lag.dist"]] * sin( xz.angle.mean[t.sel] * d2r )
    r.result[["lag.x"]] <- t.aux * sin( xy.angle.mean[t.sel] * d2r )
    r.result[["lag.y"]] <- t.aux * cos( xy.angle.mean[t.sel] * d2r )
    r.result[["lag.z"]] <- r.result[["lag.dist"]] * cos( xz.angle.mean[t.sel] * d2r )
    
  } else {
    
    t.aux <- r.result[["lag.dist"]] * sin( xz.angle.mid.class[ as.integer( r.result[["xz.angle"]] ) ] * d2r )
    r.result[["lag.x"]] <- t.aux * sin( xy.angle.mid.class[ as.integer( r.result[["xy.angle"]] ) ] * d2r )
    r.result[["lag.y"]] <- t.aux * cos( xy.angle.mid.class[ as.integer( r.result[["xy.angle"]] ) ] * d2r )
    r.result[["lag.z"]] <- r.result[["lag.dist"]] * 
      cos( xz.angle.mid.class[ as.integer( r.result[["xz.angle"]] ) ] * d2r )
    
  }
  
  class( r.result) <-  c( "sample.variogram", "data.frame" )
  
  attr( r.result, "ndim")                <- ndim
  attr( r.result, "lag.dist.def")       <- lag.dist.def
  attr( r.result, "xy.angle.mid.class")  <- xy.angle.mid.class
  attr( r.result, "xz.angle.mid.class")  <- xz.angle.mid.class
  attr( r.result, "estimator" )          <- estimator
    
  return( r.result )  
  
}

## ##############################################################################

summary.sample.variogram <- 
  function( object, ... )
{
  
  ## Summary method for class summary.sample.variogram
  
  ## 2012-11-12 A. Papritz
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  ans <- object[c( "lag.dist", "npairs", "xy.angle", "xz.angle" ) ]
  attr( ans, "estimator" ) <- attr( object, "estimator" )
  class( ans ) <- "summary.sample.variogram"
  return( ans )
  
}


## ##############################################################################

print.summary.sample.variogram <- 
  function( 
    x, digits = max(3, getOption("digits") - 3), ...
  )
{
  
  ## Print method for class summary.sample.variogram
  
  ## 2012-11-12 A. Papritz
  ## 2012-12-18 AP invisible(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  cat( "\nSample variogram estimator:", attr( x, "estimator" ), "\n" )
  
  cat( "\nSummary of lag distances\n" )
  print( summary( x[["lag.dist"]] ) )
  
  cat( "\nSummary of number of pairs per lag and distance classes\n" )
  print( summary( x[["npairs"]] ) )
  
  cat( "\nAngle classes in xy-plane:", levels( x[["xy.angle"]] ), "\n" )
  cat( "Angle classes in xz-plane:", levels( x[["xz.angle"]] ), "\n" )
  
  invisible( x )
  
}

## ##############################################################################

plot.sample.variogram <- 
  function(
    x,
    type = "p", add = FALSE, 
    xlim = c( 0., max( x[["lag.dist"]] ) ),
    ylim = c( 0, 1.1 * max( x[["gamma"]] ) ),
    col,
    pch,
    cex = 0.8,
    xlab = "lag distance", ylab = "semivariance",
    annotate.npairs = FALSE,
    npairs.pos = 3, npairs.cex = 0.7,
    legend = nlevels( x[["xy.angle"]] ) > 1 || nlevels( x[["xz.angle"]] ) > 1,
    legend.pos = "topleft",
    ...
  )
{
  
  ## Plot method for class sample.variogram
  
  ## 2012-12-12 A. Papritz
  ## 2012-12-21 AP correction for using col and pch
  ## 2013-05-12 AP correction for using ...
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  if( !add ) plot( 
    gamma ~ lag.dist, x, type = "n",
    xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...
  )
  
  if( missing( col ) ){
    col <- 1:nlevels( x[["xy.angle"]] )
  } else if( length( col ) < nlevels( x[["xy.angle"]] ) ) stop(
    "number of colors less than number of directions in x-y-plane for which semivariances were computed"
  )
  if( missing( pch ) ){
    pch <- 1:nlevels( x[["xz.angle"]] )
  } else if( length( pch ) < nlevels( x[["xz.angle"]] ) ) stop(
    "number of colors less than number of directions in x-z-plane for which semivariances were computed"
  )
  
  tapply( 
    1:nrow( x ),
    list( x[["xy.angle"]], x[["xz.angle"]] ),
    function( i, x, type, col, pch, cex ){
      points( 
        gamma ~ lag.dist, x, subset = i, 
        type = type,
        col = col[as.numeric(x[i, "xy.angle"])],
        pch = pch[as.numeric(x[i, "xz.angle"])], 
        cex = cex
      )
      
    },
    x = x, type = type, col = col, pch = pch, cex = cex
  )
  
  if( annotate.npairs ){
    with(
      x, 
      text( lag.dist, gamma, npairs, pos = npairs.pos, cex = npairs.cex, col = col ) 
    )
  }
  
  if( legend ){
    legend(
      x = legend.pos, bty = "n", 
      col = c( 
        1:nlevels( x[["xy.angle"]] ), 
        if( nlevels( x[["xz.angle"]] ) > 1 ) rep( 1, nlevels( x[["xz.angle"]] ) ) 
      ),
      pch = c( 
        rep( 1, nlevels( x[["xy.angle"]] ) ), 
        if( nlevels( x[["xz.angle"]] ) > 1 ) 1:nlevels( x[["xz.angle"]] ) 
      ),
      legend = c( 
        paste( "xy.angle:", levels( x[["xy.angle"]] ) ), 
        if( nlevels( x[["xz.angle"]] ) > 1 ) paste( "xz.angle:", levels( x[["xz.angle"]] ) )
      ),
      pt.cex = cex
    )
    
  }
  
  invisible( x )
  
}

## ##############################################################################

fit.variogram.model <-
  function(
    sv,
    variogram.model = c( "RMexp", "RMaskey", "RMbessel", "RMcauchy", 
      "RMcircular", "RMcubic", "RMdagum", "RMdampedcos", "RMdewijsian", "RMfbm",
      "RMgauss", "RMgencauchy", "RMgenfbm", "RMgengneiting", "RMgneiting", "RMlgd",
      "RMmatern", "RMpenta", "RMqexp", "RMspheric", "RMstable",
      "RMwave", "RMwhittle"
    ), 
    param, fit.param = default.fit.param()[names(param)],
	aniso = default.aniso(), fit.aniso = default.fit.aniso(),
    max.lag = max( sv[["lag.dist"]] ),
    min.npairs = 30,
    weighting.method = c( "cressie", "equal", "npairs" ),
    hessian = TRUE,
    verbose = 0,
    ...
  )
{
  
  ## Function to fit a variogram model to a sample variogram generated by
  ## sample.variogram
  
  ## 2012-12-10 A. Papritz
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  ## 2015-11-27 AP checking mandatory arguments, issuing warnings
  ## 2016-02-08 AP correcting error in setting default values for fit.param
  
  ## auxiliary function called by optim to compute objective function
  
  f.aux <- 
    function( 
      adjustable.param, envir, variogram.model, fixed.param, param.name, aniso.name,
      isotropic, param.tf, bwd.tf, lag.vectors, gamma, npairs, 
      weighting.method, d2r, verbose
    )
  {
    
    ## transform variogram and anisotropy parameters back to original scale
    
    param <- c( adjustable.param, fixed.param )[param.name]
    
    param <- sapply(
      param.name,
      function( x, bwd.tf, param.tf, param ) bwd.tf[[param.tf[x]]]( param[x] ),
      bwd.tf = bwd.tf,
      param.tf = param.tf,
      param = param
    )
    names( param ) <- param.name
    
    aniso <- c( adjustable.param, fixed.param )[aniso.name]
    
    aniso <- sapply(
      aniso.name,
      function( x, bwd.tf, param.tf, param ) bwd.tf[[param.tf[x]]]( param[x] ),
      bwd.tf = bwd.tf,
      param.tf = param.tf,
      param = aniso
    )
    names( aniso ) <- aniso.name
    
    ## check whether extra variogram parameters are within allowed bounds and
    ## return an error otherwise
    
    ep <- param.names( model = variogram.model )
    param.bounds <- param.bounds( variogram.model, NCOL( lag.vectors ) )
    ep.param <- param[ep]
    
    if( !is.null( param.bounds ) ) t.bla <- sapply(
      1:length( ep.param ),
      function( i, param, bounds ){
        if( param[i] < bounds[[i]][1] || param[i] > bounds[[i]][2] ) cat(
          "value of parameter '", names( param[i] ), "' outside of allowed range", sep = "" 
        )
        return( NA_real_ )
      }, 
      param = ep.param,
      bounds = param.bounds
    )
    
    ## update param.aniso.item
    
    param.aniso.item <- get( "param.aniso.item", pos = as.environment( envir ) )
    
    param.aniso.item[["param"]] <- param
    param.aniso.item[["aniso"]][["aniso"]] <- aniso
    param.aniso.item[["aniso"]][["sincos"]] <- list(
      co = unname( cos( aniso["omega"] ) ),
      so = unname( sin( aniso["omega"] ) ),
      cp = unname( cos( aniso["phi"] ) ),
      sp = unname( sin( aniso["phi"] ) ),
      cz = unname( cos( aniso["zeta"] ) ),
      sz = unname( sin( aniso["zeta"] ) )
    )
    param.aniso.item[["aniso"]][["rotmat"]] <- with( 
      param.aniso.item[["aniso"]][["sincos"]],
      rbind(
        c(             sp*so,             sp*co,       cp ),
        c( -cz*co + sz*cp*so,  co*sz*cp + cz*so,   -sp*sz ),
        c( -co*sz - cz*cp*so, -cz*co*cp + sz*so,    cz*sp )
      )
    )
    param.aniso.item[["aniso"]][["sclmat"]] <- with(
      param.aniso.item[["aniso"]],
      1. / c( 1., aniso[ c("f1", "f2") ] )
    )
    
    ## output parameters
    
    t.param <- param
    if( !isotropic ) t.param <- c( 
      t.param, param.aniso.item[["aniso"]][["aniso"]] / c( rep( 1., 2 ), rep( d2r, 3 ) )
    )
    if( verbose > 0 ) {
      cat(
        format(
          signif( t.param[names( adjustable.param)], digits = 7 ),
          scientific = TRUE, width = 14 ),
        sep = ""
      )
    }
    
    ## compute the semivariance
        
    t.model <- f.aux.gamma(
      lag.vectors, 
      variogram.model, param, param.aniso.item[["aniso"]] 
    )
    
    if( identical( class( t.model ), "try-error" ) || any( is.na( t.model ) ) ){
      warning( "there were errors: call function with argument 'verbose' > 1" )
      if( verbose > 1 ) cat( "\nan error occurred when computing semivariances\n" )
      return( NA )
    }
    
    ## compute weights
    
    t.weights <- rep( 1., nrow( lag.vectors ) )
    if( weighting.method != "equal" ) t.weights <- npairs
    if( weighting.method == "cressie" ) t.weights <- t.weights / t.model^2
    
    ##  ... and compute the objective function that is minimized
    
    t.res <- gamma - t.model
    sse <- sum( t.weights * t.res^2 )
    
    if( verbose > 0 ) cat(
      format(
        signif( sse, digits = 7 ),
        scientific = TRUE, width = 14
      ), "\n", sep = ""
    )
    
    ## store copies of model semivariance, residuals and weights
    
    param.aniso.item[["fitted"]] <- t.model
    param.aniso.item[["residuals"]] <- t.res
    param.aniso.item[["weights"]] <- t.weights
    
    ## store param.aniso.item
    
    assign( "param.aniso.item", param.aniso.item, pos = as.environment( envir ) )
    
    return( sse )
    
  }
  
  ## begin of main body of function
  
  ## check whether all mandatory arguments have been provided
  
  if( missing( sv ) || missing( param ) ) stop( 
    "some mandatory arguments are missing" 
  )
  
  d2r <- pi / 180.
      
  ## match arguments
  
  cl <- match.call()

  variogram.model <- match.arg( variogram.model )  
  weighting.method = match.arg( weighting.method )
  
  ## match names of param, aniso, fit.param, fit.aniso
  
  tmp <- names( param )
  tmp <- sapply(tmp, function(x, choices){
      match.arg(x, choices)
    },
    choices = names( default.fit.param() )
  )
  names( param ) <- tmp
  
  if( !missing( fit.param ) ){
    tmp <- names( fit.param )
    tmp <- sapply(tmp, function(x, choices){
        match.arg(x, choices)
      },
      choices = names( default.fit.param() )
    )
    names( fit.param ) <- tmp
    fit.param <- fit.param[names( fit.param ) %in% names( param )]
  }
  
  if( !missing( aniso ) ){
    tmp <- names( aniso )
    tmp <- sapply(tmp, function(x, choices){
        match.arg(x, choices)
      },
      choices = names( default.aniso() )
    )
    names( aniso ) <- tmp
  }
  
  if( !missing( fit.aniso ) ){
    tmp <- names( fit.aniso )
    tmp <- sapply(tmp, function(x, choices){
        match.arg(x, choices)
      },
      choices = names( default.aniso() )
    )
    names( fit.aniso ) <- tmp
  }

  ## set snugget to zero if snugget has not been specified or if there are
  ## no replicated observations
  
  if( !"snugget" %in% names( param ) ){
    param["snugget"] <- 0.
    fit.param["snugget"] <- FALSE
  }
  
  ##  check whether fitting of chosen variogram model is implemented and
  ##  return names of extra parameters (if any)
  
  ep <- param.names( model = variogram.model )
  
  ## check names of initial variogram parameters and flags for fitting
  
  param.name <- c( "variance", "snugget", "nugget", "scale", ep )
  
  if( !all( names( param ) %in% param.name ) ) stop( 
    "error in names of initial values of variogram parameters" 
  )
  
  if( !all( param.name  %in% names( param ) ) ) stop( 
    "no initial values provided for parameter(s) '", 
    param.name[ !param.name %in% names( param ) ], "'"
  )
  
  if( !all( names( fit.param ) %in% param.name ) ) stop( 
    "error in names of control flags for fitting variogram parameters" 
  )
  
  if( length( param ) != length( fit.param ) || 
    !all( names( fit.param ) %in% names( param ) )
  ) stop( 
    "names of variogram parameters and control flags for fitting do not match" 
  )
  
  if( !all( is.numeric( param ) ) ) stop(
    "initial values of variogram parameters must be of mode 'numeric'"
  )
  if( !all( is.logical( fit.param ) ) ) stop(
    "fitting control flags of variogram parameters must be of mode 'logical'"
  )
  
  ##  rearrange initial variogram parameters
  
  param <- param[param.name]
  
  ## check whether intitial values of variogram parameters are valid
  
  if( param["variance"] < 0. ) stop("initial value of 'variance' must be positive" )
  if( param["snugget"] < 0. )  stop("initial value of 'snugget' must be positive" )
  if( param["nugget"] < 0. ) stop("initial value of 'nugget' must be positive" )
  if( param["scale"] <= 0. ) stop("initial value of 'scale' must be positive" )
  
  param.bounds <- param.bounds( variogram.model, attr( sv, "ndim" ) )
  ep.param <- param[ep]
  
  if( !is.null( param.bounds ) ) t.bla <- sapply(
    1:length( ep.param ),
    function( i, param, bounds ){
      if( param[i] < bounds[[i]][1] || param[i] > bounds[[i]][2] ) stop(
        "initial value of parameter '", names( param[i] ), "' outside of allowed range" 
      )
    }, 
    param = ep.param,
    bounds = param.bounds
  )
  
  
  ##  rearrange and check flags controlling variogram parameter fitting 
  
  fit.param <- fit.param[param.name]
  
  if( 
    variogram.model %in% (t.models <- c( "RMfbm" ) ) && 
    ( 
      all( fit.param[c( "variance", "snugget", "scale" ) ] ) ||
      all( fit.param[c( "variance", "scale" ) ] ) 
    )
  ) stop( 
    "'variance', 'scale' (and 'snugget') cannot be fitted simultaneously for variograms ",
    paste( t.models, collapse = " or "), "; \n  'scale' parameter must be fixed"
  )
  
  ##  preparation for variogram parameter transformations
  
  all.param.tf <- param.transf()
  fwd.tf       <- fwd.transf()  
  bwd.tf       <- bwd.transf()  
  
  t.sel <- match( param.name, names( all.param.tf ) )
  
  if( any( is.na( t.sel ) ) ){
    stop( "transformation undefined for some variogram parameters" )
  } else {
    param.tf <- all.param.tf[t.sel]
  }
  param.tf <- sapply(
    param.tf,
    function( x ) if( length(x) > 1L ) x[variogram.model] else x
  )
  names( param.tf ) <- param.name 
  
  ##  transform initial variogram parameters
  
  transformed.param <- sapply(
    param.name,
    function( x, fwd.tf, param.tf, param ) fwd.tf[[param.tf[x]]]( param[x] ),
    fwd.tf = fwd.tf,
    param.tf = param.tf,
    param = param
  )
  
  names( transformed.param ) <- param.name 
  
  ## check whether isotropic variogram is fitted
  
  isotropic <- missing( aniso ) && missing( fit.aniso )
  
  ## check names of initial anisotropy parameters and flags for fitting
  
  aniso.name <- c( "f1", "f2", "omega", "phi", "zeta" )
  
  if( !all( names( aniso ) %in% aniso.name ) ) stop( 
    "error in names of initial values of anisotropy parameters" 
  )
  
  if( !all( aniso.name  %in% names( aniso ) ) ) stop( 
    "no initial values provided for parameter(s) '", 
    aniso.name[ !aniso.name %in% names( aniso ) ], "'"
  )
  
  if( !all( names( fit.aniso ) %in% aniso.name ) ) stop( 
    "error in names of control flags for fitting  anisotropy parameters"
  )
  
  if( length( aniso ) != length( fit.aniso ) || 
    !all( names( fit.aniso ) %in% names( aniso ) )
  ) stop( 
    "names of anisotropy parameters and control flags for fitting do not match" 
  )
  
  if( !all( is.numeric( aniso ) ) ) stop(
    "initial values of anisotropy parameters must be of mode 'numeric'"
  )
  if( !all( is.logical( fit.aniso ) ) ) stop(
    "fitting control flags of anisotropy parameters must be of mode 'logical'"
  )
  
  ##  rearrange initial anisotropy parameters
  
  aniso <- aniso[aniso.name]
  
  ## check whether intitial values of anisotropy parameters are valid
  
  if( aniso["f1"] < 0. ||  aniso["f1"] > 1. ) stop(
    "initial value of parameter 'f1' must be in [0, 1]" 
  )
  if( aniso["f2"] < 0. ||  aniso["f1"] > 1. ) stop(
    "initial value of parameter 'f2' must be in [0, 1]" 
  )
  if( aniso["omega"] < 0. ||  aniso["omega"] > 180. ) stop(
    "initial value of parameter 'omega' must be in [0, 180]" 
  )
  if( aniso["phi"] < 0. ||  aniso["phi"] > 180. ) stop(
    "initial value of parameter 'phi' must be in [0, 180]" 
  )
  if( aniso["zeta"] < -90. ||  aniso["zeta"] > 90. ) stop(
    "initial value of parameter 'zeta' must be in [-90, 90]" 
  )

  ## adjust default initial values of anisotropy parameters if these are
  ## fitted
  
  if( fit.aniso["omega"] && identical( aniso["f1"], 1. ) ) aniso["f1"] <- aniso["f1"] - sqrt( .Machine$double.eps )
  if( fit.aniso["phi"] ){
    if( identical( aniso["f1"], 1. ) ) aniso["f1"] <- aniso["f1"] - 0.0001
    if( identical( aniso["f2"], 1. ) ) aniso["f2"] <- aniso["f2"] - 0.0001
  }
  if( fit.aniso["zeta"] && identical( aniso["f2"], 1. ) ) aniso["f2"] <- aniso["f2"] - 0.0001

  ##  rearrange and check flags controlling anisotropy parameter fitting 
  
  fit.aniso <- fit.aniso[aniso.name]
  
  ##  preparation for anisotropy parameter transformations
  
  t.sel <- match( aniso.name, names( all.param.tf ) )
  
  if( any( is.na( t.sel ) ) ){
    stop( "transformation undefined for some anisotropy parameters" )
  } else {
    aniso.tf <- all.param.tf[t.sel]
  }
  aniso.tf <- sapply(
    aniso.tf,
    function( x ) if( length(x) > 1L ) x[variogram.model] else x
  )
  names( aniso.tf ) <- aniso.name 
  
  ##  convert angles from degrees to radian
  
  aniso[c("omega", "phi", "zeta" )] <- aniso[c("omega", "phi", "zeta" )] * d2r
  
  ##  transform initial anisotropy parameters
  
  transformed.aniso <- sapply(
    aniso.name,
    function( x, fwd.tf, param.tf, param ){
      fwd.tf[[param.tf[x]]]( param[x] )
    },
    fwd.tf = fwd.tf,
    param.tf = aniso.tf,
    param = aniso
  )
  names( transformed.aniso ) <- aniso.name 
  
  param.tf <- c( param.tf, aniso.tf )
  
  ## scaling factor for parameters
  
  ## select lag distances that are used for fitting
  
  t.lag.select <- sv[["lag.dist"]] <= max.lag & sv[["npairs"]] >= min.npairs
  
  if( verbose > 0 ) cat( 
    format(
      c( param.name[fit.param], aniso.name[fit.aniso], "sse" ),
      width = 14, justify= "right"
    ), "\n", sep = ""
  )
  
  ##  create environment to store items required to compute likelihood and
  ##  estimating equations that are provided by
  ##  likelihood.calculations
  
  envir <- new.env()
  
  ##  initialize values of variogram parameters stored in the environment
  
  param.aniso.item <- list(
    param = rep( -1., length( param.name ) ),
    aniso = list(
      isotropic = isotropic,
      aniso = rep( -1., length( aniso.name ) ),
      fit.aniso = fit.aniso
    )
  )  
  assign( "param.aniso.item", param.aniso.item, pos = as.environment( envir ) )
  
  ## fit the model
  
  r.fit <- optim(
    par = c( 
      transformed.param[ fit.param ], 
      transformed.aniso[ fit.aniso ] 
    ),
    fn = f.aux,
    hessian = hessian,
    ...,
    envir = envir,
    variogram.model = variogram.model,
    fixed.param = c( 
      transformed.param[ !fit.param ], 
      transformed.aniso[ !fit.aniso ]
    ),
    param.name = param.name, 
    aniso.name = aniso.name,
    isotropic = isotropic,
    param.tf = param.tf,
    bwd.tf = bwd.tf,
    lag.vectors = as.matrix(
      sv[t.lag.select, c( "lag.x", "lag.y", "lag.z" )]
    ),
    gamma = sv[["gamma"]][t.lag.select],
    npairs = sv[["npairs"]][t.lag.select],
    weighting.method = weighting.method,
    d2r = d2r,
    verbose = verbose
  )
  
  ## get param.aniso.item
  
  param.aniso.item <- get( "param.aniso.item", pos = as.environment( envir ) )
  param.aniso.item[["aniso"]][["aniso"]] <- 
    param.aniso.item[["aniso"]][["aniso"]] / c( rep( 1., 2 ), rep( d2r, 3 ) )
  
  ## collect results
  
  r.result <- list(
    sse = r.fit[["value"]],
    variogram.model = variogram.model,
    param = param.aniso.item[["param"]], fit.param = fit.param,
    aniso = param.aniso.item[["aniso"]],
    param.tf = param.tf,
    fwd.tf = fwd.tf,
    bwd.tf = bwd.tf,
    converged = if( sum( c( fit.param, fit.aniso ) ) == 0 ){ 
      NA
    } else {
      r.fit[["convergence"]] == 0
    },
    convergence.code = r.fit[["convergence"]],      
    iter = r.fit[["counts"]],
    call = cl,
    residuals = param.aniso.item[["residuals"]],
    fitted = param.aniso.item[["fitted"]],
    weights = param.aniso.item[["weights"]]
  )
  if( hessian ) r.result[["hessian"]] <- r.fit[["hessian"]]
  
  class( r.result ) <- "fitted.variogram"
  
  attr( r.result, "lag.dist.def" )      <- attr( sv, "lag.dist.def" )
  attr( r.result, "xy.angle.mid.class" ) <- attr( sv, "xy.angle.mid.class" )
  attr( r.result, "xz.angle.mid.class" ) <- attr( sv, "xz.angle.mid.class" )
  
  return( r.result )
}

## ##############################################################################

print.fitted.variogram <- 
  function( 
    x, digits = max(3, getOption("digits") - 3), ...
  )
{
  
  ## print method for fitted.variogram
  
  ## 2012-04-13 A. Papritz
  ## 2012-12-18 AP invisible(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  cat("\n")
  cat( "Variogram: ", x[["variogram.model"]], "\n" )
  param <- x[["param"]]
  names( param ) <- ifelse(
    x[["fit.param"]],
    names( param ),
    paste( names( param ), "(fixed)", sep = "" )
  )
  print( 
    format( param, digits = digits ), print.gap = 2, 
    quote = FALSE
  )
  
  ## print anisotropy parameters
  
  if( !x[["aniso"]][["isotropic"]] ){
    
    cat("\n")
    cat( "Anisotropy parameters: ", "\n" )
#     aniso <- x[["aniso"]][["aniso"]] * c( rep(1, 2), rep( 180./pi, 3 ) )
    aniso <- x[["aniso"]][["aniso"]]
    names( aniso ) <- ifelse(
      x[["aniso"]][["fit.aniso"]],
      names( aniso ),
      paste( names( aniso ), "(fixed)", sep = "" )
    )
    print( 
      format( aniso, digits = digits ), print.gap = 2, 
      quote = FALSE
    )
    
  }
  
  invisible( x )
  
}

## ##############################################################################

summary.fitted.variogram <- 
  function (
    object, correlation = FALSE,
    signif = 0.95,
    ...
  )
{
    
  ## summary method for fitted.variogram
  
  ## 2012-12-10 A. Papritz
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-04-07 AP changes for fitting anisotropic variograms

  ans <- object[c(
    "call", "residuals", "weights", "converged", "convergence.code", 
    "iter", "sse", "variogram.model", "fit.param"
  )]
  
  if( !object[["aniso"]][["isotropic"]] ) ans[["fit.param"]] <- c( 
    ans[["fit.param"]], object[["aniso"]][["fit.aniso"]]
  )
      
  ans[["param"]] <- as.matrix( object[["param"]], ncol = 1 )
  
  if( !object[["aniso"]][["isotropic"]] ) ans[["param"]] <- rbind( 
    ans[["param"]],
    as.matrix( object[["aniso"]][["aniso"]], ncol = 1 )
  )
  
  colnames( ans[["param"]] ) <- "Estimate"
  
  ## compute confidence intervals of variogram parameters from observed
  ## Fisher information matrix (Gaussian REML only)
  
  if( !is.null( object[["hessian"]] ) ){
    
    ## initialization
    
    cor.tf.param <- cov.tf.param <- matrix( 
      NA, nrow = nrow( object[["hessian"]] ), ncol = nrow( object[["hessian"]] ),
      dimnames = dimnames( object[["hessian"]] )
    )
    
    se <- rep( NA, nrow( object[["hessian"]] ) )
    names( se ) <- rownames( object[["hessian"]])
    
    ci <- matrix( NA, nrow = nrow( ans[["param"]] ), ncol = 2 )
    colnames( ci ) <- c( "Lower", "Upper" )
    rownames( ci ) <- rownames( ans[["param"]] )
    
    ## select parameters that are not on boundary of parameter space
    
    sr  <- !apply( object[["hessian"]], 1, function( x ) all( is.na( x ) ) )
    
    if( sum( sr ) > 0 ){
      
      t.chol <- try( chol( object[["hessian"]][sr, sr] ), silent = TRUE )
      
      if( !identical( class( t.chol ), "try-error" ) ){
        
        ## compute covariance matrix of fitted transformed parameters
        
        cov.tf.param[sr, sr] <- chol2inv( t.chol )
        
        ## correlation matrix and standard errors of fitted transformed
        ## parameters
        
        cor.tf.param[sr, sr] <- cov2cor( cov.tf.param[sr, sr] )
        
        se[sr] <- sqrt( diag( cov.tf.param )[sr] )
        
        ## compute confidence interval on original scale of parameters
        
        sel.names <- names( object[["param"]][object[["fit.param"]]] )
        if( !object[["aniso"]][["isotropic"]] ) sel.names <- c( 
          sel.names,
          names( object[["aniso"]][["aniso"]][object[["aniso"]][["fit.aniso"]]] )
        )
        sel.names <- sel.names[sr]
        
        ci[sel.names, ] <- t( 
          sapply(
            sel.names,
            function( x, param, se, param.tf, trafo.fct, inv.trafo.fct ){
              inv.trafo.fct[[param.tf[x]]]( 
                trafo.fct[[param.tf[x]]]( param[x] ) + 
                c(-1, 1) * se[x] * qnorm( (1-signif)/2., lower.tail = FALSE ) 
              )
            },
            param         = c( object[["param"]], object[["aniso"]][["aniso"]] ),
            se            = se,
            param.tf      = object[["param.tf"]],
            trafo.fct     = object[["fwd.tf"]],
            inv.trafo.fct = object[["bwd.tf"]]
          )
        )
        ans[["param"]] <- cbind( ans[["param"]], ci )
        if( correlation ) ans[["cor.tf.param"]] <- cor.tf.param
        
      } else {
        warning(
          "Hessian not positive definite:",
          "\nconfidence intervals of variogram parameters cannot be computed"        
        )
      }
    } 
  }
    
  class( ans ) <- c( "summary.fitted.variogram" )
  
  ans
}

## ##############################################################################

print.summary.fitted.variogram <- 
  function (
    x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"),
    ...
  ) 
{
    
  ## print.summary method for fitted.variogram
  
  ## 2012-04-13 A. Papritz
  ## 2012-12-18 AP invisible(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  
  cat("\nCall:")
  cat( paste( deparse(x[["call"]]), sep = "\n", collapse = "\n"),  "\n", sep = "" )
  
  if( is.na( x[["converged"]] ) ){
    cat( "\nEstimation with fixed variogram parameters\n" )
    
  } else {
    
    if(!(x[["converged"]])) {
      cat( 
        "\nAlgorithm did not converge, diagnostic code: ", 
        x[["convergence.code"]], "\n"
      )
    } else {
      cat(
        "\nConvergence in", x[["iter"]][1], "function and", 
        x[["iter"]][2], "Jacobian/gradient evaluations\n"
      )
    }
    
    cat(
      "\nResidual Sum of Squares:", 
      x[["sse"]], "\n"
    )
    
  }
    
  resid <- x[["residuals"]]
  cat( "\nResiduals (epsilon):\n")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <- structure( quantile(resid), names = nam )
  print( rq, digits = digits, ...)
  
  cat( "\nVariogram: ", x[["variogram.model"]], "\n" )
  rownames( x[["param"]] ) <- ifelse(
    x[["fit.param"]],
    rownames( x[["param"]] ),
    paste( rownames( x[["param"]] ), "(fixed)", sep = "" )
  )
  printCoefmat(
    x[["param"]], digits = digits, signif.stars = FALSE, ...
  )
  
  if( !is.null( x[["cor.tf.param"]] ) ){
    
    correl <- x[["cor.tf.param"]]
    p <- NCOL(correl)
    if( p > 1 ){
      cat("\nCorrelation of (transformed) variogram parameters:\n")
      correl <- format(round(correl, 2), nsmall = 2, 
        digits = digits)
      correl[!lower.tri(correl)] <- ""
      print(correl[-1, -p, drop = FALSE], quote = FALSE)
    }
    
  }
    
  invisible( x )
}

## ##############################################################################
plot.georob <- 
  function(
    x, what = c( "variogram", "covariance", "correlation", 
      "ta", "sl", "qq.res", "qq.ranef" ),
    add = FALSE,
    lag.dist.def, 
    xy.angle.def = c( 0., 180. ),
    xz.angle.def = c( 0., 180. ),
    max.lag = Inf,
    estimator = c( "mad", "qn", "ch", "matheron" ),
    mean.angle = TRUE,
    level = what != "ta", 
    smooth = what == "ta" || what == "sl",
    id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
    label.pos = c(4,2),    
    col, pch, xlab, ylab, main, lty = "solid", 
    ...    
  )
{
  
  ## Function plots the graph of a variogram fitted by f.georob
  
  ## 2012-12-11 A. Papritz
  ## 2012-12-21 AP correction for using col and pch 
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-05-08 AP changes for plotting covariances and correlations
  ## 2015-11-27 AP changes for Tukey-Anscombe, QQnorm plots and for plotting variograms
  
  x[["na.action"]] <- NULL
  
  estimator <- match.arg( estimator )
  what <- match.arg( what )
  
  ## labelling of points in residual diagnostic plots (taken from plot.lmrob())
  
  n <- length(x[["fitted.values"]])
  if (is.null(id.n)){
    id.n <- 0
  } else {
    id.n <- as.integer(id.n)
    if(id.n < 0L || id.n > n) stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
  }
  if(id.n > 0L){ ## label the largest residuals
    if(is.null(labels.id))
    labels.id <- paste(1L:n)
    iid <- 1L:id.n
    ##    show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
    ## if(any(show[2L:3L]))
    ##     show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
    text.id <- function(x, y, ind, adj.x = TRUE) {
      labpos <-
      if(adj.x) label.pos[1+as.numeric(x > mean(range(x)))] else 3
      text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE,
        pos = labpos, offset = 0.25)
    }
  }
  yh <- x[["fitted.values"]]
  
  switch(
    what,
    ta = {
      
      ## Tukey-Anscombe plot
      
      if( missing( col ) ) col <- 1; if( missing( pch ) ) pch <- 1
      if( missing( xlab ) ) xlab <- "Fitted values"
      if( missing( ylab ) ) ylab <- "Residuals"
      if( missing( main ) ) main <- "Residuals vs. Fitted"
      r <- residuals( x, level = level )
      if( !add ){
        plot( yh, r, col = col, pch = pch,
          xlab = xlab, ylab = ylab, main = main, ...
        )
      } else {
        points( yh, r, col = col, pch = pch, ... )
      }
      if( smooth ){
        lines( loess.smooth( yh, r, ... ), col = "red", lty = 1, ... )
      }
      if(id.n > 0){   ## adapted from plot.lmrob()
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        y.id <- r[show.r]
        y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
        text.id(yh[show.r], y.id, show.r)
      }
      
    },
    sl = {
      
      ## scale.location plot
      
      if( missing( col ) ) col <- 1; if( missing( pch ) ) pch <- 1
      if( missing( xlab ) ) xlab <- "Fitted values"
      if( missing( ylab ) ) ylab <- "Sqrt of abs(Residuals)"
      if( missing( main ) ) main <- "Sqrt of abs(Residuals) vs. Fitted Values"
      r <- sqrt( abs( residuals( x, level = level ) ) )
      if( !add ){
        plot( yh, r, col = col, pch = pch,
          xlab = xlab, ylab = ylab, main = main, ...
        )
      } else {
        points( yh, r, col = col, pch = pch, ... )
      }
      if( smooth ){
        lines( loess.smooth( yh, r, ... ), col = "red", lty = 1 )
      }
      if(id.n > 0) {   ## adapted from plot.lmrob()
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        y.id <- r[show.r]
        y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
        text.id(yh[show.r], y.id, show.r)
      }
      
    },
    qq.res = {
      
      ## qqnorm of standardized residuals (eps or eps + b, depending of
      ## value of level passed in ...)
      
      if( missing( col ) ) col <- 1; if( missing( pch ) ) pch <- 1
      if( missing( xlab ) ) xlab <- "Theoretial quantiles"
      if( missing( ylab ) ) ylab <- "Standardized residuals"
      if( missing( main ) ) main <- "Normal Q-Q standardized residuals"
      r <- rstandard( x, level = level )
      tmp <- qqnorm( r, col = col, pch = pch,
        xlab = xlab, ylab = ylab, main = main, 
        plot.it = !add, ...
      )
      if( add ) points( y~x, tmp, col = col, pch = pch, ... )
      if(id.n > 0){   ## adapted from plot.lmrob()
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        text.id(tmp$x[show.r], tmp$y[show.r], show.r)
      }
    },
    qq.ranef = {
      
      ## qqnorm of standardized random effects
      
      if( missing( col ) ) col <- 1; if( missing( pch ) ) pch <- 1
      if( missing( xlab ) ) xlab <- "Theoretial quantiles"
      if( missing( ylab ) ) ylab <- "Standardized random effects"
      if( missing( main ) ) main <- "Normal Q-Q standardized random effects"
      tmp <- qqnorm( ranef( x, standard = TRUE ), 
        col = col, pch = pch,
        xlab = xlab, ylab = ylab, main = main, 
        plot.it = !add, ...
      )
      if( add ) points( y~x, tmp, col = col, pch = pch, ... )
      
    },
    {
      
      ## plotting variogram, covariance or correlation
      
      ## compute and plot sample variogram
      
      if( !missing( lag.dist.def ) ){
        
        ## compute and plot sample variogram of regression residuals
        
        r.sv <- sample.variogram(
          x, 
          lag.dist.def = lag.dist.def,
          xy.angle.def = xy.angle.def, 
          xz.angle.def = xz.angle.def,
          max.lag = max.lag,
          estimator = estimator,
          mean.angle = mean.angle
        )
        
        if( missing( col ) ){
          col <- 1:nlevels( r.sv[["xy.angle"]] )
        } else if( length( col ) < nlevels( r.sv[["xy.angle"]] ) ) stop(
          "number of colors less than number of directions in x-y-plane for which semivariances are computed"
        )
        if( missing( pch ) ){
          pch <- 1:nlevels( r.sv[["xz.angle"]] )
        } else if( length( pch ) < nlevels( r.sv[["xz.angle"]] ) ) stop(
          "number of colors less than number of directions in x-z-plane for which semivariances are computed"
        )
        if( missing( xlab ) ) xlab <- "lag distance"
        if( missing( ylab ) ) ylab <- "semivariance"
        
        plot( 
          r.sv, add = add, col = col, pch = pch, xlab = xlab, ylab = ylab, ... 
        )
        
        xmax <- max( r.sv[["lag.dist"]] )
        xy.angle.mid.class <- attr( r.sv, "xy.angle.mid.class" )
        xz.angle.mid.class <- attr( r.sv, "xz.angle.mid.class" )
        
      } else {
        
        ## setup window for plotting variogram, covariance, correlation model
        
        if( is.finite( max.lag ) ){
          xmax <- max.lag
        } else {
          xmax <- sqrt( max( rowSums( x[["locations.objects"]][["lag.vectors"]]^2 ) ) )
        }
        
        ymax <- 1.1 * switch(
          what,
          correlation = 1.,
          sum( x[["param"]][c("variance", "snugget", "nugget")] )
        )
        
        ## see sample.variogram.default
        
        # join first and last angle classes if end points match and there is
        # more than one angle
        
        # xy-plane
        
        n <- length( xy.angle.def )
        d <- diff( xy.angle.def )
        xy.angle.mid.class <- 0.5 * ( xy.angle.def[-1] + xy.angle.def[-n] )
        if( 
          n > 2 &&
          identical( xy.angle.def[1], 0. ) && 
          identical( xy.angle.def[n], 180. ) &&
          !all( d[1] == d[-1] )
        ){
          
          xy.angle.mid.class[1] <- xy.angle.mid.class[1] - (180. - xy.angle.mid.class[n-1])
          xy.angle.mid.class <- xy.angle.mid.class[-(n-1)]
        }
        
        # xz-plane
        
        n <- length(xz.angle.def)
        d <- diff( xz.angle.def )
        xz.angle.mid.class <- 0.5 * ( xz.angle.def[-1] + xz.angle.def[-n] )
        if( 
          n > 2 &&
          identical( xz.angle.def[1], 0. ) && 
          identical( xz.angle.def[n], 180. ) &&
          !all( d[1] == d[-1] )
        ){
          
          xz.angle.mid.class[1] <- xz.angle.mid.class[1] - (180. - xz.angle.mid.class[n-1])
          xz.angle.mid.class <- xz.angle.mid.class[-(n-1)]
        }
                
        if( missing( col ) ) col <- 1:length( xy.angle.mid.class )
        if( missing( pch ) ) pch <- 1:length( xz.angle.mid.class )
        if( missing( xlab ) ) xlab <- "lag distance"
        if( missing( ylab ) ) ylab <- what
        
        if( !add ){
          plot( c(0, xmax), c(0, ymax ), xlab = xlab, 
            ylab = ylab, type = "n", ... )
        }
        
      }
      
      ## add graph of fitted variogram/covariance/correlation model
      
      lines( 
        x, 
        what,
        to = xmax,
        xy.angle = xy.angle.mid.class,
        xz.angle = xz.angle.mid.class,
        col = col, pch = pch, lty = lty, ...
      )
      
    }
    
  )
  
  invisible()
  
}

 ##############################################################################

lines.georob <- lines.fitted.variogram <- 
  function( 
    x, 
    what = c("variogram", "covariance", "correlation"),
    from = 1.e-6, to, n = 501, 
    xy.angle = 90,
    xz.angle = 90,
    col = 1:length( xy.angle ), pch = 1:length( xz.angle ), lty = "solid", ...
  )
{
  
  ## Function plots the graph of a variogram fitted by f.georob
  
  ## 2012-12-12 A. Papritz
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2014-05-08 AP changes for plotting covariances and correlations
  
  d2r <- pi / 180.
  
  what <- match.arg( what )
  
  if( 
    x[["variogram.model"]] %in% control.georob()[["irf.models"]] && 
    what != "variogram" 
  ) stop(
    "stationary covariance and correlation does not exist for intrinsic variogram model 'RMfbm'"  
  )
  
  ## generate grid of angle classes
  
  angle <- expand.grid(
    xy = xy.angle * d2r,
    xz = xz.angle * d2r
  )
  
  ## set up lag distances
  
  if( missing( to ) ){
    u <- par( "usr" )[1:2]
    to <- (u[2] + 0.04/1.04*u[1]) / (1.04 - 0.04^2/1.04)
  }
  
  lag.class <- seq( from, to, length = n )
  
  ## determine number of dimensions
  
  if( identical( class( x ), "fitted.variogram" ) ){
    ndim <- 3
  } else {
    ndim <- NCOL( x[["locations.objects"]][["coordinates"]] )
  }
  
  ## loop over all angles
  
  t.bla <- lapply(
    1:NROW( angle ),
    function( 
      i, what, angle, lag.class, 
      variogram.model, param, aniso, 
      nxy, nxz, ndim, col, pch, lty, ...
    ){
      
      ## generate lag vectors
      
      t.aux <- lag.class * sin( angle[i, "xz"] )
      lag.vector <- cbind(
        t.aux * sin( angle[i, "xy"] ),
        t.aux * cos( angle[i, "xy"] ),
        lag.class * cos( angle[i, "xz"] )
      )
      
      ## drop unneeded components
      
      lag.vector <- lag.vector[, 1:ndim, drop = FALSE ]
      
      ## compute semivariance
      
      r.gamma <- f.aux.gamma(
        lag.vector, variogram.model, param, aniso
      )
      
      if( identical( class( r.gamma ), "try-error" ) || any( is.na( r.gamma ) ) ){
        stop( "\nan error occurred when computing semivariances\n" )
      }
      
      ## plot semivariance
      
      sel.pch <- ((i-1) %/% nxy) + 1
      sel.col <- i - (sel.pch-1) * nxy
      type <- if( nxz > 1 ) "o" else "l"
      
      sill <- sum( param[c("variance", "snugget", "nugget")])
      r.gamma <- switch(
        what,
        variogram = r.gamma,
        covariance = sill - r.gamma,
        correlation = (sill - r.gamma ) / sill
      )

      lines( 
        lag.class, r.gamma, 
        col = col[sel.col], pch = pch[sel.pch], 
        lty = lty, type = type, ... )
    },
    what = what,
    angle = angle,
    lag.class = lag.class,
    variogram.model = x[["variogram.model"]],
    param = x[["param"]],
    aniso = x[["aniso"]],
    nxy = length( xy.angle ),
    nxz = length( xz.angle ),
    ndim = ndim,
    col = col, pch = pch, lty = lty
  )
  
  invisible( NULL )
  
}
