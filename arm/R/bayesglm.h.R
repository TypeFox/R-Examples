## Aug 11, 2007
##   1. model.matrix.bayes, terms.bayes, contr.bayes.unordered
##       & contr.bayes.ordered are in "arm" now.
##   2. bayesglm.h now uses model.matrix.bayes2 in "arm".
#
#bayesglm.h <- function ( formula, family = gaussian, data, weights, subset, 
#    na.action, start = NULL, etastart, mustart, offset, control = glm.control(...), 
#    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, 
#    prior.mean = 0, prior.scale = 2.5, prior.df = 1, scaled = TRUE, 
#    prior.mean.for.intercept = 0, prior.scale.for.intercept = 10, prior.df.for.intercept = 1,
#    batch=0, batch.mean=NA, batch.sd=NA, 
#    batch.mean.mean=0, batch.mean.scale=prior.scale.for.intercept, batch.mean.df=prior.df, 
#    batch.sd.scale=2.5, batch.sd.df=1, 
#    n.iter = 100, drop.baseline = FALSE, separete.intercept = TRUE, 
#    keep.order=TRUE, batch.mean.known=FALSE, ... ) 
#{
#    call <- match.call()
#    if (is.character(family)) 
#        family <- get(family, mode = "function", envir = parent.frame())
#    if (is.function(family)) 
#        family <- family()
#    if (is.null(family$family)) {
#        print(family)
#        stop("'family' not recognized")
#    }
#    if (missing(data)) 
#        data <- environment(formula)
#    mf <- match.call(expand.dots = FALSE)
#    m  <- match(c("formula", "data", "subset", "weights", "na.action", "etastart", "mustart", "offset"), names(mf), 0)
#    mf <- mf[c(1, m)]
#    mf$drop.unused.levels <- TRUE
#    mf[[1]] <- as.name("model.frame")
#    mf <- eval(mf, parent.frame())
#    switch(method, model.frame = return(mf), glm.fit = 1, stop("invalid 'method': ", method))
#    mt <- attr(mf, "terms")
#    Y  <- model.response(mf, "any")
#    if (length(dim(Y)) == 1) {
#        nm <- rownames(Y)
#        dim(Y) <- NULL
#        if (!is.null(nm)) 
#            names(Y) <- nm
#    }
#    if (!drop.baseline){
#        X <- if (!is.empty.model(mt)){ 
#            #class(mt) <- c("bayesglm.h", "terms", "formula")
#            model.matrix.bayes.h( mt, mf, contrasts, keep.order=keep.order, batch=batch ) 
#            }
#        else matrix(, NROW(Y), 0)
#    }
#    else {
#      X <- if (!is.empty.model(mt)) 
#        model.matrix( mt, mf, contrasts )
#        else matrix(, NROW(Y), 0)
#    }
##    if ( length( batch ) == 1 ) { batch <- rep ( batch, ncol( X ) ) }
#    intercept <- (attr(mt, "intercept") > 0)
#    if( intercept && length(batch)==1 ){
#        batch <- c(0,rep (batch, ncol( X )-1))
#    }
#    else if (length(batch)==1 ) {
#        batch <- rep (batch, ncol( X ))
#    }
#    else if ( length( batch ) > 1 ) {
#        if( length( batch ) != (length(attr(mt,"term.labels") ))) {
#            stop( "batch is ether all 0 or must be specified for each of the variables." )
#        }
#        else {
#            assignVec <- attr( X, "assign" )
#            tb <- if ( intercept ) { 0 } else { NULL }
#            for( bi in 1:length( batch ) ){
#                tb<-c( tb,rep( batch[bi], sum( assignVec == bi ) ) )
#            }
#            batch <- tb
#        }
#    }
#
#    weights <- model.weights(mf)
#    offset  <- model.offset(mf)
#    if (!is.null(weights) && any(weights < 0)) 
#        stop("negative weights not allowed")
#    if (!is.null(offset) && length(offset) != NROW(Y)) 
#        stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
#    mustart <- model.extract(mf, "mustart")
#    etastart <- model.extract(mf, "etastart")
#
#    fit <- bayesglm.hierarchical.fit(x = X, y = Y, weights = weights, start = start, 
#            etastart = etastart, mustart = mustart, offset = offset, 
#            family = family, control = glm.control( maxit = n.iter ), 
#            intercept = intercept, prior.mean = prior.mean, 
#            prior.scale = prior.scale, 
#            prior.mean.for.intercept = prior.mean.for.intercept,
#            prior.scale.for.intercept = prior.scale.for.intercept, 
#            prior.df.for.intercept = prior.df.for.intercept,
#            prior.df = prior.df, batch = batch, batch.mean=batch.mean, batch.sd = batch.sd, 
#            batch.mean.mean = batch.mean.mean, batch.mean.scale = batch.mean.scale, batch.mean.df = batch.mean.df, 
#            batch.sd.scale = batch.sd.scale, batch.sd.df = batch.sd.df, scaled = scaled ,drop.baseline=drop.baseline,
#            batch.mean.known = batch.mean.known )
#    if (any(offset) && attr(mt, "intercept") > 0) {
#        cat("bayesglm not yet set up to do deviance comparion here\n")
#        fit$null.deviance <- bayesglm.hierarchical.fit(x = X[, "(Intercept)", drop = FALSE], 
#                y = Y, weights = weights, offset = offset, family = family, 
#                control = control, intercept = intercept, prior.mean = prior.mean, prior.scale = prior.scale, 
#                prior.mean.for.intercept = prior.mean.for.intercept, 
#                prior.scale.for.intercept = prior.scale.for.intercept, 
#                prior.df.for.intercept = prior.df.for.intercept,
#                prior.df = prior.df, batch = batch, batch.mean = batch.mean, batch.sd = batch.sd, 
#                batch.mean.mean = batch.mean.mean, batch.mean.scale = batch.mean.scale, batch.mean.df = batch.mean.df, 
#                batch.sd.scale = batch.sd.scale, batch.sd.df = batch.sd.df, scaled = scaled,drop.baseline=drop.baseline,
#                batch.mean.known = batch.mean.known )$deviance
#    }
#    if (model) 
#        fit$model <- mf
#    fit$na.action <- attr(mf, "na.action")
#    if (x) 
#        fit$x <- X
#    if (!y) 
#        fit$y <- NULL
#    fit <- c(fit, list(call = call, formula = formula, terms = mt, 
#            data = data, offset = offset, control = control, method = method, 
#            contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf)))
#    class(fit) <- c("bayesglm.h","glm", "lm")
#    fit
#}
#
#
#bayesglm.hierarchical.fit <-
#function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
#    mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
#    control = glm.control(), prior.mean = 0, prior.scale = 2.5, prior.df = 1, 
#    intercept = TRUE, 
#    prior.mean.for.intercept = 0, prior.scale.for.intercept = 10, prior.df.for.intercept = prior.df,
#    batch=0, batch.mean=NA, batch.sd=NA,
#    batch.mean.mean=0, batch.mean.scale=2.5, batch.mean.df=1,
#    batch.sd.scale=2.5, batch.sd.df=1, scaled = TRUE, drop.baseline = FALSE, batch.mean.known = TRUE ) 
#{
#    J <- NCOL(x)  
#    if(intercept && length(batch)==1 ){
#        batch <- c(0,rep (batch, J-1))
#    }
#    else if (length(batch)==1 ) {
#        batch <- rep (batch, J)
#    }
#    J.0 <- sum (batch==0)
#    if (J.0 > 0) {
#      if (length(prior.mean) == 1)  { 
#        prior.mean <- rep(prior.mean, J.0)  
#        if(intercept){
#            prior.mean[1] <- prior.mean.for.intercept
#        }
#      }
#      else if (length(prior.mean) > 1) {
#        if( length( prior.mean ) + intercept != J.0 ){
#            stop(message="You must specify the prior.mean for each of the variables")
#        }
#      }
#      if (length(prior.scale) == 1) { 
#        prior.scale <- rep(prior.scale, J.0)
#        if(intercept){
#            prior.scale[1] <- prior.scale.for.intercept
#        }
#      }
#      else if (length(prior.scale) > 1) {
#        if( length( prior.scale ) + intercept != J.0 ){
#            stop(message="You must specify the prior.scale for each of the variables")
#        }
#      }
#      if (scaled == TRUE) {
#        y.scale <- 1
#        if (family$family == "gaussian") {
#            y.obs <- y[!is.na(y)]
#            num.categories <- length(unique(y.obs))
#            if (num.categories == 2) {
#                y.scale <- max(y.obs) - min(y.obs)
#            }
#            else if (num.categories > 2) {
#                y.scale <- 2 * sd(y.obs)
#            }
#        }
#        for (j in 1:J.0) {
#            x.obs <- x[,(1:J)[batch==0][j]]
#            x.obs <- x.obs[!is.na(x.obs)]
#            num.categories <- length(unique(x.obs))
#            x.scale <- 1
#            if (num.categories == 2) {
#                x.scale <- max(x.obs) - min(x.obs)
#            }
#            else if (num.categories > 2) {
#                x.scale <- 2 * sd(x.obs)
#            }
#            prior.scale[j] <- prior.scale[j] * y.scale/x.scale
#        }
#        if (is.numeric(prior.scale.for.intercept) & intercept) {
#            prior.scale[1] <- prior.scale.for.intercept * y.scale
#        }
#      }
#      if (length(prior.df) == 1) {
#        prior.df <- rep(prior.df, J.0)
#      }
#      #### Added by Masanao Yajima 8/30
#      if (intercept){
#        prior.df[1] <- prior.df.for.intercept
#      }
#    }
#    K <- max (batch)
#    if (K > 0){
#      if ( length( batch.mean ) == 1 ) { batch.means <- rep( batch.mean, K ) }
#      if ( length( batch.sd ) == 1 ) { batch.sds <- rep( batch.sd, K ) }
#      if ( length( batch.mean.mean ) == 1 ) { batch.mean.mean <- rep( batch.mean.mean, K ) }
#      if ( length( batch.mean.scale ) == 1 ) { batch.mean.scale <- rep( batch.mean.scale, K ) }
#      if ( length( batch.mean.df ) == 1 ) { batch.mean.df <- rep( batch.mean.df, K ) }
#      if ( length( batch.sd.scale ) == 1 ) { batch.sd.scale <- rep( batch.sd.scale, K ) }
#      if ( length( batch.sd.df ) == 1 ) { batch.sd.df <- rep( batch.sd.df, K ) }
#    }
#    x <- as.matrix( x )
#    xnames <- dimnames( x )[[2]]
#    ynames <- if (is.matrix( y ) ) { rownames( y ) } else { names( y ) }
#    conv <- FALSE
#    nobs <- NROW( y )
#    nvars <- ncol(x)
#    EMPTY <- nvars == 0
#    if ( is.null( weights ) ){ weights<- rep.int( 1, nobs ) }
#    if ( is.null( offset ) ) { offset <- rep.int( 0, nobs ) }
#    variance   <- family$variance
#    dev.resids <- family$dev.resids
#    aic        <- family$aic
#    linkinv    <- family$linkinv
#    mu.eta     <- family$mu.eta
#    if ( !is.function( variance ) || !is.function( linkinv ) ) { stop( "'family' argument seems not to be a valid family object" ) }
#    valideta   <- family$valideta
#    if ( is.null(valideta)){ valideta <- function( eta ) TRUE }
#    validmu    <- family$validmu
#    if ( is.null( validmu ) ) { validmu <- function( mu ) TRUE }
#    if ( is.null( mustart ) ) { eval( family$initialize ) }
#    else {
#        mukeep  <- mustart
#        eval( family$initialize )
#        mustart <- mukeep
#    }
#    if (EMPTY) {
#        eta       <- rep.int( 0, nobs ) + offset
#        if ( !valideta( eta ) ) { stop( "invalid linear predictor values in empty model" ) }
#        mu        <- linkinv( eta )
#        if ( !validmu( mu ) ) { stop( "invalid fitted means in empty model" ) }
#        dev       <- sum( dev.resids( y, mu, weights ) )
#        w         <- ( ( weights * mu.eta( eta )^2 )/variance( mu ) )^0.5
#        residuals <- ( y - mu )/mu.eta( eta )
#        good      <- rep( TRUE, length( residuals ) )
#        boundary  <- conv <- TRUE
#        coef      <- numeric( 0 )
#        iter      <- 0
#    }
#    else {
#        coefold <- NULL
#        eta  <- if (!is.null(etastart)) { etastart }
#                else if ( !is.null( start ) ) {
#                    if ( length( start ) != nvars ) {
#                        stop( gettextf( "length of 'start' should equal %d and correspond to initial coefs for %s", 
#                                nvars, paste( deparse( xnames ), collapse = ", " ) ), domain = NA )
#                    }
#                    else {
#                        coefold <- start
#                        offset + as.vector( if ( NCOL( x ) == 1) { x * start } else { crossprod( x, start ) })
#                        #offset + as.vector( if (NCOL(x) == 1) { x * start } else { x %*% start })
#                    }
#                }
#                else {family$linkfun(mustart)}
#        mu <- linkinv( eta )
#        if ( !( validmu( mu ) && valideta( eta ) ) ) 
#            stop( "cannot find valid starting values: please specify some" )
#        devold     <- sum( dev.resids(y, mu, weights ) )
#        boundary   <- conv <- FALSE
##        prior.sd <- prior.scale
#        dispersion    <- 1
#        dispersionold <- dispersion
#        # Define s's and initialize sigma's
#        mu.0    <- prior.mean
#        s.0     <- prior.scale
#        nu.0    <- prior.df
#        sigma.0 <- s.0
#        # Count the number of batches and record where mu.batch_k and sigma.batch_k are unknown
#        sigma.batch    <- NULL
#        sigma.mu.batch <- NULL
#        if ( K > 0 ) {
#            batch.mean.unknown <- is.na( batch.mean )
#            batch.sd.unknown   <- is.na( batch.sd )
#            # Create the W matrix
#            J.plus <- sum( batch > 0 )
#            W      <- array( 0, c( J, K ) )
#            for ( k in 1:K ){
#                W[batch == k, k] <- 1
#            }
#            W.plus         <- W[batch>0, ]
#            J.batch        <- colSums( W )
#            s.batch        <- batch.sd.scale
#            nu.batch       <- batch.sd.df
#            sigma.batch    <- s.batch
#            mu.mu.batch    <- batch.mean.mean
#            s.mu.batch     <- batch.mean.scale
#            sigma.mu.batch <- s.mu.batch
#            nu.mu.batch    <- ifelse( batch.mean.df == Inf, batch.mean.scale, batch.mean.df )
#            # Prepare the subtotals for the batches with unknown means
#            x.plus         <- x[ ,batch > 0]
#            #x.star         <- rbind( cbind( x, x.plus %*% W.plus ), diag( J+K ) )
#            x.star         <- rbind( cbind( x, tcrossprod( x.plus,t( W.plus ) ) ), diag( J+K ) )
#            if ( intercept ) { x.star[NROW( x )+1, 1:J] <- colMeans( x ) }  # 17 Dec
#            dimnames( x.star )[[2]] <- c ( dimnames( x )[[2]], paste( "mu.batch.", 1:K, sep="" ) )
#            xnames         <- dimnames(x.star)[[2]]  
#        }
#        else {  # if K==0
#            x.star <- as.matrix( rbind( x, diag( J ) ) )
#        }
#        nvars.star <- ncol(x.star)
## Loop #######
#        for ( iter in 1:control$maxit ) {
#            good  <- weights > 0
#            varmu <- variance( mu )[good]
#            if ( any( is.na( varmu ) ) ) { stop( "NAs in V( mu )") }
#            if ( any( varmu == 0 ) ) { stop( "0s in V( mu )" ) }
#            mu.eta.val <- mu.eta( eta )
#            if ( any( is.na( mu.eta.val[good] ) ) ) { stop( "NAs in d( mu )/d( eta )" ) }
#            good <- ( weights > 0 ) & ( mu.eta.val != 0 )
#            if ( all( !good ) ) {
#                conv <- FALSE
#                warning( "no observations informative at iteration ", iter )
#                break
#            }
#            z <- ( eta - offset )[good] + ( y - mu )[good] / mu.eta.val[good]
#            w <- sqrt( ( weights[good] * mu.eta.val[good]^2 ) / variance( mu )[good])
#            ngoodobs <- as.integer( nobs - sum( !good ) )
#            # This is where we augment the data with the prior information
#            if ( K > 0 ){
#              # Added by Masanao Yajima 2007/07/31
#              # when there is batch 0 then 
#              if (min(batch)==0){
#                  z.star <- c( z, mu.0, rep( 0, J.plus ), mu.mu.batch )
#                  w.star <- c( w, sqrt( dispersion )*c( 1/sigma.0, 1/sigma.batch[batch[batch>0]], 1/sigma.mu.batch ) )       
#                  ngoodobs.star <- ngoodobs + NCOL( x ) + NCOL( W.plus )
#              }
#              # when there is no batch 0 then 
#              else{
#                  z.star <- c( z, rep( 0, J.plus ), mu.mu.batch )
#                  w.star <- c( w, sqrt( dispersion ) * c(  1/sigma.batch[batch[batch>0]], 1/sigma.mu.batch ) )       
#                  ngoodobs.star <- ngoodobs + NCOL( x ) + NCOL( W.plus )
#              }
#            }
#            else {
#              z.star <- c( z, mu.0 )
#              w.star <- c( w, sqrt( dispersion )/sigma.0 )
#              ngoodobs.star <- ngoodobs + NCOL( x )
#            }
#            good.star <- c(good, rep( TRUE, J + K ) )
#            nvars     <- NCOL( x.star )
#            if ( intercept ) {
#              x.star[NROW( x ) + 1, 1:NCOL( x )] <- colMeans(x)
#            }
#            fit <- .Fortran( "dqrls", qr = x.star[good.star, ] * w.star, n = ngoodobs.star, 
#                            p = nvars, y = w.star * z.star, ny = as.integer( 1 ), tol = min(1e-07, control$epsilon/1000 ),
#                            coefficients = double( nvars ), residuals = double( ngoodobs.star ), effects = double( ngoodobs.star ), 
#                            rank = integer( 1 ), pivot = 1:nvars, qraux = double( nvars ), work = double( 2 * nvars ), PACKAGE = "base" )
#            if ( any( !is.finite( fit$coefficients ) ) ) {
#                conv <- FALSE
#                warning( "non-finite coefficients at iteration ", iter )
#                break
#            }
#       
#            coefs.hat <- fit$coefficients
#            V.coefs   <- chol2inv( as.matrix(fit$qr)[1:ncol( x.star ), 1:ncol( x.star ), drop = FALSE] )
#            # Now update the prior scale
#            # Allocate the coefficients to beta.0, alpha, mu.batch
#            beta.0.index <- 1:J.0
#            beta.0.hat   <- coefs.hat[beta.0.index]
#            V.beta.0     <- diag(V.coefs)[beta.0.index]
#            # Now update the sigma_j's in batch 0
#            sigma.0      <- ifelse ( nu.0 == Inf, s.0, sqrt( ( ( beta.0.hat - mu.0 )^2 + V.beta.0 + nu.0 * s.0^2 )/( 1 + nu.0 ) ) )
#            if ( K > 0 ) {
#                alpha.index    <- ( J.0 + 1 ):J
#                mu.batch.index <- ( J + 1 ):( J + K )
#                alpha.hat      <- coefs.hat[alpha.index]
#                mu.batch.hat   <- coefs.hat[mu.batch.index]
#                V.alpha        <- diag( V.coefs )[alpha.index] *dispersion  ####
#                V.mu.batch     <- diag( V.coefs )[mu.batch.index]*dispersion  ####
#                # Now estimate the sigma.batch_k's where unknown
#                sigma.batch  <- if( batch.sd.unknown ) { 
#                                    #sqrt( ( t( W.plus ) %*% ( alpha.hat^2 + V.alpha ) + nu.batch * s.batch^2 )/( J.batch + nu.batch ) )
#                                    sqrt( ( crossprod(W.plus,( alpha.hat^2 + V.alpha ) ) + nu.batch * s.batch^2 )/( J.batch + nu.batch ) )
#                                } 
#                                else{ sigma.batch }
#
#
#                # Now estimate the sigma.mu.batch_k's where mu.batch_k's are unknown
#                sigma.mu.batch <-   if ( batch.mean.unknown ) { 
#                                        sqrt( ( ( mu.batch.hat - mu.mu.batch )^2 + V.mu.batch + nu.mu.batch * s.mu.batch^2 )/( 1 + nu.mu.batch ) )
#                                    }
#                                    else{ sigma.mu.batch}
#
#            }
#            start[fit$pivot] <- fit$coefficients      
#            #eta <- drop( as.matrix(x.star[1:nrow( x ), ]) %*% start )
#            eta <- drop( tcrossprod( t(start), as.matrix(x.star[1:nrow( x ), ]) ) )
#            #eta <- drop(x %*% start)
#            mu  <- linkinv( eta <- eta + offset )
#            dev <- sum( dev.resids( y, mu, weights ) )
#            if ( !( family$family %in% c( "poisson", "binomial" ) ) ) {
#                #mse.resid <- mean((w * (z - x %*% coefs.hat))^2)
#                #mse.resid <- mean((w * (z - as.matrix(x.star[1:nrow(x),]) %*% coefs.hat))^2)
#                mse.resid <- mean( ( w * ( z - tcrossprod( as.matrix( x.star[1:nrow(x),] ),t( coefs.hat ) ) ) )^2 )
#                #mse.uncertainty <- mean(diag(x %*% V.coefs %*% t(x))) * dispersion
#                #mse.uncertainty <- mean( rowSums( ( x.star[1:nrow(x),] %*% V.coefs ) * x.star[1:nrow(x),] ) ) * dispersion
#                mse.uncertainty <- mean( rowSums( tcrossprod( x.star[1:nrow(x),], V.coefs ) * x.star[1:nrow(x),] ) ) * dispersion
#                dispersion <- mse.resid + mse.uncertainty
#            }
#            if ( control$trace ) { cat("Deviance =", dev, "Iterations -", iter, "\n") }
#            boundary <- FALSE
#            if ( !is.finite( dev ) ) {
#                if ( is.null( coefold ) ) { 
#                    stop( "no valid set of coefficients has been found: please supply starting values", call. = FALSE ) 
#                }
#                warning( "step size truncated due to divergence", call. = FALSE )
#                ii <- 1
#                while ( !is.finite( dev ) ) {
#                    if ( ii > control$maxit ) { stop( "inner loop 1; cannot correct step size" ) }
#                    ii    <- ii + 1
#                    start <- ( start + coefold )/2
#                    #eta   <- drop( x %*% start )
#                    eta   <- drop( crossprod( x, start) )
#                    mu    <- linkinv( eta <- eta + offset )
#                    dev   <- sum( dev.resids( y, mu, weights ) )
#                }
#                boundary <- TRUE
#                if ( control$trace ){ cat( "Step halved: new deviance =", dev, "\n" ) }
#            }
#            if ( !( valideta( eta ) && validmu( mu ) ) ) {
#                if ( is.null( coefold ) ) {
#                  stop("no valid set of coefficients has been found: please supply starting values", call. = FALSE)
#                }
#                warning("step size truncated: out of bounds", call. = FALSE)
#                ii <- 1
#                while ( !(valideta( eta ) && validmu( mu ) ) ) {
#                  if ( ii > control$maxit ) {
#                    stop("inner loop 2; cannot correct step size")
#                  }
#                  ii    <- ii + 1
#                  start <- ( start + coefold )/2
#                  #eta   <- drop( x %*% start )
#                  eta   <- drop( crossprod(x, start ) )
#                  mu    <- linkinv( eta <- eta + offset )
#                }
#                boundary <- TRUE
#                dev <- sum( dev.resids( y, mu, weights ) )
#                if ( control$trace ) { cat( "Step halved: new deviance =", dev, "\n" ) }
#            }
#            # Convergence Check
#            if (iter > 1 & abs( dev - devold )/( 0.1 + abs( dev ) ) < control$epsilon 
#                & abs( dispersion - dispersionold)/( 0.1 + abs( dispersion ) ) < control$epsilon ) {
#                conv <- TRUE
#                coef <- start
#                break
#            }
#            else {
#                devold <- dev
#                dispersionold <- dispersion
#                coef <- coefold <- start
#            }
#
#        }
## End of Loop #######
#        if ( !conv ) { warning( "algorithm did not converge" ) }
#        if ( boundary ) { warning( "algorithm stopped at boundary value" ) }
#        eps <- 10 * .Machine$double.eps
#        if ( family$family == "binomial" ) {
#            if ( any( mu > 1 - eps ) || any( mu < eps ) ) { warning( "fitted probabilities numerically 0 or 1 occurred" ) }
#        }
#        if ( family$family == "poisson" ) {
#            if ( any(mu < eps ) ) { warning( "fitted rates numerically 0 occurred" ) }
#        }
#        if( drop.baseline==TRUE ){
#            if ( fit$rank < nvars ) {
#                coef[fit$pivot][seq( fit$rank + 1, nvars )] <- NA
#            }
#        }
#        xxnames   <- xnames[fit$pivot]
#        residuals <- rep.int( NA, nobs )
#        residuals[good] <- z - ( eta - offset )[good]
#        fit$qr    <- as.matrix( fit$qr )
#        nr        <- min( sum( good ), nvars )
#        if ( nr < nvars ) {
#            Rmat <- diag( nvars )
#            Rmat[1:nr, 1:nvars] <- fit$qr[1:nr, 1:nvars]
#        }
#        else {
#            Rmat <- fit$qr[1:nvars, 1:nvars]
#        }
#        Rmat <- as.matrix( Rmat )
#        Rmat[ row( Rmat ) > col( Rmat ) ] <- 0
#        names( coef )      <- xnames
#        colnames( fit$qr ) <- xxnames
#        dimnames( Rmat )   <- list( xxnames, xxnames )
#    }
#    names( residuals ) <- ynames
#    names( mu )  <- ynames
#    names( eta ) <- ynames
#    wt         <- rep.int(0, nobs)
#    wt[good]   <- w^2
#    names( wt )  <- ynames
#    names( weights ) <- ynames
#    names( y ) <- ynames
#    wtdmu    <- if ( intercept ) { sum( weights * y )/sum( weights )} else linkinv( offset )
#    nulldev <- sum( dev.resids( y, wtdmu, weights ) )
#    n.ok    <- nobs - sum( weights == 0 )
#    nulldf  <- n.ok - as.integer( intercept )
#    rank    <- if ( EMPTY ) { 0 } else { fit$rank }
#    resdf   <- n.ok - rank
#    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
#    list( coefficients = coef, residuals = residuals, fitted.values = mu, 
#        effects = if ( !EMPTY ) fit$effects, R = if ( !EMPTY ) Rmat, rank = rank, 
#        qr = if ( !EMPTY ) structure(fit[c( "qr", "rank", "qraux", "pivot", "tol" )], class = "qr" ), 
#        family = family, linear.predictors = eta, deviance = dev, aic = aic.model, 
#        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
#        df.residual = resdf, df.null = nulldf, y = y, converged = conv, boundary = boundary, 
#        prior.mean = prior.mean, prior.scale = prior.scale, 
#        prior.df = prior.df, prior.sd = sigma.0, dispersion = dispersion, 
#        batch=batch, batch.mean=batch.mean, batch.sd=batch.sd,
#        batch.mean.mean=batch.mean.mean, batch.mean.scale=batch.mean.scale, batch.mean.df =batch.mean.df, 
#        batch.sd.scale=batch.sd.scale, batch.sd.df=batch.sd.df, 
#        sigma.0=sigma.0, sigma.batch=sigma.batch, sigma.mu.batch=sigma.mu.batch )
#}
#
#setMethod("print", signature(x = "bayesglm.h"), 
#  function(x, digits=2) display(object=x, digits=2))
#setMethod("show", signature(object = "bayesglm.h"), 
#    function(object) display(object, digits=2))
