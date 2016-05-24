georob <- 
  function( 
    formula, data, subset, weights, na.action,
    model = TRUE, x = FALSE, y = FALSE, 
    contrasts = NULL, offset, 
    locations,
    variogram.model = c( "RMexp", "RMaskey", "RMbessel", "RMcauchy", 
      "RMcircular", "RMcubic", "RMdagum", "RMdampedcos", "RMdewijsian", "RMfbm",
      "RMgauss", "RMgencauchy", "RMgenfbm", "RMgengneiting", "RMgneiting", "RMlgd",
      "RMmatern", "RMpenta", "RMqexp", "RMspheric", "RMstable",
      "RMwave", "RMwhittle"
    ), 
    param, fit.param = default.fit.param()[names(param)],
	aniso = default.aniso(), fit.aniso = default.fit.aniso(),
    tuning.psi = 2, control = control.georob(), verbose = 0,
    ...
  )
{
  
  ## wrapper function to georob.fit with "standard" interface for
  ## statistical modelling
  
  ## ToDos:
  
  ## - Variante fuer Klasse SpatialPointsDataFrame{sp}
  ## - Ueberpruefung der Konsistenz des Inputs hier statt in georob.fit
  ## - ausgewaehlte Elemente von control in Resultate speichern
  
  
  ## History:
  
  ## georob:
  
  ## - Implementierung von Schaetzung von Startwerten der Kovarianzparameter mittels REML
  ##   von um Ausreisser bereinigten Datensatz
  ## - Implementierung von Schaetzung von Nugget und Varianz von mikroskaliger Variation (snugget )
  ## - Implementierung der Schaetzungen fuer replizierte Messungen an gleicher Messtelle    
  ## - neue Varianten um Startwerte von betahat zu rechnen (Ersatz von use.lm durch 
  ##   Argument initial.fixef = c( "rq", "lmrob", "rlm" ), vgl. control.georob)
  ## - Kontrolle, ob Designmatrix vollen Spaltenrang hat
  ## - Modifikation fuer Fall, dass alle Variogrammparameter fixiert sind
  ## - IRWLS Berechung von betahat und bhat entweder von Werten im initial.object 
  ##   oder von Schaetzwerten aus vorangehender Iteration
  ## - Berechung der Kovarianzen zwischen betahat und (B - bhat)
  ## - Steuerung der Berechnung der verschiedenen Kovarianzen via control.georob
  ## - vollstaendige Implementierung von Standard Interfaces fuer Input und
  ##   Output fuer statistische Modellierung (analog lm, lmrob)
  ## - korrekte Behandlung von NAs und Implementierung von subset
  ## - Verzicht auf Berechnung von initialem Variogramm
  ## - Startwerte fuer bhat alle gleich Null gesetzt
  ## - neue Struktur von initial.objects
  ## - teilweise geaenderte Namen der Argumente (model -> variogram.model)
  
  ## f.glsrob803:
  
  ## - Umbenennung einiger Argument
  
  ## f.glsrob801:
  
  ## - teilweise Implementation von Standard Interfaces fuer Input und
  ##   Output fuer statistische Modellierung 
  ## - Berechnung des Objekts mit den Startwerten
  
  
  ## 2012-04-21 AP
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-05-07 AP correction of error for constant trend
  ## 2012-05-28 AP handle missing names of coefficients after calling update
  ## 2013-04-23 AP new names for robustness weights
  ## 2013-05-23 AP correct handling of missing observations and to construct model.frame
  ## 2013-06-03 AP handling design matrices with rank < ncol(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-02 AP new transformation of rotation angles
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2013-09-06 AP exclusive use of nleqslv for solving estimating equations
  ## 2014-02-18 AP correcting error when fitting models with offset
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2014-05-22 AP correcting error when selecting initial.param
  ## 2014-06-02 AP partial matching of names of variogram parameters
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  ## 2014-08-26 AP changes to return ml.method if fitted object
  ## 2015-06-30 AP function and arguments renamed
  ## 2015-07-17 AP change of control of computation of hessian for Gaussian (RE)ML, 
  ##               changes for singular design matrices, nlminb optimizer added
  ## 2015-08-19 AP control about error families for computing covariances added
  ## 2015-08-28 AP computation of hessian suppressed; correction of error when using georob.object;
  ##               control arguments hessian, initial.param, initial.fixef newly organized
  ## 2015-11-25 AP new way to control which variogram parameters are fitted
  
  ## check validity of tuning.psi
  
  if( identical( control[["psi.func"]], "t.dist" ) && tuning.psi <= 1. ) 
    stop( "'tuning.psi' must be greater than 1 for t-dist psi-function" )
    
  ## check whether all mandatory arguments have been provided
  
  if( missing( formula ) || missing( locations ) || missing( param ) ) stop( 
    "some mandatory arguments are missing" 
  )
  
  ## check whether anisotropy parameters were passed to georob
  
  aniso.missing <- missing( aniso ) && missing( fit.aniso )
  
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
  
  ## get model frame, response vector, weights, offset and design
  ## matrix (cf.  lm, lmrob)
  
  ret.x <- x
  ret.y <- y
  
#   ## vector with row number of included observations
#   
#   in.subset <- 1:NROW( data )
#   if( !missing( subset ) ) in.subset <- in.subset[subset]
  
  ## build combined formula for fixed effects and locations
  
  extended.formula <- update( 
    formula,
    paste( 
      paste( as.character( formula )[c(2, 1, 3)], collapse = " " ),
      as.character( locations )[2], sep = " + "
    )
  )
  
  ## setting-up model frame
  
  cl <- match.call()
  mf <- match.call( expand.dots = FALSE )
  m <- match( 
    c( "formula", "data", "subset", "weights", "na.action", "offset"),
    names(mf), 0L 
  )
  mf <- mf[c(1L, m)]
  mf[["formula"]] <- extended.formula
  mf[["drop.unused.levels"]] <- TRUE
  mf[[1L]] <- as.name( "model.frame" )
  
  mf <- eval( mf, parent.frame() )
  
  ## setting-up terms objects
  
  mt     <- terms( formula )
  mt.loc <- terms( locations )
  
  ## eliminate intercept from mt.loc
  
  attr( mt.loc, "intercept" ) <- 0
    
#   ## ... and assign fixed effects terms object as attribute to model.frame
#   
#   attr( mf, "terms" ) <- mt
  
  ## check whether 'empty' models have been entered
  
  if( is.empty.model( mt ) )
    stop( "an 'empty' fixed effects model has been specified" )
  if( is.empty.model( mt.loc ) )
    stop( "an 'empty' locations model has been specified" )
  
  ## check whether fixed effects model includes an intercept if an
  ## intrinsic variogram model is used
  
  if( identical( attr( mt, "intercept" ), 0L ) && 
    variogram.model %in% control[["irf.model"]] )
  stop(
    "the fixed effects model must include an intercept ",
    "if an unbounded variogram model is used"
  )
  
  ## extract fixed effects response variable, weights, offset and design matrix
  
  y <- model.response( mf, "numeric" )
  
  w <- as.vector( model.weights( mf ) )
  if( !is.null(w) )
    stop( "weights are not yet implemented for this estimator" )
  
  offset <- as.vector( model.offset(mf) )
  if( !is.null(offset) ) {
    if( length( offset ) != NROW(y) )
      stop( gettextf(
        "number of offsets is %d, should equal %d (number of observations)", 
        length(offset), NROW(y) ), domain = NA )
  }
  
  x <- model.matrix( mt, mf, contrasts )
  
  ## check if optionally provided bhat has correct length
  
  if( !is.null( control[["bhat"]] ) && length( y ) != length( control[["bhat"]] ) ) stop(
    "lengths of response vector and 'bhat' do not match"    
  )
  
  ## adjust initial.param if all variogram parameters are fixed
  
  if( !any( c( fit.param, fit.aniso ) ) ) control[["initial.param"]] <- FALSE
    
  ## adjust choice for initial.fixef to compute regression coefficients
  
  if( tuning.psi < control[["tuning.psi.nr"]] ){
    if( control[["initial.param"]] ) control[["initial.fixef"]] <- "lmrob"
  } else {
    control[["initial.param"]] <- FALSE
  }
  
  ## check whether design matrix has full column rank
  
  col.rank.XX <- list( deficient = FALSE, rank = NCOL( x ) )
  
  sv <- svd( crossprod( x ) )[["d"]]
  min.max.sv <- range( sv )
  condnum <- min.max.sv[1] / min.max.sv[2] 
  
  if( condnum <= control[["min.condnum"]] ){
    col.rank.XX[["deficient"]] <- TRUE
    col.rank.XX[["col.rank"]] <- sum( sv / min.max.sv[2] > control[["min.condnum"]] )
    cat( 
      "design matrix has not full column rank (condition number of X^T X: ", 
      signif( condnum, 2 ), ")\ninitial values of fixed effects coefficients are computed by 'lm'\n\n"
    )
    control[["initial.fixef"]] <- "lm"
    control[["initial.param"]]  <- FALSE
    warning( 
      "design matrix has not full column rank (condition number of X^T X: ", 
      signif( condnum, 2 ), ")\ninitial values of fixed effects coefficients are computed by 'lm'"
    )
  }
  
  ## subtract offset
  
  yy <- y
  if( !is.null( offset ) ) yy <- yy - offset
  
  ## compute initial guess of fixed effects parameters (betahat)
  
  r.initial.fit <- switch(
    control[["initial.fixef"]],
    rq = {
      
      Rho <- function( u, tau) u * (tau - (u < 0))
      tau <- control[["rq"]][["tau"]]
      process <- (tau < 0 || tau > 1)
      
      f.rq.fit <- rq.fit
      formals( f.rq.fit ) <- c( alist( x=, y= ), control[["rq"]] )
      fit <- f.rq.fit( x = x, y = yy ) 
      
      if( process ) {
        rho <- list(x = fit[["sol"]][1, ], y = fit[["sol"]][3, ])
      } else {
        dimnames(fit[["residuals"]]) <- list( dimnames( x )[[1]], NULL )
        rho <- sum( Rho( fit[["residuals"]], tau ) )
      }
      if( control[["rq"]][["method"]] == "lasso" ){
        class(fit) <- c("lassorq", "rq")
      } else if( control[["rq"]][["method"]] == "scad"){
        class(fit) <- c("scadrq", "rq")
      } else {
        class(fit) <- ifelse(process, "rq.process", "rq")
      }
      fit[["na.action"]] <- attr( mf, "na.action" )
      fit[["formula"]] <- formula
      fit[["terms"]] <- mt
      fit[["xlevels"]] <- .getXlevels(mt, mf)
      fit[["call"]] <- cl
      fit[["tau"]] <- tau
      fit[["weights"]] <- w
      fit[["residuals"]] <- drop( fit[["residuals"]] )
      fit[["rho"]] <- rho
      fit[["method"]] <- control[["rq"]][["method"]]
      fit[["fitted.values"]] <- drop( fit[["fitted.values"]] )
      attr(fit, "na.message") <- attr( m, "na.message" )
      if( model ) fit[["model"]] <- mf
      fit
      
    },
    lmrob = {
      
      fit <- lmrob.fit( x, yy, control = control[["lmrob"]] )
      fit[["na.action"]] <- attr(mf, "na.action")
      fit[["offset"]] <- offset
      fit[["contrasts"]] <- attr(x, "contrasts")
      fit[["xlevels"]] <- .getXlevels(mt, mf)
      fit[["call"]] <- cl
      fit[["terms"]] <- mt
      if( control[["lmrob"]][["compute.rd"]] && !is.null(x) )
      fit[["MD"]] <- robMD( x, attr( mt, "intercept" ) )
      if( !is.null( offset ) ) fit[["fitted.values"]] + offset
      fit
      
    },
    lm = {
      
      fit <- if( is.null(w) ){
        lm.fit(x, y, offset = offset, singular.ok = TRUE )
      } else {
        lm.wfit(x, y, w, offset = offset, singular.ok = TRUE )
      }
      class(fit) <- c(if (is.matrix(y)) "mlm", "lm")
      fit[["na.action"]] <- attr(mf, "na.action")
      fit[["offset"]] <- offset
      fit[["contrasts"]] <- attr(x, "contrasts")
      fit[["xlevels"]] <- .getXlevels(mt, mf)
      fit[["call"]] <- cl
      fit[["terms"]] <- mt
      if (model) fit[["model"]] <- mf
      if (ret.x) fit[["x"]] <- x
      if (ret.y) fit[["y"]] <- y
      fit[["qr"]] <- NULL
      fit
      
    }
  )
  
  ## match variogram model
  
  variogram.model <- match.arg( variogram.model )  
  
  ## compute coordinates of locations and distance object
  
  locations.coords <- model.matrix( mt.loc, mf )
  
  if( 
    !( missing( aniso ) || missing( fit.aniso ) ) && 
    ( NCOL( locations.coords ) < 2 || NCOL( locations.coords ) > 3 )
  ) stop( 
    "anisotropic variogram models are implemented only for 2 or 3 dimensions" 
  )

  names( yy ) <- rownames( mf )
  
  ##  check whether argument "object."  has been provided in call (e.g. by
  ##  update ) and extract'locations' exists in workspace
  
  extras <- match.call( expand.dots = FALSE )$...
  georob.object <- extras[names(extras) %in% "object."]
  if( 
    length( georob.object ) && 
    exists( as.character( georob.object ), envir = parent.frame() ) #&&
#     !control[["initial.param"]]
  ){
    if( verbose > 4 ) cat(
      "\n    georob: using some components of 'object.'\n"    
    )
    georob.object <- eval( 
      georob.object[[1]], parent.frame() 
    )[c( "param", "aniso", "locations.objects", "Valphaxi.objects" )]
  } else {
    georob.object <- NULL
  }
  
  #   ##  create environment to store items required to compute likelihood and
  #   ##  estimating equations that are provided by
  #   ##  likelihood.calculations
  #   
  #   envir <- new.env()
  #   lik.item <- list()
  #   assign( "lik.item", lik.item, pos = as.environment( envir ) )
  
  ## root.finding <- match.arg( root.finding )
  
  ## compute initial values of variogram and anisotropy parameters
  
  if( tuning.psi < control[["tuning.psi.nr"]] ){
        
    if( control[["initial.param"]] ){
      
      if( verbose > 0 ) cat( "\ncomputing robust initial parameter estimates ...\n" )
      
      t.sel <- r.initial.fit[["rweights"]] > control[["min.rweight"]]
      
      if( verbose > 0 ) cat( 
        "\ndiscarding", sum( !t.sel ), "of", length( t.sel ), 
        "observations for computing initial estimates of variogram\nand anisotropy parameter by gaussian reml\n"
      )
      
      ## collect.items for initial object
      
      initial.objects <- list(
        x = as.matrix( x[t.sel, ] ),
        y = yy[t.sel],
        betahat = coef( r.initial.fit ),
        bhat = if( is.null( control[["bhat"]] ) ){
          rep( 0., length( yy ) )[t.sel]
        } else {
          control[["bhat"]][t.sel]
        },
        initial.fit = r.initial.fit,
        locations.objects = list(
          locations = locations,
          coordinates = locations.coords[t.sel, ]
        ),
        isotropic = aniso.missing
      )
            
      ## estimate model parameters with pruned data set
      
      t.georob <- georob.fit(
        ## root.finding = root.finding,
        #         slv = TRUE,
        #         envir = envir,
        initial.objects = initial.objects,
        variogram.model = variogram.model, param = param, fit.param = fit.param,
        aniso = aniso, fit.aniso = fit.aniso,
        param.tf = control[["param.tf"]],
        fwd.tf = control[["fwd.tf"]], 
        deriv.fwd.tf = control[["deriv.fwd.tf"]],
        bwd.tf = control[["bwd.tf"]], 
        georob.object = georob.object,
        safe.param = control[["safe.param"]],
        tuning.psi = control[["tuning.psi.nr"]],
        error.family.estimation = control[["error.family.estimation"]],
        error.family.cov.effects = control[["error.family.cov.effects"]],
        error.family.cov.residuals = control[["error.family.cov.residuals"]],
        cov.bhat = control[["cov.bhat"]], full.cov.bhat = control[["full.cov.bhat"]],
        cov.betahat = control[["cov.betahat"]], 
        cov.bhat.betahat = control[["cov.bhat.betahat"]],
        cov.delta.bhat = control[["cov.delta.bhat"]],
        full.cov.delta.bhat = control[["full.cov.delta.bhat"]],
        cov.delta.bhat.betahat = control[["cov.delta.bhat.betahat"]],
        cov.ehat = control[["cov.ehat"]], full.cov.ehat = control[["full.cov.ehat"]],
        cov.ehat.p.bhat = control[["cov.ehat.p.bhat"]], full.cov.ehat.p.bhat = control[["full.cov.ehat.p.bhat"]],
        aux.cov.pred.target = control[["aux.cov.pred.target"]],
        min.condnum = control[["min.condnum"]],
        col.rank.XX = col.rank.XX,
        psi.func = control[["psi.func"]],
        tuning.psi.nr = tuning.psi,
        ml.method = control[["ml.method"]],
        maximizer = control[["maximizer"]],
        reparam = control[["reparam"]],
        irwls.initial = control[["irwls.initial"]],
        irwls.maxiter = control[["irwls.maxiter"]],
        irwls.ftol = control[["irwls.ftol"]],
        force.gradient = control[["force.gradient"]],
        zero.dist = control[["zero.dist"]],
        control.nleqslv =  control[["nleqslv"]],
        control.optim = control[["optim"]],
        control.nlminb = control[["nlminb"]],
        hessian = FALSE,
        control.pmm = control[["pmm"]],
        verbose = verbose
      )
      
      param = t.georob[["param"]][names(fit.param)]
      aniso = t.georob[["aniso"]][["aniso"]][names(fit.aniso)]
      
    } 
    
  }
  
  ## collect.items for initial object
  
  initial.objects <- list(
    x = as.matrix( x ),
    y = yy,
    betahat = coef( r.initial.fit ),
    bhat = if( is.null( control[["bhat"]] ) ){
      rep( 0., length( yy ) )
    } else {
      control[["bhat"]]
    },
    initial.fit = r.initial.fit,
    locations.objects = list(
      locations = locations,
      coordinates = locations.coords
    ),
    isotropic = aniso.missing
  )
  
  ## estimate model parameters

  if( verbose > 0 ) cat( "\ncomputing final parameter estimates ...\n" )

  r.georob <- georob.fit(
    ## root.finding = root.finding,
    #     slv = TRUE,
    #     envir = envir,
    initial.objects = initial.objects,
    variogram.model = variogram.model, param = param, fit.param = fit.param,
    aniso = aniso, fit.aniso = fit.aniso,
    param.tf = control[["param.tf"]],
    fwd.tf = control[["fwd.tf"]], 
    deriv.fwd.tf = control[["deriv.fwd.tf"]],
    bwd.tf = control[["bwd.tf"]], 
    georob.object = georob.object,
    safe.param = control[["safe.param"]],
    tuning.psi = tuning.psi,
    error.family.estimation = control[["error.family.estimation"]],
    error.family.cov.effects = control[["error.family.cov.effects"]],
    error.family.cov.residuals = control[["error.family.cov.residuals"]],
    cov.bhat = control[["cov.bhat"]], full.cov.bhat = control[["full.cov.bhat"]],
    cov.betahat = control[["cov.betahat"]], 
    cov.bhat.betahat = control[["cov.bhat.betahat"]],
    cov.delta.bhat = control[["cov.delta.bhat"]],
    full.cov.delta.bhat = control[["full.cov.delta.bhat"]],
    cov.delta.bhat.betahat = control[["cov.delta.bhat.betahat"]],
    cov.ehat = control[["cov.ehat"]], full.cov.ehat = control[["full.cov.ehat"]],
    cov.ehat.p.bhat = control[["cov.ehat.p.bhat"]], full.cov.ehat.p.bhat = control[["full.cov.ehat.p.bhat"]],
    aux.cov.pred.target = control[["aux.cov.pred.target"]],
    min.condnum = control[["min.condnum"]],
    col.rank.XX = col.rank.XX,
    psi.func = control[["psi.func"]],
    tuning.psi.nr = control[["tuning.psi.nr"]],
    ml.method = control[["ml.method"]],
    maximizer = control[["maximizer"]],
    reparam = control[["reparam"]],
    irwls.initial = control[["irwls.initial"]],
    irwls.maxiter = control[["irwls.maxiter"]],
    irwls.ftol = control[["irwls.ftol"]],
    force.gradient = control[["force.gradient"]],
    zero.dist = control[["zero.dist"]],
    control.nleqslv =  control[["nleqslv"]],
    control.optim = control[["optim"]],
    control.nlminb = control[["nlminb"]],
    hessian = control[["hessian"]],
    control.pmm = control[["pmm"]],
    verbose = verbose
  )
  
  ## add offset to fitted values
  
  if( !is.null( offset ) )
    r.georob[["fitted.values"]] <- r.georob[["fitted.values"]] + offset
    
  ##
    
  r.georob[["control"]] <- control
  
  ## add remaining items to output
  
  if( control[["lmrob"]][["compute.rd"]] && !is.null( x ) )
    r.georob[["MD"]] <- robMD( x, attr(mt, "intercept") )

  if( model ) r.georob[["model"]] <- mf
  if( ret.x ) r.georob[["x"]] <- x
  if( ret.y ) r.georob[["y"]] <- y

  r.georob[["df.residual"]] <- length(yy) - col.rank.XX[["rank"]]
  r.georob[["na.action"]] <- attr(mf, "na.action")
  r.georob[["offset"]] <- offset
  r.georob[["contrasts"]] <- attr(x, "contrasts")
  r.georob[["xlevels"]] <- .getXlevels(mt, mf)
  r.georob[["rank"]] <- col.rank.XX[["rank"]]
  r.georob[["call"]] <- cl
  r.georob[["terms"]] <- mt
  
  ## set missing names of coefficients (bug of update)
  
  if( length( r.georob[["coefficients"]] ) == 1 && is.null( names( r.georob[["coefficients"]] ) ) ){
    names( r.georob[["coefficients"]] ) <- "(Intercept)"
  }

  
  class( r.georob ) <- c( "georob" )
  
  invisible( r.georob )
  
}


##  ##############################################################################

pmm <- 
  function( 
    A, B, control = control.pmm() 
  )
{
  
  ## function for parallelized matrix multiplication inspired by function
  ## parMM{snow}
  
  ## 2014-06-25 A. Papritz
  ## 2015-03-13 AP small changes in f.aux
  
  ## auxiliary function 
  
  f.aux <- function(i, s, e, A, B ) A %*% B[ , s[i]:e[i], drop = FALSE]
  
  ## determine columns indices of matrix blocks
  
  k <- control[["f"]] * control[["ncores"]]
  n <- NCOL(B)
  dn <- floor( n / k )
  s <- ( (0:(k-1)) * dn ) + 1
  e <- (1:k) * dn
  e[k] <- n
  
  ##
  
  if( control[["ncores"]] > 1L ){
    
    if( identical( .Platform[["OS.type"]], "windows") ){
      
      ## create a SNOW cluster on windows OS
      
      if( !sfIsRunning() ){
        options( error = f.stop.cluster )
        junk <- sfInit( parallel = TRUE, cpus = control[["ncores"]] )
      }
      
      res <- sfLapply( 1:k, f.aux, s = s, e = e, A = A, B = B )
      
      if( control[["sfstop"]] ){
        junk <- sfStop()
        options( error = NULL )
      }
      
    } else {
      
      res <- mclapply( 
        1:k, f.aux,s = s, e = e, A = A, B = B,  mc.cores = control[["ncores"]] 
      )
      
    }
    
    matrix( unlist(res), nrow = nrow( A ) )
    
  } else A %*% B    
  
}

#  ##############################################################################

control.georob <- 
  function(
    ml.method = c( "REML", "ML" ),
    reparam = TRUE,
    maximizer = c( "nlminb", "optim" ),
    initial.param = TRUE,
    initial.fixef = c("lmrob", "rq", "lm"),
    bhat = NULL,
    min.rweight = 0.25,
    param.tf = param.transf(),
    fwd.tf = fwd.transf(), 
    deriv.fwd.tf = dfwd.transf(), 
    bwd.tf = bwd.transf(),
    safe.param = 1.e12,
    psi.func = c( "logistic", "t.dist", "huber" ),
    tuning.psi.nr = 1000,
    irwls.initial = TRUE,
    irwls.maxiter = 50, irwls.ftol = 1.e-5,
    force.gradient = FALSE,
    min.condnum = 1.e-12,
    zero.dist = sqrt( .Machine[["double.eps"]] ),
    error.family.estimation    = c( "gaussian", "long.tailed" ),
    error.family.cov.effects   = c( "gaussian", "long.tailed" ),
    error.family.cov.residuals = c( "long.tailed", "gaussian" ),
    cov.bhat = FALSE, full.cov.bhat = FALSE,
    cov.betahat = TRUE, 
    cov.bhat.betahat = FALSE,
    cov.delta.bhat = TRUE, full.cov.delta.bhat = TRUE,
    cov.delta.bhat.betahat = TRUE,
    cov.ehat = TRUE, full.cov.ehat = FALSE,
    cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
    aux.cov.pred.target = FALSE,
    hessian = TRUE,
    rq = control.rq(),
    lmrob = lmrob.control(),
    nleqslv = control.nleqslv(),
    optim = control.optim(),
    nlminb = control.nlminb(),
    pmm = control.pmm(),
    ...
  )
{
  
  ## auxiliary function to set meaningful default values for georob

  ## Arguments: 
  
  ## 2012-04-21 A. Papritz   
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-05-04 AP modifications for lognormal block kriging
  ## 2013-04-23 AP new names for robustness weights
  ## 2013-06-12 AP changes in stored items of Valphaxi object
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-12 AP solving estimating equations by BBsolve{BB} (in addition to nleqlsv)
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2014-08-18 AP changes for parallelized computations
  ## 2014-08-18 AP changes for Gaussian ML estimation
  ## 2015-03-10 AP extended variogram parameter transformations
  ## 2015-06-30 AP function and arguments renamed
  ## 2015-07-17 AP nlminb optimizer added
  ## 2015-08-19 AP control about error families for computing covariances added
  ## 2015-08-28 AP control arguments hessian, initial.param, initial.fixef newly organized
  
  if( 
    !( all( unlist( param.tf ) %in% names( fwd.tf ) ) &&
       all( unlist( param.tf ) %in% names( deriv.fwd.tf ) ) &&
       all( unlist( param.tf ) %in% names( bwd.tf ) )
    )
  ) stop( 
    "undefined transformation of variogram parameters; extend respective function definitions" 
  )
  
  if( !irwls.initial && irwls.ftol >= 1.e-6 ) warning( 
    "'irwls.initial == FALSE' and large 'ftol' may create problems for root finding"
  )
  
  list(
    ml.method = match.arg( ml.method ), reparam = reparam,
    maximizer = match.arg( maximizer ),
    initial.param = initial.param,
    initial.fixef = match.arg( initial.fixef ),
    bhat = bhat,
    min.rweight = min.rweight,
    param.tf = param.tf, fwd.tf = fwd.tf, deriv.fwd.tf = deriv.fwd.tf, bwd.tf = bwd.tf,
    safe.param = safe.param,
    psi.func = match.arg( psi.func ), 
    tuning.psi.nr = tuning.psi.nr,
    irwls.initial = irwls.initial, irwls.maxiter = irwls.maxiter, irwls.ftol = irwls.ftol,
    force.gradient = force.gradient,
    min.condnum = min.condnum,
    zero.dist = zero.dist,
    error.family.estimation = match.arg( error.family.estimation ),
    error.family.cov.effects = match.arg( error.family.cov.effects ),
    error.family.cov.residuals = match.arg( error.family.cov.residuals ),    
    cov.bhat = cov.bhat, full.cov.bhat = full.cov.bhat,
    cov.betahat = cov.betahat, 
    cov.bhat.betahat = cov.bhat.betahat,
    cov.delta.bhat = cov.delta.bhat, full.cov.delta.bhat = full.cov.delta.bhat,
    cov.delta.bhat.betahat = cov.delta.bhat.betahat,
    cov.ehat = cov.ehat, full.cov.ehat = full.cov.ehat,
    cov.ehat.p.bhat = cov.ehat.p.bhat, full.cov.ehat.p.bhat = full.cov.ehat.p.bhat,
    aux.cov.pred.target = aux.cov.pred.target,
    hessian = hessian,
    irf.models = c( "RMdewijsian", "RMfbm", "RMgenfbm" ),
    rq = rq, lmrob = lmrob, nleqslv = nleqslv, 
    optim = optim, 
    nlminb = nlminb,
    pmm = pmm
  )
  
}

## ======================================================================
param.transf <-
  function( 
    variance = "log", snugget = "log", nugget = "log", scale = "log", 
    alpha = c( 
      RMaskey = "log", RMdewijsian = "logit2", RMfbm = "logit2", RMgencauchy = "logit2", 
      RMgenfbm = "logit2", RMlgd = "identity", RMqexp = "logit1", RMstable = "logit2"
    ), 
    beta = c( RMdagum = "logit1", RMgencauchy = "log", RMlgd = "log" ), 
    delta = "logit1", 
    gamma = c( RMcauchy = "log", RMdagum = "logit1" ), 
    kappa = "logit3", lambda = "log", 
    mu = "log", 
    nu = "log",
    f1 = "log", f2  ="log", omega = "identity", phi = "identity", zeta = "identity"
  )
{
  
  ## function sets meaningful defaults for transformation of variogram
  ## parameters
  
  ## 2013-07-02 A. Papritz
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2015-03-10 AP extended transformation
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  
  list( 
    variance = variance, snugget = snugget, nugget = nugget, scale = scale,
    alpha = alpha, 
    beta = beta, 
    delta = delta, 
    gamma = gamma, 
    kappa = kappa, 
    lambda = lambda, 
    mu = mu, 
    nu = nu, 
    f1 = f1, f2 = f2, omega = omega, phi = phi, zeta = zeta
  )
  
}

## ======================================================================
fwd.transf <- 
  function( 
    ...
  )
{
  
  ## definition of forward transformation of variogram parameters
  
  ## 2013-07-02 A. Papritz
  ## 2015-03-10 AP extended variogram parameter transformations
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  
  list( 
    log = function(x) log(x), 
    logit1 = function(x) log( x / (1. - x) ), 
    logit2 = function(x) log( x / (2. - x) ), 
    logit3 = function(x) log( (x - 1.) / (3. - x) ), 
    identity = function(x) x, ...
  )
}

## ======================================================================
dfwd.transf<-
  function( 
    ...
  )
{
  
  ## definition of first derivative of forward transformation of variogram
  ## parameters
  
  ## 2013-07-02 A. Papritz
  ## 2015-03-10 AP extended variogram parameter transformations
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  
  list( 
    log = function(x) 1/x, 
    logit1 = function(x) 1. / (x - x^2),
    logit2 = function(x) 2. / (2.*x - x^2),
    logit3 = function(x) 2. / (4.*x - 3. - x^2),
    identity = function(x) rep(1., length(x)), ... 
  )  
  
}

## ======================================================================
bwd.transf <-
  function( 
    ...
  )
{
  
  ## definition of backward transformation of variogram parameters
  
  ## 2013-07-02 A. Papritz
  ## 2015-03-10 AP extended variogram parameter transformations
  ## 2015-04-07 AP changes for fitting anisotropic variograms
  
  list( 
    log = function(x) exp(x), 
    logit1 = function(x) exp(x) / (1. + exp(x)),
    logit2 = function(x) 2. * exp(x) / (1. + exp(x)),
    logit3 = function(x) (3. * exp(x) + 1.) / (1. + exp(x)),
    identity = function(x) x, ...
  )
  
}

## ======================================================================
control.rq <-
  function(
    ## arguments for rq
    tau = 0.5, rq.method = "br",
    ## arguments for rq.fit.br
    rq.alpha = 0.1, ci = FALSE, iid = TRUE, interp = TRUE, tcrit = TRUE,
    ## arguments for rq.fit.fnb
    rq.beta = 0.99995, eps = 1e-06,
    ## arguments for rq.fit.pfn
    Mm.factor = 0.8, max.bad.fixup = 3,
    ...
  )
{
  
  ## function sets meaningful defaults for selected arguments of function
  ## rq{quantreg}
  
  ## 2012-12-14 A. Papritz
  ## 2014-07-29 AP
  ## 2015-06-30 AP function and arguments renamed
  
  list( 
    tau = tau,  method = rq.method,
    alpha = rq.alpha, ci = ci, iid = iid, interp = interp, tcrit = tcrit,
    beta = rq.beta, eps = eps,
    Mm.factor = Mm.factor, max.bad.fixup = max.bad.fixup
  )
}


## ======================================================================
control.nleqslv <-
  function( 
    method = c( "Broyden", "Newton"), 
    global = c( "dbldog", "pwldog", "qline", "gline", "none" ),
    xscalm = c( "fixed", "auto" ),
    control = list( ftol = 1.e-4 ), 
    ...
  )
{
  
  ## function sets meaningful defaults for selected arguments of function
  ## nleqslv{nleqslv} 
  
  ## 2013-07-12 A. Papritz
  ## 2014-07-29 AP
  ## 2015-06-30 AP function and arguments renamed
  
  list( 
    method = match.arg( method ),
    global = match.arg( global ),
    xscalm = match.arg( xscalm ),
    control = control
  )
}

## ======================================================================
control.optim <-
  function( 
    method = c( "BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN", "Brent" ),
    lower = -Inf, upper = Inf,
    control = list(reltol = 1.e-5),
    ...
  )
{
  
  ## function sets meaningful defaults for selected arguments of function optim
  ## 2012-12-14 A. Papritz
  ## 2014-07-29 AP
  ## 2015-06-30 AP function and arguments renamed
  
  list( 
    method = match.arg( method ),
    lower = lower, upper = upper,
    control = control
  )
}

## ======================================================================
control.nlminb <-
  function( 
    control = list( rel.tol = 1.e-5 ),
    lower = -Inf, upper = Inf,
    ...
  )
{
  
  ## function sets meaningful defaults for selected arguments of function nlminb
  ## 2015-07-17 A. Papritz
  
  list( 
    control = control,
    lower = lower, upper = upper
  )
}

## ======================================================================
control.pmm <-
  function( 
    ncores = 1, max.ncores = detectCores(), 
    f = 2, sfstop = FALSE, allow.recursive = TRUE, 
    ... 
  )
{
  
  ## function sets meaningful defaults for parallelized computations
  
  ## 2014-07-29 AP
  ## 2015-06-30 AP function and arguments renamed
  ## 2015-07-29 AP changes for elimination of parallelized computation of gradient or estimating equations
  
  ncores <- min( ncores, max.ncores )
  
  list( 
    ncores = ncores,  max.ncores = max.ncores, 
    f = f, sfstop = sfstop, 
    allow.recursive = allow.recursive
  )

}

## ======================================================================

compress <- 
  function( m )
{
  
  ## function stores a list of or a single lower, upper triangular or
  ## symmetric matrix compactly

  ## 2013-06-12 AP substituting [["x"]] for $x in all lists

  aux <- function( x ){
    struc <- attr( x, "struc" )
    if( !is.null( struc ) ){
      switch( 
        struc,
        sym = {
          aux <- list( diag = diag( x ), tri = x[lower.tri(x)] )
          attr( aux, "struc" ) <- "sym"
          aux
        },
        lt = {
          aux <- list( diag = diag( x ), tri = x[lower.tri(x)] )
          attr( aux, "struc" ) <- "lt"
          aux
        }
        ,
        ut = {
          aux <- list( diag = diag( x ), tri = x[upper.tri(x)] )
          attr( aux, "struc" ) <- "ut"
          aux
        }
      )
    } else {
      x
    }
  }
  
  if( is.list( m ) ){
    lapply( m, aux )
  } else {
    aux ( m )
  }
  
  
  
}

## ======================================================================

expand <- 
  function( object )
{
  
  ## function expands a list of or a compactly stored lower, upper
  ## triangular or symmetric matrices

  ## 2013-06-12 AP substituting [["x"]] for $x in all lists

  aux <- function( x ){
    struc <- attr( x, "struc" )
    if( !is.null( struc ) ){
      switch( 
        struc,
        sym = {
          n <- length( x[["diag"]] )
          dn <- names( x[["diag"]] )
          aux <- matrix( 0., n, n )
          aux[lower.tri( aux )] <- x[["tri"]]
          aux <- aux + t( aux )
          diag( aux ) <- x[["diag"]]
          dimnames( aux ) <- list( dn, dn )
          attr( aux, "struc" ) <- "sym"
          aux
        },
        lt = {
          n <- length( x[["diag"]] )
          dn <- names( x[["diag"]] )
          aux <- matrix( 0., n, n )
          aux[lower.tri( aux )] <- x[["tri"]]
          diag( aux ) <- x[["diag"]]
          dimnames( aux ) <- list( dn, dn )
          attr( aux, "struc" ) <- "lt"
          aux
        }
        ,
        ut = {
          n <- length( x[["diag"]] )
          dn <- names( x[["diag"]] )
          aux <- matrix( 0., n, n )
          aux[upper.tri( aux )] <- x[["tri"]]
          diag( aux ) <- x[["diag"]]
          dimnames( aux ) <- list( dn, dn )
          attr( aux, "struc" ) <- "ut"
          aux
        }
      )
    } else {
      x
    }
  }
  
  ln <- names( object )
  if( is.list( object ) ){
    if( length( ln ) == 2 && all( ln == c( "diag", "tri" ) ) ){
      aux( object )
    } else {
      lapply( object, aux )
    }
  } else {
    object
  }  
  
}

## ======================================================================

param.names <- 
  function( model )
{
  
  ## function returns names of extra parameters of implemented variogram
  ## models (cf. Variogram{RandomFields})
  
  ## 2012-01-24 A. Papritz
  ## 2014-05-15 AP changes for version 3 of RandomFields
  
  switch(
    model,
    "RMaskey"         = "alpha",
    "RMbessel"        = "nu",
    "RMcauchy"        = "gamma",
    "RMcircular"      = NULL,
    "RMcubic"         = NULL,
    "RMdagum"         = c( "beta", "gamma" ),
    "RMdampedcos"     = "lambda",
    "RMdewijsian"     = "alpha",
    "RMexp"           = NULL,
    "RMfbm"           = "alpha",
    "RMgauss"         = NULL,
    "RMgencauchy"     = c( "alpha", "beta" ),
    "RMgenfbm"        = c( "alpha", "delta" ),
    "RMgengneiting"   = c( "kappa", "mu" ),
    "RMgneiting"      = NULL,
    "RMlgd"           = c( "alpha", "beta" ),
    "RMmatern"        = "nu",
    "RMpenta"         = NULL,
    "RMqexp"          = "alpha",
    "RMspheric"       = NULL,
    "RMstable"        = "alpha",
    "RMwave"          = NULL,
    "RMwhittle"       = "nu",
    stop( model, " variogram not implemented" )
  )
}

##  ##############################################################################

param.bounds <- 
function( model, d )
{
  
  ## function returns range of parameters for which variogram models are
  ## valid (cf.  Variogram{RandomFields})
  
  ## 2012-03-30 A. Papritz
  ## 2014-05-15 AP changes for version 3 of RandomFields
  
  switch(
    model,
    "RMaskey"         = list( alpha = c( 0.5 * (d + 1), Inf ) ),
    "RMbessel"        = list( nu = c( 0.5 * (d - 2), Inf ) ),
    "RMcauchy"        = list( gamma = c( 1.e-18, Inf ) ),
    "RMcircular"      = NULL,
    "RMcubic"         = NULL,
    "RMdagum"         = list( beta = c( 1.e-18, 1.), gamma = c( 1.e-18, 1.-1.e-18) ),
    "RMdampedcos"     = list( lambda = c( if( d > 2 ) sqrt(3.) else 1., Inf ) ),
    "RMdewijsian"     = list( alpha = c( 1.e-18, 2. ) ),
    "RMexp"           = NULL,
    "RMfbm"           = list( alpha = c( 1.e-18, 2.) ),
    "RMgauss"         = NULL,
    "RMgencauchy"     = list( alpha = c(1.e-18, 2.), beta = c(1.e-18, Inf) ),
    "RMgenfbm"        = list( alpha = c(1.e-18, 2.), delta = c(1.e-18, 1.-1.e-18) ),
    "RMgengneiting"   = list( kappa = c(1, 3), mu = c( d/2, Inf ) ),
    "RMgneiting"      = NULL,
    "RMlgd"           = list( 
                        alpha = c( 
                          1.e-18, 
                          if( d <= 3 ) 0.5 * (3-d) else stop("dimension > 3 not allowed for RMlgd model" ) 
                        ), 
                        beta = c(1.e-18, Inf)
                      ),
    "RMmatern"        = list( nu = c(1.e-18, Inf) ),
    "RMpenta"         = NULL,
    "RMqexp"          = list( alpha = c(0., 1.) ),
    "RMspheric"       = NULL,
    "RMstable"        = list( alpha = c(1.e-18, 2.) ),
    "RMwave"          = NULL,
    "RMwhittle"       = list( nu = c(1.e-18, Inf) ),
    stop( model, " variogram not implemented" )
  )
}

 
##  ##############################################################################

default.fit.param <- 
function(
  variance = TRUE, snugget = FALSE, nugget = TRUE, scale = TRUE, 
  alpha = FALSE, beta = FALSE, delta = FALSE, gamma = FALSE, 
  kappa = FALSE, lambda = FALSE, mu = FALSE, nu = FALSE )
{
  
  ## function sets default flags for fitting variogram parameters
  
  ## 2015-11-27 A. Papritz
  
  c( 
	variance = variance, snugget = snugget, nugget = nugget, scale = scale,
	alpha = alpha, beta = beta, delta = delta, gamma = gamma,
	kappa = kappa, lambda = lambda, mu = mu, nu = nu
  )
  
}

 
##  ##############################################################################

default.fit.aniso <- 
function( f1 = FALSE, f2 = FALSE, omega = FALSE, phi = FALSE, zeta = FALSE )
{
  
  ## function sets default flags for fitting anisotropy parameters
  
  ## 2015-11-27 A. Papritz
  
  c( f1 = f1, f2 = f2, omega = omega, phi = phi, zeta = zeta )
  
}

 
##  ##############################################################################

default.aniso <- 
function(
  f1 = 1., f2 = 1., omega = 90., phi = 90., zeta = 0. )
{
  
  ## function sets default values for anisotropy parameters
  
  ## 2015-11-26 A. Papritz
  
  c( f1 = f1, f2 = f2, omega = omega, phi = phi, zeta = zeta )
  
}

 
##  ##############################################################################

profilelogLik <- 
function( object, values, use.fitted = TRUE, verbose = 0, 
  ncores = min( detectCores(), NROW(values) ) ){
  
  ## function to compute (restricted) likelihood profile for a georob fit
  
  ## arguments:
  
  ## object             a fitted georob object
  ## values             a matrix or dataframe with the fixed variogram parameter 
  ## use.fitted         control whether fitted variogram parameters should be used as initial values
  ##                    values for which remaining parameters are estimated
  ## ncores             number of cores for parallelized computations
  
  ## 2015-03-18 A. Papritz
  ## 2015-04-08 AP changes in returned results
  
  ## auxiliary function to fit model and return maximized (pseudo) log-likelihood
  
  f.aux <- function( 
    i, values, object, data, 
    fixed.param, param, fit.param, 
    fixed.aniso, aniso, fit.aniso,
    reparam, verbose
  ){
    
    values <- values[i, ]
    if( length( fixed.param ) ) param[fixed.param] <- values[fixed.param]
    if( length( fixed.aniso ) ) aniso[fixed.aniso] <- values[fixed.aniso]
    
    ## set hessian equal to FALSE in control argument of georob call and update
    ## call
    
    cl <- object[["call"]]
    
    if( "control" %in% names(cl) ){
      
      ## georob called with control argument
    
      cl.control <- as.list( cl[["control"]] )
      cl <- cl[ -match( "control", names(cl) ) ]
      if( "hessian" %in% names(cl.control) ){
        cl.control["hessian"] <- list( hessian = FALSE )
      } else {
        cl.control <- c( cl.control, hessian = FALSE )
      }
      
    } else {
      
      ## georob called without control argument
      
      cl.control <- list( as.symbol("control.georob"), hessian = FALSE )
      
    }
    
    object[["call"]] <- as.call(  c( as.list(cl), control = as.call(cl.control) ) )
    object <- update( object )
    
    fit <- update( 
      object, data = data, param = param, fit.param = fit.param,
      aniso = aniso, fit.aniso = fit.aniso, verbose = verbose
    )
    
    c( 
      loglik = logLik( 
        fit, warn = FALSE, REML = identical( object[["control"]][["ml.method"]], "REML" ) 
      ), 
      fit[["param"]][fit.param],
      fit[["aniso"]][["aniso"]][fit.aniso],      
      coef( fit ),
      converged = fit[["converged"]]
    )
    
  }
  
  ## check whether all mandatory arguments have been provided
  
  if( missing(object) || missing(values) ) stop(
	"some mandatory arguments are missing" 
  )
  
  ## warning for robust fits
  
  if( object[["tuning.psi"]] < object[["control"]][["tuning.psi.nr"]] ){
    warning( 
      "likelihood approximated for robustly fitted model by likelihood of\n",
      "  equivalent Gaussian model with heteroscedastic nugget"
    )
  }
  
  if( !(is.matrix(values) || is.data.frame( values )) ) stop( 
    "'values' must be a dataframe or a matrix" 
  )
  
  ## get data.frame with required variables (note that the data.frame passed
  ## as data argument to georob must exist in GlobalEnv)
  
  data <- cbind(
    get_all_vars( 
      formula( object ), data = eval( getCall(object)[["data"]] )
    ),
    get_all_vars( 
      object[["locations.objects"]][["locations"]], eval( getCall(object)[["data"]] )
    )
  )
  
  if( identical( class( object[["na.action"]] ), "omit" ) ) data <- na.omit(data)
  
  ## select subset if appropriate
  
  if( !is.null( getCall(object)[["subset"]] ) ){
   data <- data[eval( getCall(object)[["subset"]] ), ]
  }
  
  ## check names of fixed variogram parameters
  
  fixed.param.aniso <- colnames( values )
  if( any( 
      !fixed.param.aniso %in% gsub( 
        "(fixed)", "", 
        c( names( object[["param"]] ), names( object[["aniso"]][["aniso"]] ) ),
        fixed = TRUE
      )
    ) 
  ) stop( "column names of 'values' do not match names of variogram parameters" )
  
  ## determine variogram parameters that should be kept fixed at specified
  ## values and update call of object
  
  fixed.param <- fixed.param.aniso[fixed.param.aniso %in% names( object[["param"]] )]
  if( use.fitted ){
    param <- object[["param"]]
  } else {
    param <- object[["initial.objects"]][["param"]]
  }
  fit.param <- object[["initial.objects"]][["fit.param"]]
  if( length( fixed.param ) ) fit.param[fixed.param] <- FALSE
  
  fixed.aniso <- fixed.param.aniso[fixed.param.aniso %in% names( object[["aniso"]][["aniso"]] )]
  if( use.fitted ){
    aniso <- object[["aniso"]][["aniso"]]
  } else {
    aniso <- object[["initial.objects"]][["aniso"]]
  }
  fit.aniso <- object[["initial.objects"]][["fit.aniso"]]
  if( length( fixed.aniso ) ) fit.aniso[fixed.aniso] <- FALSE
  
  ## update object call to avoid computation of covariance matrices
  ## and to set reparam = FALSE if variance parameters are fitted
  
  reparam <- !any( colnames( values ) %in% c( "variance", "snugget", "nugget" ) )
  
  cl <- object[["call"]]

  if( "control" %in% names( cl ) ){
    
    ## georob called with control argument
    
    cl.control <- as.list( cl[["control"]] )
    cl <- cl[ -match( "control", names(cl) ) ]
    cl.control <- cl.control[ !names( cl.control ) %in% c( 
      "force.gradient", "cov.bhat", "cov.betahat", "cov.bhat.betahat", 
      "cov.delta.bhat", "cov.delta.bhat.betahat", "cov.ehat", "cov.ehat.p.bhat",
      "aux.cov.pred.target", "reparam"
    )]
    cl.control <- c(
      cl.control,
      "force.gradient" = FALSE, 
      "cov.bhat" = FALSE, "cov.betahat" = FALSE, "cov.bhat.betahat" = FALSE, 
      "cov.delta.bhat" = FALSE, "cov.delta.bhat.betahat" = FALSE,
      "cov.ehat" = FALSE,  "cov.ehat.p.bhat" = FALSE, 
      "aux.cov.pred.target" = FALSE,
      reparam = reparam
    )
    object[["call"]] <- as.call( c( as.list(cl), control = as.call(cl.control) ) )
    
  } else {
    
    ## georob called without control argument
    
    cl.control <- list( 
      as.symbol("control.georob"),
      force.gradient = FALSE, cov.bhat = FALSE, 
      cov.betahat = FALSE, cov.bhat.betahat = FALSE, 
      cov.delta.bhat = FALSE, cov.delta.bhat.betahat = FALSE,
      cov.ehat = FALSE,  cov.ehat.p.bhat = FALSE, 
      aux.cov.pred.target = FALSE,
      reparam = reparam
    )
  
  }
  
  object[["call"]] <- as.call( c( as.list(cl), control = as.call(cl.control) ) )

  ## loop over all elements of values
  
  values <- as.matrix( values )
  
  if( ncores > 1 && .Platform[["OS.type"]] == "windows" ){
    
    ## create a SNOW cluster on windows OS
    
    cl <- makePSOCKcluster( ncores, outfile = "")
    
    ## export required items to workers
    
    junk <- clusterEvalQ( cl, require( georob, quietly = TRUE ) )
    
    result <- parLapply(
      cl, 
      1:NROW(values),
      f.aux, 
      values = values, object = object, data = data,
      fixed.param = fixed.param, param = param, fit.param = fit.param,
      fixed.aniso = fixed.aniso, aniso = aniso, fit.aniso = fit.aniso,
      reparam = reparam, verbose = verbose
    )
    
    stopCluster(cl)
    
  } else {
        
    ## fork child processes on non-windows OS
    
    result <- mclapply(
      1:NROW(values),
      f.aux, 
      values = values, object = object, data = data, 
      fixed.param = fixed.param, param = param, fit.param = fit.param,
      fixed.aniso = fixed.aniso, aniso = aniso, fit.aniso = fit.aniso,
      reparam = reparam, verbose = verbose,
      mc.cores = ncores,
      mc.allow.recursive = FALSE
    )
    
  }
  
  ## collect results
  
  as.data.frame( cbind( values, t( simplify2array( result ) ) ) )
  
}

