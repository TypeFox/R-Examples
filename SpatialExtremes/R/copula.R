fitcopula <- function(data, coord, copula = "gaussian", cov.mod = "whitmat",
                      loc.form, scale.form, shape.form, marg.cov = NULL, temp.cov = NULL,
                      temp.form.loc = NULL, temp.form.scale = NULL,
                      temp.form.shape = NULL, ..., start,
                      control = list(maxit = 10000), method = "Nelder",
                      std.err = TRUE, warn = TRUE, corr = FALSE){

  n.site <- ncol(data)
  n.obs <- nrow(data)
  dist.dim <- ncol(coord)
  std.err2 <- std.err
  
  dist <- t(as.matrix(dist(coord, diag = TRUE)))
  dist <- dist[lower.tri(dist, diag = TRUE)]
  
  if (!(cov.mod %in% c("whitmat","cauchy","powexp","bessel","caugen")))
    stop("''cov.mod'' must be one of 'whitmat', 'cauchy', 'powexp', 'bessel', 'caugen'")

  if (!(copula %in% c("gaussian", "student")))
    stop("''copula'' must be one of 'gaussian' or 'student'")
  
  if (cov.mod == "whitmat")
    cov.mod.num <- 1
  if (cov.mod == "cauchy")
    cov.mod.num <- 2
  if (cov.mod == "powexp")
    cov.mod.num <- 3
  if (cov.mod == "bessel")
    cov.mod.num <- 4
  if (cov.mod == "caugen")
    cov.mod.num <- 5

  if (copula == "gaussian")
    copula.num <- 1

  else
    copula.num <- 2

  if (missing(loc.form) && missing(scale.form) && missing(shape.form))
    fit.marge <- FALSE

  if (!missing(loc.form) && !missing(scale.form) && !missing(shape.form)){
    fit.marge <- TRUE

    if ((class(loc.form) != "formula") || (class(scale.form) != "formula") ||
        (class(shape.form) != "formula"))
      stop("''loc.form'', ''scale.form'' and ''shape.form'' must be valid R formulas")
  }

  flag <- missing(loc.form) + missing(scale.form)  + missing(shape.form)

  if (!(flag %in% c(0, 3)))
    stop("if one formula is given for the GEV parameters, then it should
be given for *ALL* GEV parameters")

  param <- c("nugget", "range", "smooth")

  if (cov.mod == "caugen")
    param <- c(param, "smooth2")

  else
    ##Fix it to 0 as it won't be used anyway
    smooth2 <- 0
  
  if (copula == "student")
    param <- c("DoF", param)

  else
    ##Fix it to 0 as it won't be used anyway
    DoF <- 0

  if (n.site != nrow(coord))
    stop("'data' and 'coord' doesn't match")

  use.temp.cov <- c(!is.null(temp.form.loc), !is.null(temp.form.scale),
                    !is.null(temp.form.shape))

  if (any(use.temp.cov) && (n.obs != nrow(temp.cov)))
    stop("'data' and 'temp.cov' doesn't match")

  if (any(use.temp.cov) && is.null(temp.cov))
    stop("'temp.cov' must be supplied if at least one temporal formula is given")

  if (fit.marge){
    ##With our notation, formula must be of the form y ~ xxxx
    loc.form <- update(loc.form, y ~ .)
    scale.form <- update(scale.form, y ~ .)
    shape.form <- update(shape.form, y ~ .)
    
    if (use.temp.cov[1])
      temp.form.loc <- update(temp.form.loc, y ~. + 0)
    
    if (use.temp.cov[2])
      temp.form.scale <- update(temp.form.scale, y ~. + 0)
    
    if (use.temp.cov[3])
      temp.form.shape <- update(temp.form.shape, y ~. + 0)
    
    if (is.null(marg.cov))
      covariables <- data.frame(coord)
    
    else
      covariables <- data.frame(coord, marg.cov)
    
    loc.model <- modeldef(covariables, loc.form)
    scale.model <- modeldef(covariables, scale.form)
    shape.model <- modeldef(covariables, shape.form)
    
    loc.dsgn.mat <- loc.model$dsgn.mat
    scale.dsgn.mat <- scale.model$dsgn.mat
    shape.dsgn.mat <- shape.model$dsgn.mat
    
    loc.pen.mat <- loc.model$pen.mat
    scale.pen.mat <- scale.model$pen.mat
    shape.pen.mat <- shape.model$pen.mat
    
    loc.penalty <- loc.model$penalty.tot
    scale.penalty <- scale.model$penalty.tot
    shape.penalty <- shape.model$penalty.tot
    
    ##The total number of parameters to be estimated for each GEV
    ##parameter
    n.loccoeff <- ncol(loc.dsgn.mat)
    n.scalecoeff <- ncol(scale.dsgn.mat)
    n.shapecoeff <- ncol(shape.dsgn.mat)
    
    ##The number of ``purely parametric'' parameters to estimate i.e. we
    ##do not consider the weigths given to each basis function
    n.pparloc <- loc.model$n.ppar
    n.pparscale <- scale.model$n.ppar
    n.pparshape <- shape.model$n.ppar
    
    loc.names <- paste("locCoeff", 1:n.loccoeff, sep="")
    scale.names <- paste("scaleCoeff", 1:n.scalecoeff, sep="")
    shape.names <- paste("shapeCoeff", 1:n.shapecoeff, sep="")
  }

  else {
    loc.dsgn.mat <- scale.dsgn.mat <- shape.dsgn.mat <- loc.pen.mat <-
      scale.pen.mat <- shape.pen.mat <- loc.names <- scale.names <-
        shape.names <- loc.form <- scale.form <- shape.form <- NULL
    n.loccoeff <- n.scalecoeff <- n.shapecoeff <- n.pparloc <-
      n.pparscale <- n.pparshape <- loc.penalty <- scale.penalty <-
        shape.penalty <- 0
  }
  
  ##Do the same for the temporal regression coefficients
  if (use.temp.cov[1]){
    temp.model.loc <- modeldef(temp.cov, temp.form.loc)
    temp.dsgn.mat.loc <- temp.model.loc$dsgn.mat
    temp.pen.mat.loc <- temp.model.loc$pen.mat
    temp.penalty.loc <- temp.model.loc$penalty.tot
    n.tempcoeff.loc <- ncol(temp.dsgn.mat.loc)
    n.ppartemp.loc <- temp.model.loc$n.ppar
    temp.names.loc <- paste("tempCoeffLoc", 1:n.tempcoeff.loc, sep="")
  }
  
  else {
    temp.model.loc <- temp.dsgn.mat.loc <- temp.pen.mat.loc <- temp.names.loc <- NULL
    n.tempcoeff.loc <- n.ppartemp.loc <- temp.penalty.loc <- 0
  }
  
  if (use.temp.cov[2]){
    temp.model.scale <- modeldef(temp.cov, temp.form.scale)
    temp.dsgn.mat.scale <- temp.model.scale$dsgn.mat
    temp.pen.mat.scale <- temp.model.scale$pen.mat
    temp.penalty.scale <- temp.model.scale$penalty.tot
    n.tempcoeff.scale <- ncol(temp.dsgn.mat.scale)
    n.ppartemp.scale <- temp.model.scale$n.ppar
    temp.names.scale <- paste("tempCoeffScale", 1:n.tempcoeff.scale, sep="")
  }
  
  else {
    temp.model.scale <- temp.dsgn.mat.scale <- temp.pen.mat.scale <- temp.names.scale <- NULL
    n.tempcoeff.scale <- n.ppartemp.scale <- temp.penalty.scale <- 0
  }
  
  if (use.temp.cov[3]){
    temp.model.shape <- modeldef(temp.cov, temp.form.shape)
    temp.dsgn.mat.shape <- temp.model.shape$dsgn.mat
    temp.pen.mat.shape <- temp.model.shape$pen.mat
    temp.penalty.shape <- temp.model.shape$penalty.tot
    n.tempcoeff.shape <- ncol(temp.dsgn.mat.shape)
    n.ppartemp.shape <- temp.model.shape$n.ppar
    temp.names.shape <- paste("tempCoeffShape", 1:n.tempcoeff.shape, sep="")
  }
  
  else {
    temp.model.shape <- temp.dsgn.mat.shape <- temp.pen.mat.shape <- temp.names.shape <- NULL
    n.tempcoeff.shape <- n.ppartemp.shape <- temp.penalty.shape <- 0
  }
  
  param <- c(param, loc.names, scale.names, shape.names, temp.names.loc, temp.names.scale,
             temp.names.shape)

  nllik <- function(x) x
  
  body(nllik) <- parse(text = paste("-.C('copula', as.integer(copula.num), as.integer(cov.mod.num),
as.double(dist), as.double(data), as.integer(n.site), as.integer(n.obs), as.integer(dist.dim), as.integer(fit.marge),
as.double(loc.dsgn.mat), as.double(loc.pen.mat), as.integer(n.loccoeff), as.integer(n.pparloc),
as.double(loc.penalty), as.double(scale.dsgn.mat), as.double(scale.pen.mat),
as.integer(n.scalecoeff), as.integer(n.pparscale), as.double(scale.penalty),
as.double(shape.dsgn.mat), as.double(shape.pen.mat), as.integer(n.shapecoeff),
as.integer(n.pparshape), as.double(shape.penalty), as.integer(use.temp.cov),
as.double(temp.dsgn.mat.loc), as.double(temp.pen.mat.loc), as.integer(n.tempcoeff.loc),
as.integer(n.ppartemp.loc), as.double(temp.penalty.loc), as.double(temp.dsgn.mat.scale),
as.double(temp.pen.mat.scale), as.integer(n.tempcoeff.scale), as.integer(n.ppartemp.scale),
as.double(temp.penalty.scale), as.double(temp.dsgn.mat.shape), as.double(temp.pen.mat.shape),
as.integer(n.tempcoeff.shape), as.integer(n.ppartemp.shape), as.double(temp.penalty.shape),",
                         paste("as.double(c(", paste(loc.names, collapse = ","), ")), "),
                         paste("as.double(c(", paste(scale.names, collapse = ","), ")), "),
                         paste("as.double(c(", paste(shape.names, collapse = ","), ")), "),
                         paste("as.double(c(", paste(temp.names.loc, collapse = ","), ")), "),
                         paste("as.double(c(", paste(temp.names.scale, collapse = ","), ")), "),
                         paste("as.double(c(", paste(temp.names.shape, collapse = ","), ")), "),
                         "as.double(DoF), as.double(nugget), as.double(range), as.double(smooth),
as.double(smooth2), dns = double(1), PACKAGE = 'SpatialExtremes')$dns"))
  
  
  ##Define the formal arguments of the function
  form.nllik <- NULL
  for (i in 1:length(param))
    form.nllik <- c(form.nllik, alist(a=))
  
  names(form.nllik) <- param
  formals(nllik) <- form.nllik
  
  if (missing(start)){

    start <- list(nugget = 0, range = 0.25 * max(dist), smooth = 1)
    
    if (copula == "student")
      start <- c(list(DoF = 1), start)
    
    if (fit.marge){
      loc <- scale <- shape <- rep(0, n.site)
    
      for (i in 1:n.site){
        gev.param <- gevmle(data[,i])
        loc[i] <- gev.param["loc"]
        scale[i] <- gev.param["scale"]
        shape[i] <- gev.param["shape"]
      }
      
      locCoeff <- loc.model$init.fun(loc)
      scaleCoeff <- scale.model$init.fun(scale)
      shapeCoeff <- shape.model$init.fun(shape)
      
      locCoeff[is.na(locCoeff)] <- 0
      scaleCoeff[is.na(scaleCoeff)] <- 0
      shapeCoeff[is.na(shapeCoeff)] <- 0
      
      ##To be sure that the scale parameter is always positive at starting
      ##values
      scales.hat <- scale.model$dsgn.mat %*% scaleCoeff
      
      if (any(scales.hat <= 0))
        scaleCoeff[1] <- scaleCoeff[1] - 1.001 * min(scales.hat)
      
      names(locCoeff) <- loc.names
      names(scaleCoeff) <- scale.names
      names(shapeCoeff) <- shape.names
      
      if (use.temp.cov[1]){
        tempCoeff.loc <- rep(0, n.tempcoeff.loc)
        names(tempCoeff.loc) <- temp.names.loc
      }
      
      else
        tempCoeff.loc <- NULL
      
      if (use.temp.cov[2]){
        tempCoeff.scale <- rep(0, n.tempcoeff.scale)
        names(tempCoeff.scale) <- temp.names.scale
      }
      
      else
        tempCoeff.scale <- NULL
      
      if (use.temp.cov[3]){
        tempCoeff.shape <- rep(0, n.tempcoeff.shape)
        names(tempCoeff.shape) <- temp.names.shape
      }
      
      else
        tempCoeff.shape <- NULL
      
      start <- c(start, as.list(c(locCoeff, scaleCoeff, shapeCoeff, tempCoeff.loc, tempCoeff.scale, tempCoeff.shape)))
    }
    
    start <- start[!(param %in% names(list(...)))]
  }

  if (!length(start)) 
    stop("there are no parameters left to maximize over")
  
  nm <- names(start)
  l <- length(nm)
  f <- formals(nllik)
  names(f) <- param
  m <- match(nm, param)
  
  if(any(is.na(m))) 
    stop("'start' specifies unknown arguments")

  formals(nllik) <- c(f[m], f[-m])
  nllh <- function(p, ...) nllik(p, ...)

  if(l > 1)
    body(nllh) <- parse(text = paste("nllik(", paste("p[",1:l,
                          "]", collapse = ", "), ", ...)"))
  
  fixed.param <- list(...)[names(list(...)) %in% param]
  
  if(any(!(param %in% c(nm,names(fixed.param)))))
    stop("unspecified parameters")
  
  start.arg <- c(list(p = unlist(start)), fixed.param)

  init.lik <- do.call("nllh", start.arg)
  if (warn && (init.lik >= 1.0e6)) 
    warning("negative log-likelihood is infinite at starting values")

  if (method == "nlminb"){
    start <- as.numeric(start)
    opt <- nlminb(start, nllh, ..., control = control, hessian = std.err)
    opt$counts <- opt$evaluations
    opt$value <- opt$objective
    names(opt$par) <- nm
  }
  
  else if (method == "nlm"){
    start <- as.numeric(start)
    opt <- nlm(nllh, start, ..., hessian = std.err)
    opt$counts <- opt$iterations
    names(opt$counts) <- "function"
    opt$value <- opt$minimum
    opt$par <- opt$estimate
    names(opt$par) <- nm

    if (opt$code <= 2){
      opt$convergence <- 0
      opt$message <- NULL
    }
    
    if (opt$code > 2){
      opt$convergence <- 1
      opt$message <- paste("nlm error code", opt$code)
    }      
  }

  else 
    opt <- optim(start, nllh, ..., method = method, control = control,
                 hessian = std.err)
  
  if ((opt$convergence != 0) || (opt$value >= 1.0e6)){
    if (warn)
      warning("optimization may not have succeeded")
  }
  
  else
    opt$convergence <- "successful"

  param.names <- param
  param <- c(opt$par, unlist(fixed.param))
  param <- param[param.names]

  if (std.err2){
    ihessian <- try(solve(opt$hessian), silent = TRUE)
    
    if(!is.matrix(ihessian)){
      if (warn)
        warning("observed information matrix is singular; std. err. won't be computed")
      
      std.err2 <- FALSE
    }
    
    else{
      std.err <- diag(ihessian)
      std.idx <- which(std.err <= 0)
      if(length(std.idx) > 0){
        if (warn)
          warning("Some (observed) standard errors are negative;\n passing them to NA")
        
        std.err[std.idx] <- NA
      }
            
      std.err <- sqrt(std.err)
      
      if(corr) {
        .mat <- diag(1/std.err, nrow = length(std.err))
        corr.mat <- structure(.mat %*% ihessian %*% .mat, dimnames = list(nm,nm))
        diag(corr.mat) <- rep(1, length(std.err))
      }
      
      else
        corr.mat <- NULL
      
      colnames(ihessian) <- rownames(ihessian) <- names(std.err) <- nm
    }
  }

  if (!std.err2)
    std.err <- std.err.type <- corr.mat <- var.cov <- ihessian <-
      var.score <- NULL

  if (cov.mod == "caugen")
    cov.fun <- covariance(nugget = param["nugget"], sill = 1 - param["nugget"], range = param["range"],
                          smooth = param["smooth"], smooth2 = param["smooth2"],
                          cov.mod = cov.mod, plot = FALSE)

  else
    cov.fun <- covariance(nugget = param["nugget"], sill = 1 - param["nugget"], range = param["range"],
                          smooth = param["smooth"], cov.mod = cov.mod, plot = FALSE)

  if (copula == "gaussian")
    ext.coeff <- function(h)
      rep(2, length(h))

  else
    ext.coeff <- function(h)
      2 * (1 - pt(-sqrt((param["DoF"] + 1) * (1 - cov.fun(h)) /
                        (1 + cov.fun(h))), param["DoF"] + 1))
  
  ans <- list(fitted.values = opt$par, param = param, std.err = std.err,
              counts = opt$counts, message = opt$message, coord = coord,
              logLik = -opt$value, loc.form = loc.form, scale.form = scale.form,
              shape.form = shape.form, convergence = opt$convergence,
              nllh = nllh, deviance = 2 * opt$value, var.cov = ihessian,
              data = data, fixed = unlist(fixed.param), hessian = opt$hessian,
              marg.cov = marg.cov, use.temp.cov = use.temp.cov, corr = corr.mat,
              copula = copula, ext.coeff = ext.coeff, cov.fun = cov.fun,
              cov.mod = cov.mod, fit.marge = fit.marge, iso = TRUE)
  
  class(ans) <- "copula"
  return(ans)
}

rcopula <- function(n, coord, copula = "gaussian", cov.mod = "whitmat", grid = FALSE,
                    control = list(), nugget = 0, range = 1, smooth = 1, DoF = 1){

  ## This function simulates realizations from a copula with unit
  ## Frechet margins

  if (!(cov.mod %in% c("whitmat","cauchy","powexp","bessel")))
    stop("'cov.mod' must be one of 'gauss', 'whitmat', 'cauchy', 'powexp' or 'bessel'")

  if ((nugget > 1) || (nugget < 0))
    stop("'nugget' should lie in [0, 1]")

  sill <- 1 - nugget

  gp <- rgp(n, coord, cov.mod, grid = grid, control = list(), nugget = nugget,
            sill = sill, range = range, smooth = smooth)
    
  if (copula == "gaussian")
    ans <- qgev(pnorm(gp), 1, 1, 1)

  if (copula == "student"){
    if (DoF <= 0)
      stop("'DoF' has to be positive.")

    scalings <- sqrt(DoF / rchisq(n, DoF))

    if (!grid)
      ans <- gp * scalings

    else {
      ans <- array(NA, dim = dim(gp))
      for (i in 1:n)
        ans[,,i] <- gp[,,i] * scalings[i]
    }

    ans <- qgev(pt(ans, DoF), 1, 1, 1)
  }

  return(ans)
}
   

                     

  
