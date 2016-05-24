##
## Auxiliary functions for the geoR library
## ----------------------------------------
##
## These functions are typically called by other main functions
## to perform internal calculations
##

".geoR.check.locations" <- ".check.locations" <-
  function(locations)
{
  if(all(locations == "no")) return("no")
  if(any(is.na(locations))) stop("NA's cannot be provided for the argument locations")
  if(inherits(locations, "SpatialPoints"))
    locations <- coordinates(locations)
  if(!is.list(locations) && is.vector(locations)) {
    ##
    ## Checking the spatial dimension for prediction
    ##  1 (data/prediction on a transect) or 2 (data/prediction on an area)
    ##
    if(length(locations) == 2) {
      locations <- t(as.matrix(locations))
      warning("assuming that there is only 1 prediction point")
    }
    else{
      warning("locations provided as a vector, assuming one spatial dimension")
      locations <- as.matrix(cbind(locations, 0))
    }
  }
  else{
    if(is.matrix(locations) | is.data.frame(locations)){
      if(ncol(locations) < 2)
        stop("locations must be a 2 column matrix, data-frame or a list with 2 components")
      if(ncol(locations) > 2){
        warning("locations provided with a matrix or data-frame with more than 2 columns. Only the first two columns used as coordinates")
        locations <- locations[,1:2]
      }
    }
    else{
      if(is.list(locations)){
        if(length(locations) < 2)
          stop("locations must be a 2 column matrix, data-frame or a list with 2 components")
        if(length(locations) > 2)
          warning("locations provided as a list with more than 2 components. Only the 2 first will be used as coordinates")
        locations <- matrix(unlist(locations[1:2]), ncol=2)
      }
    }
  }
  return(as.matrix(locations))
}

".solve.geoR" <-
  function (a, b = NULL, ...) 
{
  a <- eval(a)
  b <- eval(b)
#  if(exists("trySilent")){
    if (is.null(b)) res <- try(solve(a, ...), silent=TRUE)
    else res <- try(solve(a, b, ...), silent=TRUE)
#  }
#  else{
#    error.now <- options()$show.error.messages
#    if (is.null(error.now) | error.now) 
#      on.exit(options(show.error.messages = TRUE))
#    options(show.error.messages = FALSE)
#    if (is.null(b)) res <- try(solve(a, ...))
#    else res <- try(solve(a, b, ...))
#  }
  if (inherits(res, "try-error")) {
    test <- all.equal.numeric(a, t(a), 100 * .Machine$double.eps)
    if(!(is.logical(test) && test)){
      ##      options(show.error.messages = TRUE)
      stop("matrix `a' is not symmetric")
    }
    t.ei <- eigen(a, symmetric = TRUE)
#    if(exists("trySilent")){
      if (is.null(b))
        res <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)), silent=TRUE)
      else
        res <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)) %*% b, silent=TRUE)
#    }
#    else{
#      if (is.null(b)) res <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)))
#      else res <- try(crossprod(t(t.ei$vec)/sqrt(t.ei$val)) %*% b)
#    }
    if (any(is.na(res)) | any(is.nan(res)) | any(is.infinite(res))) 
      oldClass(res) <- "try-error"
  }
  if (inherits(res, "try-error")) 
    stop("Singular matrix. Covariates may have different orders of magnitude.")
  return(res)
}

"coords.aniso" <- 
  function(coords, aniso.pars, reverse=FALSE)
{
  coords <- as.matrix(coords)
  n <- nrow(coords)
  if(length(aniso.pars) != 2)
    stop("argument aniso.pars must be a vector with 2 elementsm the anisotropy angle and anisotropy ratio, respectively")
  psiA <- aniso.pars[1]
  psiR <- aniso.pars[2]
  if(psiR < 1){
    psiR <- round(psiR, digits=8)
    if(psiR < 1)
      stop("anisotropy ratio must be greater than 1")
  }
  rm <- matrix(c(cos(psiA), -sin(psiA), sin(psiA), cos(psiA)), ncol = 2)
  tm <- diag(c(1, 1/psiR))
  if(reverse)
    coords.mod <- coords %*% solve(rm %*% tm)
  else
    coords.mod <- coords %*% rm %*% tm
  return(coords.mod)
}

#"dist0.krige" <-
#function (x0, coords) 
#{
#  if (length(x0) != 2) 
#    stop(paste("length of x0 is", length(x0), "(it must be 2)"))
#  coords[, 1] <- coords[, 1] - x0[1]
#  coords[, 2] <- coords[, 2] - x0[2]
#  return(sqrt(coords[, 1]^2 + coords[, 2]^2))
#}


"set.coords.lims" <-
  function(coords, borders = coords, xlim, ylim, ...)
{
  if(!is.null(borders)){
    if(ncol(borders) != 2)
      stop("argument borders must be an object with 2 columns with the XY coordinates of the borders of the area")
    coords <- rbind(coords, as.matrix(as.data.frame(borders)))
  }
  coords.lims <- apply(coords, 2, range, na.rm=TRUE)
  if(!missing(xlim) && mode(xlim) == "numeric")
    coords.lims[,1] <- xlim[order(xlim)]
  if(!missing(ylim) && mode(ylim) == "numeric")
    coords.lims[,2] <- ylim[order(ylim)]
  coords.diff <- diff(coords.lims)
  if (coords.diff[1] != coords.diff[2]) {
    coords.diff.diff <- abs(diff(as.vector(coords.diff)))
    ind.min <- which(coords.diff == min(coords.diff))
    coords.lims[, ind.min] <- coords.lims[, ind.min] + c(-coords.diff.diff, coords.diff.diff)/2
  }
  return(coords.lims)
}


"dinvchisq" <-
  function(x, df, scale=1/df, log = FALSE)
{
  if(df <= 0)
    stop("df must be greater than zero")
  if(scale <= 0)
    stop("scale must be greater than zero")
  nu <- df/2
  if(log)
    return(ifelse(x > 0, nu*log(nu) - log(gamma(nu)) + nu*log(scale) -
                  (nu+1)*log(x) - (nu*scale/x), NA))
  else
    return(ifelse(x > 0,
                  (((nu)^(nu))/gamma(nu)) * (scale^nu) *
                  (x^(-(nu+1))) * exp(-nu*scale/x), NA))
}


"rinvchisq" <- 
  function (n, df, scale = 1/df)
{
  if((length(scale)!= 1) & (length(scale) != n))
    stop("scale should be a scalar or a vector of the same length as x")
  if(df <= 0)
    stop("df must be greater than zero")
  if(any(scale <= 0))
    stop("scale must be greater than zero")
  return((df*scale)/rchisq(n, df=df)) 
}

".check.coords" <-
  function(x)
{
  xname <- deparse(substitute(x))
  if(inherits(x, "SpatialPoints")) x <- coordinates(x)
  if(is.matrix(x)){
    locnames <- colnames(x)[1:2]
    x <- x[,1:2]
    attr(x, "type") <- "matrix"
  }
  else{
    if(is.data.frame(x)){
      locnames <- names(x)[1:2]
      x <- as.matrix(x[,1:2])
      attr(x, "type") <- "dataframe"
    }
    else
      if(is.list(x)){ 
        locnames <- names(x)[1:2]
        x <- matrix(unlist(x[1:2]), ncol=2)
        attr(x, "type") <- "list"
      } 
  }
  if(is.vector(x))
    stop(paste(xname, "must be a matrix, data-frame or list with 2D coordinates"))
  colnames(x) <- locnames
  if(ncol(x) != 2)
    stop(paste(xname, "must be a matrix, data-frame or a list with 2D coordinates"))
  return(x)
}

".check.borders" <-
  function(x)
{
  xname <- deparse(substitute(x))
  if(!inherits(x, "SpatialPolygons")) {
    if (is.matrix(x)) { 
      x <- x[,1:2]
    } else if(is.data.frame(x)) { 
      x <- as.matrix(x[,1:2])
    } else if(is.list(x)) {
      x <- matrix(unlist(x[1:2]), ncol=2) 
    } else stop(paste(xname, "must be a SpatialPolygons object or a matrix"))
    if(nrow(x) < 3) stop("borders must have at least 3 points")
    if(!identical(x[1,], x[nrow(x),])) x <- rbind(x, x[1,])
    x <- SpatialPolygons(list(Polygons(list(Polygon(coords=x)), ID="borders")))
  }
  return(x)
}

  
"locations.inside" <-
  function(locations, borders, as.is = TRUE, ...)
{
  locations <- .check.coords(locations)
  borders <- .check.borders(borders)
  res <- .geoR_pip(pts = locations, poly = borders, ...)
  if(as.is){
    if(attr(locations, "type") == "dataframe")
      res <- as.data.frame(res)
    if(attr(locations, "type") == "list")
      res <- as.list(as.data.frame(res))
  }
  return(res)
}

"polygrid" <- 
  function(xgrid, ygrid, borders, vec.inout = FALSE, ...)
{
  ## checking input
  if(!is.list(xgrid) && is.vector(drop(xgrid))){
    if(missing(ygrid))
      stop("xgrid must have x and y coordinates or a vector must be provided for ygrid")
    if(!is.vector(ygrid)) stop("ygrid must be a vector")
    xygrid <- expand.grid(x = xgrid, y = ygrid)
  }
  if(is.matrix(xgrid) || is.data.frame(xgrid)){
    if(ncol(xgrid) != 2)
      stop("xgrid must be a vector or a 2 column matrix or data-frame")
    xygrid <- xgrid
    if(!missing(ygrid)) warning("xgrid has 2 column, ygrid was ignored")
  }
  else
    if(is.list(xgrid)){
      if(length(xgrid) != 2)
        stop("if xgrid is a list it must have 2 elements")
      xygrid <- expand.grid(x = xgrid[[1]], y = xgrid[[2]])
      if(!missing(ygrid)) warning("xgrid is a list, ygrid was ignored")
    }
  borders <- .check.borders(borders)
#  if(nrow(borders) < 3) stop("borders must have at least 3 points")
#  if(!identical(borders[1,], borders[nrow(borders),]) 
#    borders <- rbind(borders, borders[1,])
  ind <- as.vector(.geoR_inout(pts=xygrid, poly=borders, ...))
  xypoly <- xygrid[ind == TRUE,  ]
  if(vec.inout == FALSE)
    return(xypoly)
  else return(list(xypoly = xypoly, vec.inout = ind))
}

"trend.spatial" <-
  function (trend, geodata, add.to.trend) 
{  
  if(!missing(geodata)){
    # if(any(class(geodata) %in% ls(pattern=glob2rx("Spatial*DataFrame"), pos="package:sp")))
    if(any(class(geodata) %in% c("SpatialGridDataFrame","SpatialLinesDataFrame","SpatialPolygonsDataFrame","SpatialPointsDataFrame","SpatialPixelsDataFrame")))
        geodata <- geodata@data
    attach(geodata, pos=2, warn.conflicts=FALSE)
    if(!is.null(geodata$covariate)){
      attach(geodata$covariate, pos=3, warn.conflicts=FALSE)
      on.exit(detach("geodata$covariate"), add=TRUE)
    }
    on.exit(detach("geodata"), add=TRUE)
  }
  if(inherits(trend, "formula")) {
                                        #    require(methods)
                                        #    if(exists("trySilent")){
    trend.mat <- try(model.matrix(trend), silent=TRUE)
                                        #    }
                                        #    else{
                                        #      error.now <- options()$show.error.messages
                                        #      if (is.null(error.now) | error.now) 
    ##      options(show.error.messages = FALSE)
                                        #     trend.mat <- try(model.matrix(trend))
                                        #   }    
    if (inherits(trend.mat, "try-error")) 
      stop("\ntrend elements not found")
  }
  else {
    if(mode(trend) == "numeric")
      trend.mat <- unclass(trend)
    else if (trend == "cte"){
      if(missing(geodata))
        stop("argument geodata must be provided with trend=\"cte\"")
      trend.mat <- as.matrix(rep(1, nrow(geodata$coords)))
    }
    else if (trend == "1st"){
      if(missing(geodata))
        stop("argument geodata must be provided with trend=\"1st\"")
      trend.mat <- cbind(1, geodata$coords)
    }
    else if (trend == "2nd"){ 
      if(missing(geodata))
        stop("argument geodata must be provided with trend=\"2nd\"")
      trend.mat <- cbind(1, geodata$coords, geodata$coords[,1]^2,
                         geodata$coords[,2]^2,
                         geodata$coords[,1] * geodata$coords[,2])
    }
    else stop("external trend must be provided for data locations to be estimated using the arguments trend.d and trend.l. Allowed values are the strings \"cte\", \"1st\", \"2nd\" or  a model formula")
  }
  trend.mat <- as.matrix(trend.mat)
  if(!missing(add.to.trend)){
    if(missing(geodata))
      trend.mat <- cbind(trend.mat, trend.spatial(add.to.trend)[,-1])
    else
      trend.mat <- cbind(trend.mat, trend.spatial(add.to.trend, geodata = geodata)[,-1])
  }
  dimnames(trend.mat) <- list(NULL, NULL)
  oldClass(trend.mat) <- "trend.spatial"
  return(trend.mat)
}

".nlmP" <- function(objfunc, params, lower = rep( -Inf, length(params)),
                   upper = rep(+Inf, length(params)), ... )
{
  ## minimizer, using nlm with transformation of variables
  ## to allow for limits for the parameters   
  ##
  ## objfunc is a function to be optimised
  ## params is a starting value for the parameters
  ##
  ## NOTE: this function was used before optim() becomes available for R
  ##       It has limited usage now.
  ##
  ## Adapted from a function from Patrick E. Brown, Lancaster University 
  ##
  Nparams <- length(params)
  if(length(lower) != Nparams)
    stop(" length of lower boundaries differs from length of params")
  if(length(upper) != Nparams)
    stop(" upper boundry different length than params")
  checklimits <- upper - lower
  if(any(checklimits <= 0))
    stop(" bad boundries")
  if(any(params < lower))
    stop(" starting params too low")
  if(any(params > upper))
    stop(" starting params too high")
  
  bothlimQQ <- (lower != (-Inf)) & (upper != +Inf)
  loweronlyQQ <- (lower != (-Inf)) & (upper == +Inf)
  upperonlyQQ <- (lower == (-Inf)) & (upper != +Inf)
  ubothQQ <- upper[bothlimQQ]
  lbothQQ <- lower[bothlimQQ]
  dbothQQ <- ubothQQ - lbothQQ
  loneQQ <- lower[loweronlyQQ]
  uoneQQ <- upper[upperonlyQQ]
  
  .bounds.list <- list(bothlimQQ = bothlimQQ, 
                       loweronlyQQ = loweronlyQQ,
                       upperonlyQQ = upperonlyQQ,
                       ubothQQ = ubothQQ,
                       lbothQQ = lbothQQ,
                       dbothQQ = dbothQQ,
                       loneQQ = loneQQ,
                       uoneQQ = uoneQQ)
  
  assign(".objfuncQQ", objfunc, pos=1)
  assign(".bounds.list", .bounds.list, pos=1)
  
  ## reduce the parameter space by a scale to keep parameters
  ## away from the boundries
  
  normaltomad <- function(normalparamsX)
    {
      madparamsX <- normalparamsX
      if(any(.bounds.list$bothlimQQ)) {
        noughtone <- (normalparamsX[.bounds.list$bothlimQQ] -
                      .bounds.list$lbothQQ)/.bounds.list$dbothQQ
        madparamsX[.bounds.list$bothlimQQ] <- log(noughtone/(1 - noughtone))
      }
      
      if(any(.bounds.list$loweronlyQQ))
        madparamsX[.bounds.list$loweronlyQQ] <-
          log(normalparamsX[.bounds.list$loweronlyQQ] - .bounds.list$loneQQ)
      
      if(any(.bounds.list$upperonlyQQ))
        madparamsX[.bounds.list$upperonlyQQ] <-
          log(.bounds.list$uoneQQ - normalparamsX[.bounds.list$upperonlyQQ])
      
      return(madparamsX)
    }
#  madtonormalQQ <<- function(madparamsX)
  "madtonormalQQ" <- function(madparamsX)
  {
      normalparamsX <- madparamsX
      
      if(any(.bounds.list$bothlimQQ)) {
###        madparamsX[((.bounds.list$bothlimQQ) & (madparamsX > 300))] <- 300
        emad <- exp(madparamsX[.bounds.list$bothlimQQ])
        normalparamsX[.bounds.list$bothlimQQ] <-
          .bounds.list$dbothQQ * (emad/(1 + emad)) + .bounds.list$lbothQQ
      }
      
      if(any(.bounds.list$loweronlyQQ)){
        normalparamsX[.bounds.list$loweronlyQQ] <-
          exp(madparamsX[.bounds.list$loweronlyQQ]) + .bounds.list$loneQQ
      }
      
      if(any(.bounds.list$upperonlyQQ))
        normalparamsX[.bounds.list$upperonlyQQ] <-
          - exp(madparamsX[.bounds.list$upperonlyQQ]) + .bounds.list$uoneQQ
      
      if(exists(get(".ind.prof.phi", pos=1)))
        if(is.nan(normalparamsX[get(".ind.prof.phi", pos=1)]))
          normalparamsX[get(".ind.prof.phi", pos=1)] <- 0
      
      return(normalparamsX)
    }
  
  newobjfunc <- function(madparams) {
    normalparams <-  get("madtonormalQQ", pos=1)(madparams)
    
    get(".objfuncQQ", pos=1)(normalparams)
    
  }
  assign("madtonormalQQ", madtonormalQQ, pos=1)

  startmadparams <- normaltomad(params)
  result <- nlm(newobjfunc, startmadparams, ...)
  result$madestimate <- result$estimate
  result$estimate <- get("madtonormalQQ", pos=1)(result$madestimate)
  remove(".bounds.list", pos=1, inherits=TRUE)
  remove(".objfuncQQ", pos=1, inherits=TRUE)
  remove("madtonormalQQ", pos=1, inherits=TRUE)
  
###  return(result, madtonormalQQ(normaltomad(params)),params)
  return(result)
}


"pars.limits" <-
  function(phi = c(lower=0, upper=+Inf),
           sigmasq = c(lower=0, upper=+Inf),
           nugget.rel = c(lower=0, upper=+Inf),
           kappa = c(lower=0, upper=+Inf),
           kappa2 = c(lower=0, upper=+Inf),
           lambda = c(lower=-3, upper=3),
           psiR = c(lower=1, upper=+Inf),
           psiA = c(lower=0, upper=2*pi),
           tausq.rel = nugget.rel
           )
{
  if(length(phi) != 2)
    stop("phi must be a 2 components vector with lower and upper limits for the parameter phi") 
  if(length(sigmasq) != 2)
    stop("sigmasq must be a 2 components vector with lower and upper limits for the parameter sigmasq") 
  if(length(tausq.rel) != 2)
    stop("tausq.rel must be a 2 components vector with lower and upper limits for the parameter tausq.rel") 
  if(length(kappa) != 2)
    stop("kappa must be a 2 components vector with lower and upper limits for the parameter kappa") 
  if(length(kappa2) != 2)
    stop("kappa must be a 2 components vector with lower and upper limits for the parameter kappa") 
  if(length(lambda) != 2)
    stop("lambda must be a 2 components vector with lower and upper limits for the parameter lambda")
  if(length(psiR) != 2)
    stop("psiR must be a 2 components vector with lower and upper limits for the parameter psiR") 
  if(length(psiA) != 2)
    stop("psiA must be a 2 components vector with lower and upper limits for the parameter psiA") 
  if(phi[1] >= phi[2])
    stop("parameter phi: lower limit greater or equal upper limit")
  if(sigmasq[1] >= sigmasq[2])
    stop("parameter sigmasq: lower limit greater or equal upper limit")
  if(tausq.rel[1] >= tausq.rel[2])
    stop("parameter tausq.rel: lower limit greater or equal upper limit")
  if(kappa[1] >= kappa[2])
    stop("parameter kappa: lower limit greater or equal upper limit")
  if(kappa2[1] >= kappa2[2])
    stop("parameter kappa: lower limit greater or equal upper limit")
  if(lambda[1] >= lambda[2])
    stop("parameter lambda: lower limit greater or equal upper limit")
  if(psiR[1] >= psiR[2])
    stop("parameter psiR: lower limit greater or equal upper limit")
  if(psiA[1] >= psiA[2])
    stop("parameter psiA: lower limit greater or equal upper limit")
  names(phi) <- names(sigmasq) <- names(tausq.rel) <- names(kappa) <- 
  names(kappa2) <- names(lambda) <- names(psiR) <- names(psiA) <- c("lower", "upper")
  return(list(phi = phi, sigmasq = sigmasq,
              tausq.rel = tausq.rel, kappa = kappa, kappa2 = kappa2,
              lambda = lambda, psiR = psiR, psiA = psiA))
}


"legend.krige" <-
  function(x.leg, y.leg, values, scale.vals, vertical = FALSE,
           offset.leg = 1, ...)
{
  values <- values[!is.na(values)]
  if(length(x.leg) != 2 | length(y.leg) != 2)
    stop("x.leg and y.leg require a vector with 2 elements")
  v.r <- range(values[is.finite(values)], na.rm = TRUE)
  lags.x <- function(xs, nl){
    xs.r <- 0.5 * diff(xs/(nl-1))
    return(seq(xs[1]+xs.r, xs[2]-xs.r, l=nl))
  }
  leg.l <- list(...)
  if(is.null(leg.l$br))
    nc <- ifelse(is.null(leg.l$col), 12, length(leg.l$col))
  else
    nc <- length(leg.l$breaks) - 1
  if(is.null(leg.l$col)) leg.l$col <- heat.colors(nc)
  if(is.null(leg.l$zl)) leg.l$zlim <- c(v.r[1], v.r[2])
  if(vertical){
    xy <- list(x=x.leg, y=lags.x(xs=y.leg, nl=nc))
    if(is.null(leg.l$br))
      image(x=xy$x, y=xy$y,
            z=matrix(seq(leg.l$zlim[1], leg.l$zlim[2], l=nc), nrow=1),
            add=TRUE, xaxs = "i", yaxs = "i", xlab="", ylab="",
            zlim = leg.l$zlim, col=leg.l$col)
    else
      image(x=xy$x, y=xy$y,
            z=matrix(seq(leg.l$zlim[1], leg.l$zlim[2], l=nc), nrow=1),
            add=TRUE, xaxs = "i", yaxs = "i", xlab="", ylab="",
            zlim = leg.l$zlim, col=leg.l$col, breaks = leg.l$br)
  }
  else{
    xy <- list(x=lags.x(xs=x.leg, nl=nc), y=y.leg)
    if(is.null(leg.l$br))
      image(x=xy$x, y=xy$y,
            z=matrix(seq(leg.l$zlim[1], leg.l$zlim[2], l=nc), ncol=1),
            add=TRUE, xaxs = "i", yaxs = "i", xlab="", ylab="",
            zlim = leg.l$zlim, col=leg.l$col)
    else
      image(x=xy$x, y=xy$y,
            z=matrix(seq(leg.l$zlim[1], leg.l$zlim[2], l=nc), ncol=1),
            add=TRUE, xaxs = "i", yaxs = "i", xlab="", ylab="",
            zlim = leg.l$zlim, col=leg.l$col, breaks = leg.l$br)
  }
  leg.poly <- rbind(c(x.leg[1], y.leg[1]), c(x.leg[2], y.leg[1]),
                    c(x.leg[2], y.leg[2]), c(x.leg[1], y.leg[2]),
                    c(x.leg[1], y.leg[1]))
  polygon(leg.poly)
#  if(is.null(leg.l$cex)) leg.l$cex <- par()$cex
  if(is.null(leg.l$cex)) leg.l$cex <- 0.8
  if(missing(scale.vals))
    scale.vals <- pretty(c(values,leg.l$zlim), n=5, min.n=4)
  scale.vals <- scale.vals[scale.vals > leg.l$zlim[1] &
                           scale.vals < leg.l$zlim[2]]
  if(vertical){
    y.r <- range(lags.x(xs=y.leg,nl=nc))
    y.text <- y.r[1] + ((scale.vals - leg.l$zlim[1]) * diff(y.r))/diff(leg.l$zlim)
    text((max(x.leg)+ offset.leg * diff(x.leg)), y.text,
         lab=scale.vals, col=1, cex=leg.l$cex)
  }
  else{
    x.r <- range(lags.x(xs=x.leg,nl=nc))
    x.text <- x.r[1] + ((scale.vals - leg.l$zlim[1]) * diff(x.r))/diff(leg.l$zlim)
    text(x.text, (max(y.leg)+ offset.leg * (diff(y.leg)/2)), lab=scale.vals, col=1, cex=leg.l$cex)
  }
  return(invisible())
}

"plot.1d" <-
  function(x, x1vals, ...)
{
  #cat("data in 1-D\n")
  if(length(x1vals) == 1) col.ind <- 2
  else col.ind <- 1
  order.it <- order(x$coords[,col.ind])
  pty.prev <- par()$pty
  par(pty="m")
  plot(x$coords[order.it,col.ind], x$data[order.it], ...)
  par(pty=pty.prev)
  return(invisible())
}

".ldots.set" <-
  function(ldots,
           type=c("persp", "image", "plot.1d", "contour", "filled.contour"),
           data = c("simulation", "prediction"))
{
  type <- match.arg(type)
  if(!is.null(ldots)){
#    all.dots <- c(names(formals(plot.default)),
#                  names(formals(title)),names(par()))
    fct <- switch(type,
                  persp = getS3method("persp", "default"),
                  image = image.default,
                  plot.1d = plot.default,
                  contour = contour.default,
                  filled.contour = filled.contour)
    ind <- pmatch(names(ldots), names(formals(fct)))
    #if(any(is.na(ind)))
#      ldots[is.na(ind)] <- NULL
#    if(!is.null(ldots))
#      names(ldots) <- names(formals(fct))[ind[!is.na(ind)]]
      ind <- pmatch(names(ldots), names(formals(fct)))
      names(ldots) <- ifelse(is.na(ind), names(ldots), names(formals(fct))[ind])
  }
  if(data == "simulation"){
    if(type == "plot.1d"){
      if(is.null(ldots$xlab)) ldots$xlab <- "x"
      if(is.null(ldots$ylab)) ldots$ylab <- "Y(x)"
    }
    else{
      if(is.null(ldots$xlab)) ldots$xlab <- "X Coord"
      if(is.null(ldots$ylab)) ldots$ylab <- "Y Coord"
      if(type == "persp")
        if(is.null(ldots$zlab)) ldots$zlab <- "Data"
    }
  }
  if(data == "prediction"){
    if(type == "plot.1d"){
      if(is.null(ldots$xlab)) ldots$xlab <- "x"
      if(is.null(ldots$ylab)) ldots$ylab <- "value"
    }
    else{
      if(is.null(ldots$xlab)) ldots$xlab <- "X Coord"
      if(is.null(ldots$ylab)) ldots$ylab <- "Y Coord"
      if(type == "persp")
        if(is.null(ldots$zlab)) ldots$zlab <- "value"
    }
  }
  if(any(type == c("image","contour","filled.contour")))
    if(is.null(ldots$asp)) ldots$asp=1
  return(ldots)
}

".prepare.graph.kriging" <-
  function (locations, borders, borders.obj=NULL, values, xlim, ylim, ...) 
{
  ind <- order(locations[, 2], locations[, 1])
  locations <- locations[ind, ]
  ## added 23/02/2009:
  if(nrow(locations) == length(borders)) values <- values[ind]
##  xx <- as.numeric(levels(as.factor(round(locations[, 1], digits = 8))))
  xx <- sort(unique(locations[, 1]))
  nx <- length(xx)
##  yy <- as.numeric(levels(as.factor(round(locations[, 2], digits = 8))))
  yy <- sort(unique(locations[, 2]))
  ny <- length(yy)
  ##
  ##  if(is.null(borders) && (nx*ny) > length(values))
  ##    borders <- locations[chull(locations),]
  ##
  values.loc <- rep(NA, nrow(locations))
  if(length(values.loc) == length(values))
    values.loc <- values
##    values.loc <- values[ind]
  if(!is.null(borders.obj)){
    borders.obj <- as.matrix(as.data.frame(borders.obj))
    dimnames(borders.obj) <- list(NULL, NULL)
    inout.vec <- as.vector(.geoR_inout(pts = locations,
                                       poly = borders.obj, ...))
    values.loc[inout.vec] <- values
    rm("inout.vec")
  }
  if (!is.null(borders)){
    borders <- as.matrix(as.data.frame(borders))
    dimnames(borders) <- list(NULL, NULL)
    if(!(!is.null(borders.obj) && identical(borders,borders.obj))){
      inout.vec <- as.vector(.geoR_inout(pts = locations,
                                         poly = borders, ...))
      if(length(values.loc[inout.vec]) == length(values))
        values.loc[inout.vec] <- values
      values.loc[!inout.vec] <- NA
      rm("inout.vec")
    }
  }
  ##
  if (missing(xlim) || is.null(xlim))
    if(is.null(borders)) xlim <- NULL
    else xlim <- range(borders[,1]) 
  if (missing(ylim) || is.null(ylim))
    if(is.null(borders)) ylim <- NULL
    else ylim <- range(borders[,2])
  coords.lims <- set.coords.lims(coords = locations,
                                 xlim = xlim, 
                                 ylim = ylim)
  coords.lims[, 1] <- coords.lims[, 1] + c(-0.025, 0.025) * 
    diff(coords.lims[, 1])
  coords.lims[, 2] <- coords.lims[, 2] + c(-0.025, 0.025) * 
    diff(coords.lims[, 2])
  return(list(x = xx, y = yy,
              values = matrix(values.loc, ncol = ny), 
              coords.lims = coords.lims))
}

".geoR_inout" <- function(pts, poly, ...) {
  ## inout returns logical vector
  ## poly <- .check.borders(poly)
  ## pts <- SpatialPoints(coords=pts)
  ## res <- over(pts, poly)
  ## !is.na(res)
  !is.na(sp::over(SpatialPoints(coords=pts), .check.borders(poly)))
}

".geoR_pip" <- function(pts, poly, ...) {
  ## pip returns the points matrix
  poly <- .check.borders(poly)
  pts <- sp::SpatialPoints(coords=pts)
  res <- sp::over(pts, poly)
  opts <- sp::coordinates(pts)[!is.na(res),]
  .check.locations(opts)
}

"pred_grid" <- function(coords, y.coords = NULL,
                        ... , y.by = NULL, y.length.out = NULL,
                        y.along.with = NULL){
  #, round.minmax=FALSE, digits){
  if(is.list(coords)){
    x.coords <- range(coords[[1]])
    y.coords <- range(coords[[2]])
  }
  else{
    if(is.matrix(coords) || is.data.frame(coords)){
      if(ncol(coords) != 2)
        stop("coords must be a two column matrix with xy-coordinates
or a vector with x-coordinates")
      x.coords <- range(coords[,1])
      y.coords <- range(coords[,2])
    }
    else{
      if(is.null(y.coords))
        stop("if a vector is provided in coords, y.coords must also be provided as a vector")
      x.coords <- range(coords)
      y.coords <- range(y.coords)
    }
  }
#  if(round.minmax){
#    if(missing(digits))
#      stop("argument \"digits\" must be provided if round.minmax=TRUE")
#    x.coords[1] <- x.coords[1]
#    x.coords[2] <- x.coords[2]
#    y.coords[1] <- y.coords[1]
#    y.coords[2] <- y.coords[2]
#  }
  gx <- seq(x.coords[1], x.coords[2], ...)
  ldots <- list(...)
  if(!is.null(ldots))
    names(ldots) <- sapply(names(ldots), function(x) match.arg(x, names(formals(seq.default))))
  if(is.null(y.length.out)) y.length.out <- ldots$length.out
  if(is.null(y.along.with)) y.along.with <- ldots$along.with
  if(is.null(y.by))
    y.by <- ifelse(is.null(ldots$by), ((y.coords[2] - y.coords[1])/y.length.out - 1), ldots$by)
  gy <- seq(from=y.coords[1], to=y.coords[2], by = y.by)
  res <- expand.grid(gx,gy)
  attr(res, "x") <- gx
  attr(res, "y") <- gy
  return(res)
}

