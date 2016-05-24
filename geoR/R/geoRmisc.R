##
## Miscelaneous geoR functions
##

"jitterDupCoords" <-
  function(x, ...)
{
  UseMethod("jitterDupCoords")
}

"jitterDupCoords.default" <-
  function(x, ...){
    ## checking input
    if(!is.matrix(x) && !is.data.frame(x) || ncol(x) != 2)
      stop("jitterCoords.default: coords must be a matrix or data-frame with two columns")
    ind <- dup.coords(x, simplify=FALSE, USE.NAMES=FALSE)
    i <- 1
    while(i <= length(ind)){
      x[ind[[i]],] <- jitter2d(coords=x[ind[[i]],], ...)
      i <- i+1
    }
    return(x)
  }

"jitterDupCoords.geodata" <-
  function(x, ...){
    x$jitter.Random.seed <- .Random.seed
    x$coords <- jitterDupCoords.default(x=x$coords, ...)
    return(x)
  }


"jitter2d" <-
  function(coords, max, min = 0.2*max,
           fix.one = TRUE, which.fix=c("random", "first", "last")){
    ## checking input
    if(!is.matrix(coords) && !is.data.frame(coords) || ncol(coords) != 2)
      stop("jitterCoords: coords must be a matrix or data-frame with two columns")
    nc <- nrow(coords)
    if(mode(which.fix) != "numeric"){
      which.fix <- match.arg(which.fix)
      which.fix <- switch(which.fix,
                          random = sample(1:nc, 1),
                          first = 1,
                          last = nc)
    }
    if(length(which.fix) > 1)
      stop("jitterCoords: which.fix must be a single element vector")
    ## number of coordinates to be jittered
    if(fix.one) nc <- nc-1
    ## random displacements of the coordinates
    angle <- runif(nc, min = 0, max = 2 * pi)
    d <- sqrt(runif(nc, min = min^2, max = max^2))
    coords[-which.fix,] <- coords[-which.fix,] + cbind(d * cos(angle), d * sin(angle))
    return(coords)
  }



"globalvar" <- function(geodata, locations, coords = geodata$coords, krige)
{
  ##
  ## computes the variance of a global mean (as if a block size was equal
  ## to the study area) cf. Isaaks & Srivastava, pag. 508
  ##
  if(ncol(coords) != 2)
    stop("function only available for 2-D coordinates")
  if(ncol(locations) != 2)
    stop("function only available for 2-D locations")
  ##
  ## reading input
  ##
  if(missing(krige))
    krige <- krige.control()
  else{
    ##    if(is.null(class(krige)) || class(krige) != "krige.geoR"){
    if(length(class(krige)) == 0 || class(krige) != "krige.geoR"){
      if(!is.list(krige))
        stop("krige.conv: the argument krige only takes a list or an output of the function krige.control")
      else{
        krige.names <- c("type.krige","trend.d","trend.l","obj.model",
                         "beta","cov.model", "cov.pars",
                         "kappa","nugget","micro.scale","dist.epsilon",
                         "lambda","aniso.pars")
        krige.user <- krige
        krige <- list()
        if(length(krige.user) > 0){
          for(i in 1:length(krige.user)){
            n.match <- match.arg(names(krige.user)[i], krige.names)
            krige[[n.match]] <- krige.user[[i]]
          }
        }
        if(is.null(krige$obj.model)) krige$obj.model <-  NULL
        if(is.null(krige$cov.model)) krige$cov.model <- "matern"  
        if(is.null(krige$cov.pars))
          stop("covariance parameters (sigmasq and phi) should be provided in cov.pars")
        if(is.null(krige$kappa)) krige$kappa <-  0.5
        if(is.null(krige$nugget)) krige$nugget <-  0
        if(is.null(krige$micro.scale)) krige$micro.scale <- 0  
        if(is.null(krige$dist.epsilon)) krige$dist.epsilon <-  1e-10
        if(is.null(krige$aniso.pars)) krige$aniso.pars <- NULL  
        if(is.null(krige$lambda)) krige$lambda <- 1 
        krige <- krige.control(obj.model = krige$obj.model,
                               cov.model = krige$cov.model,
                               cov.pars = krige$cov.pars,
                               kappa = krige$kappa,
                               nugget = krige$nugget,
                               micro.scale = krige$micro.scale,
                               dist.epsilon = krige$dist.epsilon, 
                               aniso.pars = krige$aniso.pars,
                               lambda = krige$lambda)   
      }
    }
  }
  ## mean of the values in the prediction grid covariance matrix
  krige$coords <- locations
  CAA <- mean(do.call("varcov.spatial", krige)$varcov)
  ## \Sigma^{-1/2}
  krige$coords <- coords
  krige$inv <- TRUE
  invS <- do.call("varcov.spatial", krige)$inverse
  
  ## mean of the data/location covariance matrix
  CiA <- rowMeans(cov.spatial(obj = loccoords(coords=coords,
                             locations=locations),
                           cov.model = krige$cov.model,
                           cov.pars=krige$cov.pars,
                           kappa=krige$kappa))
  ##
  oneSinvone <- sum(invS)
  w <- colSums(invS)/oneSinvone
  return(CAA + (1/oneSinvone) - 2 * sum(w * CiA))
}

"nearloc" <- function(points, locations, positions=FALSE)
{
  if(!is.data.frame(points) && is.list(points)){
    points <- matrix(unlist(points[1:2]), ncol=2)
    backtolist <- TRUE
    backnames <- names(points[1:2])
  }
  else backtolist <- FALSE
  if(is.matrix(locations))
    locations <- locations[,1:2]
  if(is.data.frame(locations))
    locations <- as.matrix(locations[,1:2])
  if(is.list(locations))
    locations <- matrix(unlist(locations[1:2]), ncol=2)
  if(ncol(locations) != 2)
    stop("locations must be a matrix, data.frame, or list")
  ind <- apply(loccoords(points, locations), 1, which.min)
  if(positions)
    return(ind)
  res <- locations[apply(loccoords(points, locations), 1, which.min),]
  rownames(res) <- ind
  if(backtolist){
    res <- list(res[,1], res[,2])
    names(res) <- backnames
  }
  return(res)
}


"coords2coords" <-
  function(coords, xlim, ylim, xlim.ori, ylim.ori)
{
  if(missing(ylim.ori)) xlim.ori <- range(coords[,1], na.rm=TRUE)
  if(missing(ylim.ori)) ylim.ori <- range(coords[,2], na.rm=TRUE)
  coords[,1] <- xlim[1] + (coords[,1] - xlim.ori[1]) * diff(xlim)/diff(xlim.ori)
  coords[,2] <- ylim[1] + (coords[,2] - ylim.ori[1]) * diff(ylim)/diff(ylim.ori)
  return(coords)
}

"zoom.coords" <-
  function(x, ...)
{
  UseMethod("zoom.coords")
}

"zoom.coords.default" <-
    function(x, xzoom, yzoom=xzoom, xlim.ori, ylim.ori, xoff=0, yoff=0, ...)
{
  if(missing(ylim.ori)) xlim.ori <- range(x[,1], na.rm=TRUE)
  if(missing(ylim.ori)) ylim.ori <- range(x[,2], na.rm=TRUE)
  xlim <- xlim.ori + c(-1,1) * (diff(xlim.ori)/2) * (xzoom - 1)
  ylim <- ylim.ori + c(-1,1) * (diff(ylim.ori)/2) * (yzoom - 1)
  res <- coords2coords(x, xlim=xlim, ylim=ylim, xlim.ori = xlim.ori, ylim.ori=ylim.ori)
  res[,1] <- res[,1] + xoff
  res[,2] <- res[,2] + yoff
  return(res)
}

"zoom.coords.geodata" <-
  function(x, ...)
{
  x$coords <- zoom.coords.default(x$coords, ...)
  if(!is.null(x$borders)) x$borders <- zoom.coords.default(x$borders, ...)
  if(!is.null(x$subarea.lims)) x$subarea.lims <- zoom.coords.default(x$subarea.lims, ...)
  return(x)
}

"rect.coords" <-
  function(coords, xzoom = 1, yzoom=xzoom, add.to.plot=TRUE, quiet=FALSE, ...)
{
  rx <- range(coords[,1], na.rm=TRUE)
  ry <- range(coords[,2], na.rm=TRUE)
  res <- cbind(c(rx,rev(rx)), rep(ry,c(2,2)))
  res <- zoom.coords(res, xzoom=xzoom, yzoom=yzoom)
  if(add.to.plot) rect(res[1,1], res[1,2], res[3,1], res[3,2], ...)
  if(quiet)
    return(invisible())
  else
    return(res)  
}

"dup.coords" <-
  function(x, ...)
{
  UseMethod("dup.coords")
}

"dup.coords.default" <-
  function(x, ...)
{
  ap1 <- unclass(factor(paste("x",x[,1],"y",x[,2], sep="")))
  ap2 <- table(ap1)
  ap2 <- ap2[ap2 > 1]
  takecoords <- function(n){
    if(!is.null(rownames(x))) rownames(x[ap1 == n,])
    else (1:nrow(x))[ap1 == n]
    }
  res <- sapply(as.numeric(names(ap2)), takecoords, ...)
  if(length(res) == 0) res <- NULL
  if(!is.null(res)) class(res) <- "duplicated.coords"
  return(res)
}

"dup.coords.geodata" <- "duplicated.geodata" <-
  function(x, incomparables, ...)
{
  xdf <- as.data.frame.geodata(x)
  dc <- dup.coords.default(x$coords)
  if(is.null(dc)) return(dc)
  else{
    return(data.frame(dup=factor(rep(1:length(dc), sapply(dc, length))),
                      xdf[unlist(dc),]))
  }
}


"subarea" <-
  function(geodata, xlim, ylim, ...)
{
  if(class(geodata) != "geodata")
    stop("an object of the class geodata must be provided")
  if(missing(xlim) & !missing(ylim)) xlim <- c(-Inf, +Inf)
  if(!missing(xlim) & missing(ylim)) ylim <- c(-Inf, +Inf)
  if(missing(xlim) & missing(ylim)){
    cat("Enter 2 points defining the corners of the subarea\n")
    pt <- locator(2)
    xlim <- sort(pt[[1]])
    ylim <- sort(pt[[2]])
  }
  if(!is.vector(xlim) || length(xlim) != 2)
    stop("xlim must be a vector with 2 elements")
  if(!is.vector(ylim) || length(ylim) != 2)
    stop("ylim must be a vector with 2 elements")
  geo <- geodata
  ind <- (geodata$coords[,1] > xlim[1] & geodata$coords[,1] < xlim[2] &  
          geodata$coords[,2] > ylim[1] & geodata$coords[,2] < ylim[2])  
  geo$coords <- geodata$coords[ind,]
  xlim.all <- c(xlim, range(geo$coords[,1]))
  ylim.all <- c(ylim, range(geo$coords[,2]))
  if(is.vector(geodata$data)) geo$data <- geodata$data[ind]
  else geo$data <- geodata$data[ind,]
  geo$units.m <- geodata$units.m[ind]
  if(!is.null(geodata$covariate)){
    if(is.vector(geodata$covariate))
      geo$covariate <- geodata$covariate[ind]
    if(is.matrix(geodata$covariate) |  is.data.frame(geodata$covariate))   
      geo$covariate <- geodata$covariate[ind,]
  }
  if(!is.null(geodata$borders)){
    geo$borders <- geodata$borders[(geodata$borders[,1]>xlim[1] & geodata$borders[,1]<xlim[2] &  
                                    geodata$borders[,2]>ylim[1] & geodata$borders[,2] < ylim[2]),]
    xlim.all <- c(xlim.all, range(geo$borders[,1]))
    ylim.all <- c(ylim.all, range(geo$borders[,2]))
  }
  geo$subarea.lims <- cbind(range(xlim.all[is.finite(xlim.all)]),
                            range(ylim.all[is.finite(ylim.all)]))
  return(geo)
}

