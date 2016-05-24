##
## "wrappers" for pieces of C code in geoR package
## -----------------------------------------------
##
".Ccor.spatial" <- function(x, phi, kappa, cov.model)
{
  res <- rep(0.0, length(x))
  .C("veccorrval",
     as.double(phi),
     as.double(kappa),
     as.double(x),
     as.integer(length(x)),
     cornr = as.integer(.cor.number(cov.model)),
     out = as.double(res), PACKAGE = "geoR")$out
}

".bilinearformXAY" <-
  function(X, lowerA, diagA, Y)
  {
    nA <- length(diagA)
    nX <- length(X)/nA
    nY <- length(Y)/nA
    if(length(lowerA) != (nA * (nA -1)/2))
      stop("lowerA and diagA have incompatible dimensions")
    #temp.na <- c(any(is.na(lowerA)),any(is.na(diagA)),any(is.na(X)),
    #          any(is.na(Y)),any(is.na(nX)),any(is.na(nY)),
    #          any(is.na(nA)))
    #names(temp.na) <- c("lowerA","diagA","X","Y","nX","nY","nA")
    #print(temp.na)
    #temp.nan <- c(any(is.nan(lowerA)),any(is.nan(diagA)),any(is.nan(X)),
    #          any(is.nan(Y)),any(is.nan(nX)),any(is.nan(nY)),
    #          any(is.nan(nA)))
    #names(temp.nan) <- c("lowerA","diagA","X","Y","nX","nY","nA")
    #print(temp.nan)
    #temp.inf <- c(any(!is.finite(lowerA)),any(!is.finite(diagA)),any(!is.finite(X)),
    #          any(!is.finite(Y)),any(!is.finite(nX)),any(!is.finite(nY)),
    #          any(!is.finite(nA)))
    #names(temp.inf) <- c("lowerA","diagA","X","Y","nX","nY","nA")
    #print(temp.inf)
    out <- .C("bilinearform_XAY",
              as.double(as.vector(lowerA)),
              as.double(as.vector(diagA)),
              as.double(as.vector(X)),
              as.double(as.vector(Y)),
              as.integer(nX),
              as.integer(nY),
              as.integer(nA),
              res=as.double(rep(0,(nX*nY))),
              PACKAGE = "geoR")$res
    attr(out, "dim") <- c(nX, nY)
    return(out)
  }

".diagquadraticformXAX" <-
  function(X, lowerA, diagA)
  {
    nA <- length(diagA)
    nX <- length(X)/nA
    if(length(lowerA) != (nA * (nA -1)/2))
      stop("lowerA and diagA have incompatible dimensions")
    out <- .C("diag_quadraticform_XAX",
              as.double(as.vector(lowerA)),
              as.double(as.vector(diagA)),
              as.double(as.vector(X)),
              as.integer(nX),
              as.integer(nA),
              res = as.double(rep(0,nX)),
              PACKAGE = "geoR")$res
    return(out)
  }

"loccoords" <-
  function(coords, locations)
  {
    ## Computes a matrix for which each row has the distances between
    ## each point in 'locations' to all the points in 'coords'
    coords <- as.matrix(coords)
    locations <- as.matrix(locations)
    dimc <- dim(coords)
    diml <- dim(locations)
    if((dimc[2] != 2) | (diml[2] != 2))
      stop("coords and locations must have two columns")
    nc <- dimc[1]
    nl <- diml[1]
    out <- .C("loccoords",
       as.double(as.vector(locations[,1])),
       as.double(as.vector(locations[,2])),
       as.double(as.vector(coords[,1])),
       as.double(as.vector(coords[,2])),
       as.integer(nl),
       as.integer(nc),
       out = as.double(rep(0, nc*nl)),
       PACKAGE = "geoR")$out
    dim(out) <- c(nc, nl)
    return(out)
  }

"diffpairs" <-
  function(coords, data)
  {
    ## returns:
    ##   - the lower triangle of the matrix with euclidean distances
    ##     between pairs of points,
    ##   - the difference between data values at these locations
    ##
    coords <- as.matrix(coords)
    data <- as.vector(data)
    if(length(data) != nrow(coords)) stop('incompatible dimensions between data and coords')
    dimc <- dim(coords)
    if(dimc[2] == 1 & dimc[1] == 2) return(0)
    else{
      if(dimc[2] != 2)
        stop("coords must have two columns")
      nc <- dimc[1]
      out <- as.double(rep(0, (nc * (nc-1))/2))
      res <- .C("diffpairs",
                as.double(coords[,1]),
                as.double(coords[,2]),
                as.double(data),
                as.integer(nc),
                dist = out, diff = out,
                PACKAGE = "geoR")[c("dist","diff")]
      return(res)
    }
  }

".corr.diaglowertri" <-
  function(coords, cov.model, phi, kappa)
{
  cov.model <- match.arg(cov.model, choices = .geoR.cov.models)
  if(any(cov.model == c("cauchy", "gencauchy", "matern", "powered.exponential", "power", "gneiting.matern"))){
    if(missing(kappa))
      stop("argument kappa is needed for this choice of correlation function")
  }
  else kappa <- 1
  coords <- as.matrix(coords)
  dimc <- dim(coords)
  if(dimc[2] == 1 & dimc[1] == 2)
    return(0)
  else{
    if(dimc[2] != 2)
      stop("coords must have two columns")
    nc <- dimc[1]
    .C("cor_diag",
       as.double(coords[,1]),
       as.double(coords[,2]),
       as.integer(nc),
       as.integer(.cor.number(cov.model)),
       as.double(phi),
       as.double(kappa),
       out = as.double(rep(0, (nc * (nc+1)/2))), PACKAGE = "geoR")$out
  }
}

".cond.sim" <-
  function(env.loc, env.iter, loc.coincide, coincide.cond, tmean,
           Rinv, mod, vbetai, fixed.sigmasq)
{
  NTOT <- mod$nloc * mod$Nsims
  if(fixed.sigmasq)
    invchisc <- rep(1, NTOT)
  else
    invchisc <- sqrt(mod$df.model/rchisq(mod$Nsims, df=mod$df.model))
  ##
  if(mod$beta.size == 1){
    Blower <- 0
    Bdiag <- vbetai
  }
  else{
    Blower <- vbetai[lower.tri(vbetai)]
    Bdiag <- diag(vbetai)
  }
  ##
  if((length(tmean) %% mod$nloc) > 0)
    stop(".cond.sim: wrong size of tmean")
  tmean <- matrix(tmean, nrow = mod$nloc)
  ncol.tmean <- ncol(tmean)
  if(ncol(tmean) > 1){
    if(ncol.tmean != mod$Nsims)
      stop(".cond.sim: size of tmean does not matches with Nsims")

    diff.mean <- as.integer(1)
  }
  else
    diff.mean <- as.integer(0)
  if(coincide.cond) loccoin <- -loc.coincide
  else loccoin <- TRUE
  if(mod$phi < 1e-24 && mod$s2 < 1e-24)
    normalsc <- rnorm(NTOT, mean=tmean, sd=sqrt(mod$nugget+mod$s2))
  else{
    normalsc <- rnorm(NTOT)
    normalsc <- .C("kb_sim_new",
                   as.double(as.vector(tmean)),
                   out = as.double(normalsc),
                   as.double(as.vector(Rinv$lower)),
                   as.double(as.vector(Rinv$diag)),
                   as.double(as.vector(get("v0", envir=env.iter))),
                   as.integer(mod$nloc),
                   as.integer(mod$n),
                   as.double(mod$Dval),
                   as.integer(mod$Nsims),
                   as.double(invchisc),
                   as.double(mod$s2),
                   as.double(Blower),
                   as.double(Bdiag),
                   as.double(as.vector(get("b", envir=env.iter))),
                   as.integer(mod$beta.size),
                   as.double(get("locations", envir=env.loc)[loccoin,1]),
                   as.double(get("locations", envir=env.loc)[loccoin,2]),
                   as.integer(mod$cov.model.number),
                   as.double(mod$phi),
                   as.double(mod$kappa),
                   as.integer(diff.mean),
                   PACKAGE = "geoR")$out
  }
  attr(normalsc, "dim") <- c(mod$nloc, mod$Nsims)
  return(normalsc)
}

