##
## Basic data manipulation for the geoR package
## ---------------------------------------------------------------
##
## Functions for reading data and basic exploratory analysis. 
## These functions include
##    - creating objets of the (S3) class "geodata"
##    - (S3) methods for this class
##

"names.geodata" <-
  function(x)
{
  res <- list(coords = colnames(x$coords))
  if(is.vector(x$data)) res$data <- "data"
  else res$data <- colnames(x$data)
  if(!is.null(x$units.m)) res$units.m <- "units.m"
  if(!is.null(x$covariate)) res$covariate <- colnames(x$covariate)
  if(!is.null(x$realisations)) res$realisations <- "realisations"
  if(!is.null(x$borders)) res$borders <- "borders"
  res$other <- names(unclass(x))[!(names(unclass(x)) %in% c("coords","data","units.m","covariate","realisations","borders"))]
  if(length(res$other) == 0) res$other <- NULL
  return(res)
}


"subset.geodata" <-
  function (x, ..., other = TRUE) 
{
  xdf <- as.data.frame(x)
  attribs <- list(ncol.data = attributes(xdf)$ncol.data,
                            ncol.covariate = attributes(xdf)$ncol.covariate,
                            ncol.units.m = attributes(xdf)$ncol.units.m,
                            ncol.realisations = attributes(xdf)$ncol.realisations)
  xdf <- subset.data.frame(xdf, ...)
  attributes(xdf) <- c(attributes(xdf), attribs)
  xdf <- as.geodata(xdf)
  if (other) 
    xdf <- c(xdf, x[names(unclass(x))[!(names(unclass(x)) %in% 
                                        c("coords", "data", "units.m", "covariate", "realisations"))]])
  oldClass(xdf) <- "geodata"
  return(xdf)
}

"as.data.frame.geodata" <-
  function (x, ..., borders = TRUE) 
{
  if (is.vector(x$data)) 
    xdf <- data.frame(x$coords, data = x$data)
  else xdf <- data.frame(x$coords, x$data)
  nc.data <- ncol(xdf)
  nc.units <- nc.covariate <- nc.realisations <- NULL
  if (!is.null(x$units.m)) {
    xdf$units.m <- x$units.m
    nc.units <- ncol(xdf)
  }
  if (!is.null(x$realisation)){
    xdf$realisations <- x$realisations
    nc.realisations <- ncol(xdf)
  }
  if (!is.null(x$covariate)){
    nc <- ncol(xdf)+1
    xdf <- data.frame(xdf, x$covariate)
    nc.covariate <- nc:ncol(xdf)
  }
  ## setting columns ID's as attributes
  attr(xdf, "ncol.data") <- 3:nc.data
  attr(xdf, "ncol.units.m") <- nc.units
  attr(xdf, "ncol.realisations") <- nc.realisations
  attr(xdf, "ncol.covariate") <- nc.covariate
  ##
  if (borders && !is.null(x$borders)) 
    attr(xdf, "borders") <- x$borders
  oldClass(xdf) <- c("geodata.frame", "data.frame")
  return(xdf)
}

"is.geodata" <- function (x){
  inherits(x, "geodata") && (!is.null(x$coords)) &&
  (!is.null(x$data)) && ncol(x$coords) == 2 &&
  (nrow(x$coords) == nrow(as.matrix(x$data)))
}

"read.geodata" <-
  function(file, header = FALSE, coords.col= 1:2, data.col = 3,
           data.names = NULL, covar.col = NULL,
           covar.names = "header",
           units.m.col = NULL, realisations = NULL,
           na.action = c("ifany", "ifdata", "ifcovar", "none"),
           rep.data.action, rep.covar.action, rep.units.action, ...)
{
  call.fc <- match.call()
  ##
  obj <- read.table(file = file, header = header, ...)
  if(all(covar.names == "header")){
    if(!is.null(covar.col)){
      col.names <- names(obj)
      covar.names <- col.names[covar.col]
    }
    else covar.names <- NULL
  }
  ##
  if(missing(rep.data.action)) rep.data.action <- "none"
  if(!is.function(rep.data.action))
    rep.data.action <- match.arg(rep.data.action, choices = c("none", "first")) 
  if(missing(rep.covar.action)) rep.covar.action <- rep.data.action
  if(!is.function(rep.covar.action))
    rep.covar.action <- match.arg(rep.covar.action, choices = c("none", "first")) 
  ##
  res <- as.geodata(obj = obj, coords.col = coords.col, data.col = data.col,
                    covar.col = covar.col, covar.names = covar.names,
                    realisations = realisations,
                    rep.data.action = rep.data.action,
                    rep.covar.action = rep.covar.action)
  res$call <- call.fc
  return(res)
}

"as.geodata" <-
  function (obj, ...)
  UseMethod("as.geodata")

"as.geodata.geodata.frame" <-
  function(obj, ...)
{
  res <- as.geodata.default(obj, coords.col = 1:2,
                            data.col = attributes(obj)$ncol.data,
                            units.m.col = obj$units.m,
                            covar.col = attributes(obj)$ncol.covariate,
                            realisations = obj$realisation)
  res$borders <- attributes(obj)$borders
  return(res)
}

"as.geodata.SpatialPointsDataFrame" <-
  function(obj, data.col = 1, ...)
{
  df <- as.data.frame(coordinates(obj))
  if(ncol(df) > 2){
    if(is.null(list(...)$coords.col))
      stop("geoR only works with 2D coordinates. Pass extra argument coords.col to specify the positions of first and second coordinate")
    else
      df <- df[,list(...)$coords.col]
  }
  df <- cbind(as.data.frame(obj@data), df)
  nc <- ncol(df)    
  return(as.geodata.default(df, coords.col = (nc-1):nc,
                            data.col = data.col, ...)) 
}

"as.geodata.default" <-
  function(obj, coords.col = 1:2, data.col = 3, data.names = NULL, 
           covar.col = NULL, covar.names = "obj.names",
           units.m.col = NULL, realisations = NULL,
           na.action = c("ifany", "ifdata", "ifcovar", "none"),
           rep.data.action, rep.covar.action, rep.units.action, ...)
{
  ## converts a simulation generated by grf
  if(class(obj)[1] == "grf"){
    res <- list(coords=obj$coords, data=obj$data)
    oldClass(res) <- "geodata"
    return(res)
  }
  ## checking input
  if(!is.matrix(obj) & !is.data.frame(obj))
    stop("object must be a matrix or data.frame")
  if(!is.null(data.names) & length(data.col) < 2)
    stop("data.names allowed only if there is more than 1 column of data")
  res <- list()
  ##
  ## testing for NA's setting the coordinates of the data locations
  ##
  if(any(is.na(obj[,coords.col]))){
    warning("NA's not allowed in the coordinates")
    obj <- obj[complete.cases(obj),,drop = FALSE]
    warning("eliminating rows with NA's")
  }
  res$coords <- as.matrix(obj[,coords.col])
  if(any(!is.numeric(res$coords)))
    stop("non-numerical values in the coordinates")
  if(is.null(colnames(res$coords)))
     colnames(res$coords) <- paste("Coord", c("1","2"), sep="")
  ##
  ## setting the data
  ##
  res$data <- as.matrix(obj[,data.col])
  if(length(data.col) == 1) res$data <- as.vector(res$data)
  else if(!is.null(data.names)) colnames(res$data) <- data.names
  ##
  ## adding optional elements to the list
  ##
  pos <- 3
  ##
  ## setting the covariates, if the case 
  ##
  if(!is.null(covar.col)){
    pos.covar <- pos
    res[[pos]] <- as.data.frame(obj[,covar.col])
    if(mode(covar.col) == "character") covar.names <- covar.col
    if(all(covar.names == "obj.names")){
      if(is.matrix(obj))      
        col.names <- dimnames(obj)[2]
      if(is.data.frame(obj))      
        col.names <- names(obj)
    }
    names(res)[pos] <- "covariate"
    if(all(covar.names == "obj.names"))
      if(is.null(col.names))
        names(res[[pos]]) <- paste("covar", 1:length(covar.col), sep="")
      else  names(res[[pos]]) <- col.names[covar.col]
    else
      names(res[[pos]]) <- covar.names
    covar.names <- names(res[[pos]])
    pos <- pos + 1
  }
  ##
  ## setting units.m, if the case 
  ##
  if(!is.null(units.m.col)){
    pos.units <- pos
    if(length(units.m.col) > 1){
      if(length(units.m.col) != nrow(res$coords))
        stop("units.m.col must be of length 1 or the same length as the coordinates")
      else res[[pos]] <- units.m.col
    }
    else
      res[[pos]] <- obj[,units.m.col,drop=TRUE]
    if(!all(res[[pos]] > 0))
      stop("all values of units.m must be greater than zero")
    names(res)[pos] <- "units.m"
  }
  ##
  ## Dealing with NA's
  ##
  na.action <- match.arg(na.action)
  if(na.action != "none"){
    if(na.action == "ifany")
      na.data <- na.covar <- TRUE
    if(na.action == "ifdata")
      {na.data <- TRUE; na.covar <- FALSE}
    if(na.action == "ifcovar")
      {na.data <- FALSE; na.covar <- TRUE}
    not.na <- function(x) !any(is.na(x))
    if(na.data){
      ind <- apply(as.matrix(res$data), 1, not.na)
      if(!is.null(units.m.col)) ind <- (ind & !is.na(res$units.m))
      if(!all(ind)){
        res$coords <- res$coords[ind,]
        res$data <- drop(as.matrix(res$data)[ind,])
        if(!is.null(covar.col))
          res[[pos.covar]] <- drop(as.matrix(res[[pos.covar]][ind,]))
        cat(paste("as.geodata:", sum(!ind),
  "points removed due to NA in the data\n"))
        if(!is.null(units.m.col))
          res$units.m <- res$units.m[ind]
      }
    }
    if(!is.null(covar.col) && na.covar){
      ind <- apply(as.matrix(res[[pos.covar]]), 1, not.na)
      if(!all(ind)){
        res$coords <- res$coords[ind,]
        res$data <- drop(as.matrix(res$data)[ind,])
        if(!is.null(covar.col))
          res[[pos.covar]] <- drop(res[[pos.covar]][ind,])
        cat(paste("as.geodata:", sum(!ind), "points removed due to NA in the covariate(s)\n")) 
        if(!is.null(units.m.col))
          res$units.m <- res$units.m[ind]
      }
    }
  }
  ##
  ## Checking whether there are data from different realisations
  ##
  if(is.null(realisations))
    realisations <- as.factor(rep(1, nrow(res$coords)))
  else{
    if(mode(realisations) == "numeric" && length(realisations) == 1)
      realisations <- as.factor(obj[,realisations])
    res$realisations <- realisations
  }
  if(length(realisations) != nrow(res$coords))
    stop("realisations and coords have incompatible dimensions")
  ##
  ## Checking whether there are data at coincident locations
  ## and dealing with this according to the value of the argument
  ## rep.data.action 
  ##
  if(missing(rep.data.action)) rep.data.action <- "none"
  if(!is.function(rep.data.action))
    rep.data.action <- match.arg(rep.data.action,
                                 choices = c("none", "first")) 
  if(missing(rep.covar.action)) rep.covar.action <- rep.data.action
  if(!is.function(rep.covar.action))
    rep.covar.action <- match.arg(rep.covar.action,
                                  choices = c("none", "first")) 
  if(missing(rep.units.action)) rep.units.action <- rep.data.action
  if(!is.function(rep.units.action))
    rep.units.action <- match.arg(rep.units.action,
                                  choices = c("none", "first")) 
  ## check also whether this should be checked within each realisation 
  ##  if(sum(rep.dup) > 0){
  if(any(duplicated(res$coords, MAR=1))){
    op.dig <- options()$digits
    options(digits=12)
    rep.lev <- as.character(paste("x",res$coords[,1],"y",
                                  res$coords[,2], sep=""))
    options(digits=op.dig)
    rep.dup <- duplicated(rep.lev)
    cat(paste("as.geodata:", sum(rep.dup), "replicated data locations found. \n Consider using jitterDupCoords() for jittering replicated locations. \n"))
    if(is.function(rep.data.action) || rep.data.action == "first"){
      res$coords <- res$coords[!rep.dup,]
      measure.var.f <- function(x)
        return(summary(lm(x ~ as.factor(rep.lev)))$sigma^2)
      res$m.var <- drop(apply(as.matrix(res$data),2,measure.var.f))
      rep.action.f <- function(x, rep.action){
        if(!is.function(rep.action) && rep.action == "first")
          return(x[!rep.dup])
        else
          return((as.vector(by(x, rep.lev, rep.action))[unclass(factor(rep.lev))])[!rep.dup])
      }
      res$data <- drop(apply(as.matrix(res$data), 2, rep.action.f, rep.action=rep.data.action))
      if(!is.null(covar.col))
        res[[pos.covar]] <- drop(apply(res[[pos.covar]], 2, rep.action.f, rep.action=rep.covar.action))
      if(!is.null(units.m.col))
        res$units.m <- drop(apply(as.matrix(res$units.m), 2 , rep.action.f, rep.action=rep.units.action))
      if(!is.null(res$realisations))
        res$realisations <- res$realisations[!rep.dup]
    }
    else{
      any.coincide <- as.vector(unlist(by(as.data.frame(res$coords), list(realisations),
                                function(x) any(duplicated(x)))))
     ## check.coincide <- function(x){sum(dist(x) < 1e-16) > 0}
     ## any.coincide <- lapply(split(as.data.frame(res$coords), realisations), check.coincide)
     ## any.coincide <- as.vector(unlist(any.coincide))
      if(sum(any.coincide) > 0)
        cat("WARNING: there are data at coincident or very closed locations, some of the geoR's functions may not work.\n Use function dup.coords() to locate duplicated coordinates.\n Consider using jitterDupCoords() for jittering replicated locations \n")      
    }
  }
  ##
  if(!is.null(covar.col)){
    res[[pos.covar]] <- as.data.frame(res[[pos.covar]])
    names(res[[pos.covar]]) <- covar.names
  }
  oldClass(res) <- "geodata"
  return(res)
}

##"summary.geodata" <- function(object, trend="cte", lambda=1,by.realisations=TRUE, ...)
"summary.geodata" <-
  function(object, lambda=1, add.to.data = 0, by.realisations=TRUE,  ...)
{
  res <- list()
  ##
  ## data transformation (Box-Cox)
  ##
  if (lambda != 1) Tdata <- BCtransform(x=object$data, lambda=lambda,
        add.to.data = add.to.data)$data
  ##
  ## trend "removal" (getting residuals from a linear model)
  ##
#  xmat <- unclass(trend.spatial(trend = trend, geodata = x))
#  if (nrow(xmat) != nrow(coords)) 
#    stop("coords and trend have incompatible sizes")
#  if (trend != "cte") {
#    data <- lm(data ~ xmat + 0)$residuals
#    names(data) <- NULL
#  }
  if(!is.null(object$realisations) && by.realisations && length(unique(object$realisations)) > 1){
    by2matrix <- function(x){
      x <- unclass(x)
      attributes(x) <- NULL
      x.names <- names(x[[1]])
      x.n <- length(x)
      x <- matrix(unlist(x), nrow=x.n, byrow=TRUE)
      colnames(x) <- x.names
      rownames(x) <- real.names
      return(x)
    }
    real <- as.factor(object$realisations)
    real.names <- paste("realisation.", levels(real), sep="")
    res$n <- as.vector(by(object$coords, real, nrow))
    names(res$n) <- real.names
    res$coords.summary <- by(object$coords, real,
                             function(x) {res <- apply(x,2,range);
                                          rownames(res) <- c("min", "max");
                                          if(is.null(colnames(res))) colnames(res) <- c("Coord.X", "Coord.Y"); 
                                          res})
    attr(res$coords.summary, "dimnames") <- list(realisation = levels(real)) 
    res$distances.summary <- by2matrix(by(object$coords, real,
                                function(x){ res <- range(dist(x));
                                             names(res) <- c("min", "max");  
                                             res}))
    res$data.summary <- by(object$data, real,
                           function(x) drop(apply(as.matrix(x), 2, summary)))
    if(ncol(as.matrix(object$data)) == 1)
      res$data.summary <- by2matrix(res$data.summary)
    if(lambda != 1){
      res$Tdata.summary <- by(Tdata, real,
                              function(x) drop(apply(as.matrix(x), 2, summary)))
      if(ncol(as.matrix(Tdata)) == 1)
        res$Tdata.summary <- by2matrix(res$Tdata.summary)
    }
    if(!is.null(object$units.m)){
      res$units.m.summary <- by(object$units.m, real,
                                function(x) drop(apply(as.matrix(x), 2, summary)))
      if(ncol(as.matrix(object$units.m)) == 1)
        res$units.m.summary <- by2matrix(res$units.m.summary)
    }
    if(!is.null(object$covariate)){
      res$covariate.summary <- by(object$covariate, real, summary)
      attr(res$covariate.summary, "dimnames") <- list(realisation = levels(real)) 
    }
    res$duplicated.coords <- lapply(sort(unique(object$realisations)), function(x) dup.coords(eval(substitute(subset(object, realisations == z), list(z=x)))))
    if(all(sapply(res$duplicated.coords, is.null))) res$duplicated.coords <- NULL
  }
  else{
    res$n <- nrow(object$coords)
    res$coords.summary <- apply(object$coords, 2, range)
    rownames(res$coords.summary) <- c("min", "max")
    if(is.null(colnames(object$coords)))
      colnames(res$coords.summary) <- c("Coord.X", "Coord.Y")
    res$distances.summary <- range(dist(object$coords))
    names(res$distances.summary) <- c("min", "max")  
    res$data.summary <- drop(apply(as.matrix(object$data), 2, summary))
    if(lambda != 1)
      res$Tdata.summary <- drop(apply(as.matrix(Tdata), 2, summary))
    if(!is.null(object$units.m))
      res$units.m.summary <- drop(apply(as.matrix(object$units.m), 2, summary))
    if(!is.null(object$covariate))
      res$covariate.summary <- summary(object$covariate)
    res$duplicated.coords <- dup.coords(object)
  }
  if(!is.null(object$borders)){
    res$borders.summary <- apply(object$borders, 2, range)
    rownames(res$borders.summary) <- c("min", "max")
  }
  others.ind <- is.na(match(names(unclass(object)),
                            c("coords","data","covariate","borders",
                              "realisations", "units.m")))
  if(sum(others.ind) > 0)
    res$others <- names(unclass(object[others.ind]))
  oldClass(res) <- "summary.geodata"
  return(res)
}

"print.summary.geodata" <- function(x, ...)
{
  if(length(x$n) == 1)
    cat(paste("Number of data points:",x$n,"\n"))
  else
    cat( paste(names(x$n),": number of data points:",x$n,"\n"))
  cat("\nCoordinates summary\n")
  print(x$coords.summary)
  cat("\nDistance summary\n")
  print(x$distances.summary)
  if(!is.null(x$borders.summary)){
    cat("\nBorders summary\n")
    print(x$borders.summary)
  }
  cat("\nData summary\n")
  print(x$data.summary)
  if(!is.null(x$Tdata.summary)){
    cat("\nTransformed Data summary\n")
    print(x$Tdata.summary)
  }
  if(!is.null(x$units.m.summary)){
    cat("\nOffset variable summary\n")
    print(x$units.m.summary)
  }
  if(!is.null(x$covariate.summary)){
    cat("\nCovariates summary\n")
    print(x$covariate.summary)
  }
  if(!is.null(x$duplicated.coords)){
    cat("\nDuplicated Coordinates\n") 
    print(x$duplicated.coords)
  }
  if(!is.null(x$others)){
    cat("\nOther elements in the geodata object\n") 
    print(x$others)
  }
  return(invisible())
}

"points.geodata" <-
  function(x, coords = x$coords, data = x$data, 
           data.col = 1, borders = x$borders,
           pt.divide = c("data.proportional",
             "rank.proportional", "quintiles",
             "quartiles", "deciles", "equal"),
           lambda=1, trend="cte", abs.residuals=FALSE,
           weights.divide = "units.m",
           cex.min, cex.max, cex.var,
           pch.seq, col.seq, add.to.plot = FALSE,
           x.leg, y.leg = NULL, dig.leg = 2, 
           round.quantiles = FALSE, permute = FALSE, ...) 
{
  ldots <- list(...)
  TYPE <- ldots$type
  if(length(ldots) > 0){
    ind <- pmatch(names(ldots), names(formals(plot.default)))
    names(ldots) <- ifelse(is.na(ind), names(ldots), names(formals(plot.default))[ind])
#    names(ldots) <- match.arg(names(ldots), names(formals(plot.default)), several.ok=TRUE)
  }
  ##
  ## Checking input
  ##
  if(missing(x)) x <- list(coords = coords, data = data)
  ## This is for compatibility with previously used argument pt.sizes
  if(!is.null(list(...)$pt.s)) pt.divide <- list(...)$pt.s
  ##
  if(mode(pt.divide) != "numeric") pt.divide <- match.arg(pt.divide)
  if(!is.vector(data)) data <- (as.data.frame(data))[,data.col]
  npts <- nrow(coords)
  if(nrow(coords) != length(data))
    stop("coords and data have incompatible sizes")
  if (!is.null(weights.divide)) {
    if(all(weights.divide == "units.m")){
      if(!is.null(x$units.m) && identical(data, x$data)){
        if(!all(x$units.m > 0))
          stop("all values in units.m must be greater than zero")
        data <- data/x$units.m
      }
    }
    else{
      if((length(weights.divide) != 1) &&
         (length(weights.divide) != length(data))) 
        stop("length of weights.divide must be equal to the length of data")
      data <- data/weights.divide
    }
  }
  if (missing(cex.min)) cex.min <- 0.5
  if (missing(cex.max)) cex.max <- 1.5
  ##
  ## data transformation (Box-Cox)
  ##
  if (lambda != 1) data <- BCtransform(x=data, lambda=lambda)$data
  ##
  ## trend "removal" (getting residuals from a linear model)
  ##
  xmat <- unclass(trend.spatial(trend = trend, geodata = x))
  if (nrow(xmat) != nrow(coords)) 
    stop("coords and trend have incompatible dimensions")
  if (trend != "cte") {
    data <- lm(data ~ xmat + 0)$residuals
    if(abs.residuals) abs.res <- abs(data)
    names(data) <- NULL
  }
  ##
  ## permuting data locations
  ##
  if(permute) data <- sample(data)
  ##
  ## symbols size proportional to the data or to an "external" variable
  ##
  if(missing(cex.var)){
    if(trend != "cte" && abs.residuals) cex.var.data <- abs.res
    else cex.var.data <- data
  }
  else{
    cex.var.data <- cex.var
    if(length(cex.var) != npts)
      stop("length of cex.var must be the same as the number of data locations")
  }
  ind <- order(cex.var.data)
  ind.order <- order(ind)
  r.y <- range(cex.var.data)
  size <- cex.min + ((cex.var.data[ind] - r.y[1]) *
                     (cex.max - cex.min))/(r.y[2] - r.y[1])
  ##
  attach(x, pos=2, warn.conflicts=FALSE)
  eval(borders, envir=as.environment(2))
  detach(2)
  if (!add.to.plot) {
    coords.lims <- set.coords.lims(coords=coords, borders = borders, ...)
    par(pty = "s")
    toplot <- apply(coords, 2, range)
    colnames(toplot) <- c("X Coord", "Y Coord")
    ldots$xlim <- coords.lims[,1]
    ldots$ylim <- coords.lims[,2]
    ldots$type= "n"
    do.call("plot", c(list(x=toplot), ldots))
  }
  if(!is.null(TYPE) && ldots$type == "n") return(invisible())
  if(!is.null(borders)) polygon(borders)
  graph.list <- list()
  ##
  ##
  if(mode(pt.divide) == "numeric" ||
     any(pt.divide %in% c("quintiles", "quartiles", "deciles"))) {
    if (all(pt.divide == "quintiles")) {
      n.quant <- 5
      if (missing(col.seq)) 
        col.seq <- c("blue", "green", "yellow", "orange3", "red2")
    }
    if (all(pt.divide == "quartiles")) {
      n.quant <- 4
      if (missing(col.seq)) 
        col.seq <- c("blue", "green", "yellow", "red") 
    }
    if (all(pt.divide == "deciles")) {
      n.quant <- 10
      if (missing(col.seq)) 
        col.seq <- terrain.colors(46)[seq(1,46,by=5)]
    }
    if(mode(pt.divide) == "numeric"){
      if(length(pt.divide) == 1){
        n.quant <- pt.divide
        data.quantile <- quantile(data, probs = seq(0, 1, by = (1/n.quant)))
        if (missing(col.seq))
          col.seq <- "gray"
      }
      else{
        if(length(pt.divide) <= length(data)){
          data.quantile <- pt.divide
          n.quant <- length(pt.divide) - 1
        }
        else
          stop("length of pt.divide cannot be greater than length of the data")
      }
    }
    else
      data.quantile <- quantile(data, probs = seq(0, 1, by = (1/n.quant)))
    if(!missing(col.seq) && all(col.seq == "gray"))
      col.seq <- gray(seq(1, 0, l=n.quant))
    if(length(col.seq) > n.quant)
      col.seq <- col.seq[round(seq(1, length(col.seq), length=n.quant))]
    if(missing(pch.seq)) pch.seq <- rep(21, n.quant)
    if(round.quantiles) {
      data.quantile[1] <- floor(data.quantile[1])
      data.quantile[n.quant + 1] <- ceiling(data.quantile[n.quant + 1])
      data.quantile <- round(data.quantile)
    }
    graph.list$quantiles <- data.quantile
    graph.list$data.group <- cut(data, breaks=data.quantile, include.l=TRUE)
    if(missing(cex.var))
      size <- seq(cex.min, cex.max, l = n.quant)[as.numeric(graph.list$data.group)]
    else size <- size[ind.order]
    graph.list$cex <- size
    graph.list$pch <- pch.seq
    graph.list$col <- col.seq[as.numeric(graph.list$data.group)]
    if (add.to.plot) 
      points(coords, pch = pch.seq, cex = size, bg = graph.list$col, ...)
    else
      points(coords, pch = pch.seq, cex = size, bg = graph.list$col)
    ##
    ## Adding legend
    ##
    if(!missing(x.leg)){
      nums <- formatC(graph.list$quantiles, digits=dig.leg, format="f")
      textleg <- expression()
      for (i in 1:(length(graph.list$quantiles)-1))
        textleg <- c(textleg,
                     substitute(group("[", list(a, b), ")"),
                                list(a=nums[i], b=nums[i+1])))
      legend(x=x.leg, y=y.leg, legend=textleg,
             pt.bg=col.seq, col=col.seq, pch=graph.list$pch)
    }
  }
  else {
    n <- length(data)
    if (missing(pch.seq)) pch.seq <- 21
    ##    coords.order <- coords[ind, ]
    ind.d <- order(data)
    ind.d.order <- order(ind.d)
    if (pt.divide == "rank.proportional") {
      if(missing(cex.var))
        size <- seq(cex.min, cex.max, l = n)
    }
    if (pt.divide == "data.proportional") {
      if(missing(cex.var)){
        r.y <- range(data)
        size <- cex.min + ((sort(data) - r.y[1]) *
                           (cex.max - cex.min))/(r.y[2] - r.y[1])
      }
    }
    if (pt.divide == "equal") size <- cex.max
    else size <- size[ind.order]
    if (missing(col.seq)) col.seq <- 0
    if(all(col.seq == "gray")) col.seq <- gray(seq(1,0.1, l=n))
    if (length(col.seq) == 1) col.seq <- rep(col.seq, n)
    if (length(col.seq) != n)
      col.seq <- col.seq[round(seq(1,length(col.seq),length=n))]
    col.seq <- col.seq[ind.d.order]
    graph.list$cex <- size
    if(mode(pch.seq) == "numeric")
      graph.list$pch <- unique(range(pch.seq))
    else
      graph.list$pch <- pch.seq
    graph.list$col <- col.seq
    ##
    if (add.to.plot) 
      points(coords, cex = size, pch = pch.seq, bg = col.seq, ...)
    else points(coords, cex = size, pch = pch.seq, bg = col.seq)
    if(!missing(x.leg) && !missing(y.leg))
      warning(paste('arguments x.leg and y.leg are ignored when pt.divide = ', pt.divide,'\n'))
  }
  return(invisible(graph.list))
}

"plot.geodata" <-
  function (x, coords = x$coords, data = x$data, borders = x$borders, 
    trend = "cte", lambda = 1, col.data = 1, weights.divide = "units.m", 
    lowess = FALSE, scatter3d = FALSE, density = TRUE, rug = TRUE, qt.col, ...) 
{
  if(missing(x)) x <- list()
  x$coords <- coords
  x$data <- data
  if(missing(qt.col)) qt.col <- c("blue", "green", "yellow2", "red")
  if(length(qt.col) == 1) qt.col <- rep(qt.col, 4)
  par.ori <- par(no.readonly = TRUE)
  on.exit(par(par.ori))
  coords <- as.matrix(coords)
  data <- as.matrix(data)
  data <- data[, col.data]
  if(missing(borders)) borders <- x$borders
  attach(x, pos=2, warn.conflicts=FALSE)
  eval(borders, envir=as.environment(2))
  detach(2)
  if (!is.null(weights.divide)) {
    if(all(weights.divide == "units.m")){
      if(!is.null(x$units.m) && identical(data, x$data)){
        if(!all(x$units.m > 0))
          stop("all values in units.m must be greater than zero")
        data <- data/x$units.m
      }
    }
    else{
      if((length(weights.divide) != 1) &&
         (length(weights.divide) != length(data))) 
        stop("length of weights.divide must be equals to the length of data")
      data <- data/weights.divide
    }
  }
  if (lambda != 1) {
    if (lambda == 0) 
      data <- log(data)
    else data <- ((data^lambda) - 1)/lambda
  }
  xmat <- unclass(trend.spatial(trend = trend, geodata = x))
  if (nrow(xmat) != nrow(coords)) 
    stop("coords and trend have incompatible sizes")
  if (trend != "cte") {
    data <- lm(data ~ xmat + 0)$residuals
    names(data) <- NULL
    data.lab <- "residuals"
  }
  else data.lab <- "data"
  par(mfrow = c(2, 2), mar = c(4, 4, 0, 0.5), mgp=c(2,.8,0))
  if (is.null(borders)) 
    coords.lims <- set.coords.lims(coords = coords, ...)
  else {
    borders <- as.matrix(as.data.frame(borders))
    if (ncol(borders) != 2) 
      stop("argument \"borders\" must be a 2 columns object with coordinates of the borders of the study area")
    coords.lims <- set.coords.lims(coords = rbind(as.matrix(coords), as.matrix(borders)), ...)
  }
  par(pty = "s")
  plot(coords, xlab = "X Coord", ylab = "Y Coord ", type = "n", 
       xlim = coords.lims[, 1], ylim = coords.lims[, 2])
  if (!is.null(borders)) polygon(borders)
  data.breaks <- unique(quantile(data))
  data.cut <- cut(data, breaks = data.breaks, include.l = TRUE, labels = FALSE)
  points(coords, pch = (1:4)[data.cut], col = qt.col[data.cut])
  plot(data, coords[, 2], ylab = "Y Coord", xlab = data.lab, cex = 1, ylim = coords.lims[, 2])
  if(lowess){
    foo <- lowess(data ~ coords[,2])
    lines(foo[[2]], foo[[1]])
  }
  par(mar = c(5, 5, 1, 0.5))
  plot(coords[, 1], data, xlab = "X Coord", ylab = data.lab, cex = 1, xlim = coords.lims[, 1], )
  if(lowess) lines(lowess(data ~ coords[,1]))
  par(pty = "m")
  par(mar = c(4, 4, 1, 1))
  if (scatter3d && !requireNamespace("scatterplot3d", quietly=TRUE)){
    scatter3d <- FALSE
    cat("plot.geodata: the argument scatter3d=TRUE requires the package scatterplot3d\n which is not available, will plot an histogram instead")
  }
  if(scatter3d)
    scatterplot3d::scatterplot3d(x = coords[, 1],
                                  y = coords[, 2], z = data,
                                  box = FALSE, type = "h", pch = 16,
                                  xlab= "X Coord", ylab = "Y Coord", ...)
  else{
    hist(data, main="", xlab= data.lab, prob=TRUE, ...)
    if(density) lines(density(data))
    if(rug) rug(data)
  }
  ##  else xyzplot(coords = coords, data = data, ...)
  return(invisible())
}

"sample.geodata" <-
  function(x, size, replace = FALSE, prob = NULL,
           coef.logCox, external)
{
  if(all(class(x) != "geodata") & all(class(x) != "grf"))
    stop("object must be of the class \"geodata\" or \"grf\"")
  if(all(class(x) == "grf"))
    x <- as.geodata(x)
  if(is.null(prob)){
    if(!missing(coef.logCox)){
      if(!is.numeric(coef.logCox) || length(coef.logCox) > 1)
        stop("coef.logCox must be an scalar")
      if(missing(external))
        prob <- exp(coef.logCox*x$data)
      else{
        if(!is.numeric(coef.logCox) && !is.vector(external))
          stop("external must be a numeric vector")
        if(length(external) != nrow(x$coords))
          stop("length of external must be of the same as the number
 of coordinates")
        else
          prob <- exp(coef.logCox*external)
      }
    }
  }
  ind <- sample(1:nrow(x$coords), size = size, prob=prob)
  as.geodata(as.data.frame(x)[ind,])
}




