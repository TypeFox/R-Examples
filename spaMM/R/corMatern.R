#require(nlme)
#require(HL)
#dyn.load("corMatern.so")

corMatern <- function(value = c(1,0.5), form = ~1, nugget=FALSE, nuScaled=FALSE, metric = c("euclidean", "maximum", "manhattan"), fixed = FALSE)
{
    attr(value, "formula") <- form
    attr(value, "nugget") <- nugget
    attr(value, "nuScaled") <- nuScaled
    attr(value, "metric") <- match.arg(metric)
    attr(value, "fixed") <- fixed
    class(value) <- c("corMatern", "corStruct") ## corSpatial removed
    value
}

getCovariate.corMatern <-
  function(object, form = formula(object), data)
{
  if (is.null(covar <- attr(object, "covariate"))) { # need to calculate it
    if (missing(data)) {
      stop("need data to calculate covariate")
    }
    covForm <- getCovariateFormula(form)
    if (length(all.vars(covForm)) > 0) { # covariate present
      if (attr(terms(covForm), "intercept") == 1) {
	covForm <-
          eval(parse(text = paste("~", DEPARSE(covForm[[2]]),"-1",sep="")))
      }
      covar <-
          as.data.frame(unclass(model.matrix(covForm,
                                             model.frame(covForm, data,
                                                         drop.unused.levels = TRUE))))
    } else {
      covar <- NULL
    }

    if (!is.null(getGroupsFormula(form))) { # by groups
      grps <- getGroups(object, data = data)
      if (is.null(covar)) {
	covar <- lapply(split(grps, grps),
                        function(x) as.vector(proxy::dist(1:length(x))))
      } else {
	covar <- lapply(split(covar, grps),
			function(el, metric) {
                          el <- as.matrix(el)
                          if (nrow(el) > 1) {
                            as.vector(proxy::dist(el, metric))
                          } else {
                            numeric(0)
                          }
			}, metric = attr(object, "metric"))
      }
      covar <- covar[sapply(covar, length) > 0]  # no 1-obs groups
    } else {				# no groups
      if (is.null(covar)) {
	covar <- as.vector(proxy::dist(1:nrow(data)))
      } else {
	covar <- as.vector(proxy::dist(as.matrix(covar),
                                method = attr(object, "metric")))
      }
    }
    if (any(unlist(covar) == 0)) {
      stop("cannot have zero distances in \"corMatern\"")
    }
  }
  covar
}

Initialize.corMatern <- function(object, data, ...)
{
    if (!is.null(attr(object, "minD"))) { #already initialized
      return(object)
    }
    
    form <- formula(object)
    ## obtaining the groups information, if any
    if (!is.null(getGroupsFormula(form))) {
      attr(object, "groups") <- getGroups(object, form, data = data)
      attr(object, "Dim") <- Dim(object, attr(object, "groups"))
    } else {
      # no groups
      attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
    }
    ## obtaining the covariate(s)
    ## is this where the distance matrix is actually computed ?
    attr(object, "covariate") <- getCovariate(object, data = data)

    nug  <- attr(object, "nugget")
    #nusc <- attr(object, "nuScaled")
    
    val <- as.vector(object)            # how many parameters?
    if (length(val) > 0) {		# is initialized
    if (val[1] <= 0) {                  # test for values of range
      stop("'range' must be > 0 in \"corMatern\" initial value")
    }
    if (val[2] <= 0) {                  # test for values of nu
      stop("'nu' must be > 0 in \"corMatern\" initial value")
    }    
    if (nug) {				# with nugget effect
      if (length(val) == 2) {		# assuming nugget effect not given
	    val <- c(val, 0.1)		# setting it to 0.1
      } else {
	if (length(val) != 3) {
	  stop("initial value for \"corMatern\" parameters of wrong dimension")
	}
      }
      if ((val[3] <= 0) || (val[3] >= 1)) {
	stop("initial value of nugget ratio must be between 0 and 1")
      }
    } else {				# only range and nu
      if (length(val) != 2) {
	stop("initial value for \"corMatern\" parameters of wrong dimension")
      }
    }
  } else { ## val of length 0
    val <- min(unlist(attr(object, "covariate"))) * 0.9
    val <- c(val,0.5) ## FR 18/03/13
    if (nug) val <- c(val, 0.1)
  }
  val[1] <- log(val[1]) ## range
  val[2] <- log(val[2]) ## nu    
  if (nug) val[3] <- log(val[3]/(1 - val[3]))
  oldAttr <- attributes(object)
  object <- val
  attributes(object) <- oldAttr
  attr(object, "minD") <- min(unlist(attr(object, "covariate")))
  attr(object, "factor") <- corFactor(object)
  attr(object, "logDet") <- -attr(attr(object, "factor"), "logDet")

  object
}

logDet.corMatern <-
  function(object, covariate = getCovariate(object), ...)
{
  if (!is.null(aux <- attr(object, "logDet"))) {
    return(aux)
  }
  if (is.null(aux <- attr(object, "factor"))) {
    ## getting the transpose sqrt factor
    aux <- corMatrix(object, covariate = covariate, corr = FALSE)
  }
  if (is.null(aux1 <- attr(aux, "logDet"))) {
    ## checking for logDet attribute; if not present, get corr matrix
    aux <- corMatrix(object, covariate)
    if (data.class(aux) == "list") {    # by group
      sum(log(abs(unlist(lapply(aux, function(el) svd(el)$d)))))/2
    } else {
      sum(log(abs(svd(aux)$d)))/2
    }
  } else {
    -aux1
  }
}

corMatrix.corMatern <-
    function(object, covariate = getCovariate(object), 
             corr = TRUE, ## returns matrix, else returna cholesky factor + the logDet  
             ...)
{

  if (data.class(covariate) == "list") { ## groups are defined
    if (is.null(names(covariate))) {
      names(covariate) <- 1:length(covariate)
    }
    corD <- Dim(object, rep(names(covariate),
			    unlist(lapply(covariate,
                                          function(el) round((1 + sqrt(1 + 8 * length(el)))/2)))))
  } else { ## no groups
    corD <- Dim(object, rep(1, round((1 + sqrt(1 + 8* length(covariate)))/2)))
  }

  if (corr) {
    val <- .C("matern_matList",
	      as.double(as.vector(object)),
	      as.integer(attr(object, "nugget")),
              as.integer(attr(object, "nuScaled")),
	      as.double(unlist(covariate)), #distances ?
	      as.integer(unlist(corD)),
	      as.double(attr(object, "minD")),
	      mat = double(corD[["sumLenSq"]]))[["mat"]]
    lD <- NULL
  } else {
    val <- .C("matern_factList",
              as.double(as.vector(object)),
              as.integer(attr(object, "nugget")),
              as.integer(attr(object, "nuScaled")),
              as.double(unlist(getCovariate(object))),
              as.integer(unlist(corD)),
              as.double(attr(object, "minD")),
              factor = double(corD[["sumLenSq"]]),
              logDet = double(1))[c("factor", "logDet")]
    lD <- val[["logDet"]]
    val <- val[["factor"]]
  }
  if (corD[["M"]] > 1) {
    val <- split(val, rep(1:corD[["M"]], (corD[["len"]])^2))
    val <- lapply(val, function(el) {
      nel <- round(sqrt(length(el)))
      array(el, c(nel, nel))
    })
    names(val) <- names(corD[["len"]])
    val <- as.list(val)
  } else {
    val <- array(val, c(corD[["N"]], corD[["N"]]))
  }
  attr(val, "logDet") <- lD
  val
}

coef.corMatern <-
  function(object, unconstrained = TRUE, ...)
{
  ##cat("tata",object[1],object[2],"\n")
  if (attr(object, "fixed") && unconstrained) {
    return(numeric(0))
  }
  val <- as.vector(object)
  if (length(val) == 0) {               # uninitialized
    return(val)
  }
  if (!unconstrained) {
    val <- exp(val)
    if (attr(object, "nugget")) val[3] <- val[3]/(1+val[3])
  }
  if (attr(object, "nugget")) names(val) <- c("range", "nu", "nugget")
  else names(val) <- c("range","nu")
  val
}

"coef<-.corMatern" <-
  function(object, ..., value)
{
  if (length(value) != length(object)) {
    stop("cannot change the length of the parameter after initialization")
  }
#  cat("tutu",object[1],"->",value[1],object[2],"->",value[2],"\n")
  object[] <- value
  corD <- attr(object, "Dim")
  ## updating the factor list and logDet
  aux <- .C("matern_factList",
	    as.double(as.vector(object)),
	    as.integer(attr(object, "nugget")),
            as.integer(attr(object, "nuScaled")),
	    as.double(unlist(getCovariate(object))),
	    as.integer(unlist(corD)),
	    as.double(attr(object, "minD")),
	    factor = double(corD[["sumLenSq"]]),
	    logDet = double(1))[c("factor", "logDet")]
  attr(object, "factor") <- aux[["factor"]]
  attr(object, "logDet") <- -aux[["logDet"]]
  object
}


corFactor.corMatern <-
  function(object, ...)
{
  corD <- Dim(object)
  val <- .C("matern_factList",
	    as.double(as.vector(object)),
	    as.integer(attr(object, "nugget")),
            as.integer(attr(object, "nuScaled")),
	    as.double(unlist(getCovariate(object))),
	    as.integer(unlist(corD)),
	    as.double(attr(object, "minD")),
	    factor = double(corD[["sumLenSq"]]),
	    logDet = double(1))[c("factor", "logDet")]
  lD <- val[["logDet"]]
  val <- val[["factor"]]
  attr(val, "logDet") <- lD
  val
}

Dim.corMatern <-
  function(object, groups, ...)
{
  if (missing(groups)) return(attr(object, "Dim"))
  ugrp <- unique(groups)
  groups <- factor(groups, levels = ugrp)
  len <- table(groups)
  val <- list(N = length(groups),
              M = length(len),
              maxLen = max(len),
              sumLenSq = sum(len^2),
              len = len,
              start = NA,
              spClass=0)
  val[["start"]] <- c(0, 4, cumsum(val[["len"]] * (val[["len"]] - 1)/2)[-val[["M"]]])
  val
}

recalc.corMatern <-
  function(object, conLin, ...)
{
  ##print(sys.function(1))
  val <-
    .C("matern_recalc",
       Xy = as.double(conLin[["Xy"]]),
       as.integer(unlist(Dim(object))),
       as.integer(ncol(conLin[["Xy"]])),
       as.double(as.vector(object)),
       as.double(unlist(getCovariate(object))),
       as.double(attr(object, "minD")),
       as.integer(attr(object, "nugget")),
       as.integer(attr(object, "nuScaled")),
       logLik = double(1))[c("Xy", "logLik")]
  conLin[["Xy"]][] <- val[["Xy"]]
  conLin[["logLik"]] <- conLin[["logLik"]] + val[["logLik"]]
  conLin
}

Variogram.corMatern <-
  function(object, distance = NULL, sig2 = 1, length.out = 50, FUN, ...)
{
  if (is.null(distance)) {
    rangeDist <- range(unlist(getCovariate(object)))
    distance <- seq(rangeDist[1], rangeDist[2], length = length.out)
  }
  params <- coef(object, unconstrained = FALSE)
  if (length(params) == 2) {            # no nugget effect
    range  <- params[1]
    nu   <- params[2]
    nugg <- 0
  } else {                              # nugget effect
    range  <- params[1]
    nu   <- params[2]
    nugg <- params[3]
  }
  cat(range,nu,nugg,"\n")
  f <- function(distance) {
    z <- rep(0,length(distance))
    for(i in 1:length(distance)) {
      z[i] <- .C("matern_cor",
         c(range,nu,nugg),
         as.double(distance[i]),
         as.integer(attr(object, "nugget")),
         as.integer(attr(object, "nuScaled")),
         double(1))[[5]]
    }
    return(z)
  }
  ## static void matern(double *par, double dist, longint *nug, longint *nusc,
  ##          double *cor)

  zz <- f(distance)
  print(zz)
  val <- data.frame(variog = sig2 * (nugg + (1 - nugg) * zz),
                    dist = distance)
  class(val) <- c("Variogram", "data.frame")
  val
}
