
setGeneric(name = "transPars", def = function(x, 
    type = c("square", "StructTS", "exp", "exp10sq"), gradient = FALSE, hessian = FALSE, 
    rp, sclrho = 1.7, sclomega = 1.7, ftrans = NULL, ...) {
      standardGeneric("transPars") })

setMethod(f = "transPars", signature = "numeric", 
  definition = function(x, 
    type = eval(formals(stsm.class::transPars)$type), gradient = FALSE, hessian = FALSE,
    rp, sclrho = 1.7, sclomega = 1.7, ftrans = NULL, ...) {

  tpStructTS <- function(pars, rp)
  {
    if (length(idvar) > 0) {
      tpars[idvar] <- pars[idvar] * rp
    } #else warning("No changes done by 'transPars'.")
    # changes may be done in other parameters, e.g. 'phi1', 'phi2'

    if (gradient)
    {
      if (length(idvar) > 0)
        d1[idvar] <- rp
    }
    #if (hessian)
    #  d2[] <- 0 #array(0, dim = c(p, p))
    list(pars = tpars, d1 = d1, d2 = d2)
  }

  tpsq <- function(pars)
  {
    if (length(idvar) > 0) {
      tpars[idvar] <- pars[idvar]^2
    } #else warning("No changes done by 'transPars'.")

    if (gradient)
    {
      if (length(idvar) > 0)
        d1[idvar] <- 2 * pars[idvar]
    }
    if (hessian) {
      #d2[] <- 0 #array(0, dim = c(p, p)) #diag(2, p, p)
      diag(d2)[idvar] <- 2
    }
    list(pars = tpars, d1 = d1, d2 = d2)
  }

  tpexp <- function(pars)
  {
    if (length(idvar) > 0) {
      tpars[idvar] <- exp(pars[idvar])
    } #else warning("No changes done by 'transPars'.")

    if (gradient)
    {
      if (length(idvar) > 0)
        d1[idvar] <- tpars[idvar]
    }
    if (hessian) {
      #d2[] <- 0 #array(0, dim = c(p, p)) #diag(2, p, p)
      diag(d2)[idvar] <- tpars[idvar]
    }
    list(pars = tpars, d1 = d1, d2 = d2)
  }

  tpexp10sq <- function(pars)
  {
    if (length(idvar) > 0) {
      tpars[idvar] <- (exp(-pars[idvar])/10)^2
    } #else warning("No changes done by 'transPars'.")

##FIXME TODO
    if (gradient)
    {
      d1 <- NULL
    }
    if (hessian) {
      d2 <- NULL
    }
    list(pars = tpars, d1 = d1, d2 = d2)
  }

  if (is.function(ftrans))
    return(ftrans(x = x, gradient = gradient, hessian = hessian, ...))

  type <- match.arg(type)[1]
  if (type == "StructTS")
    stopifnot(!is.null(rp))
  tpars <- pars <- x

  label <- "var"
  p <- length(pars)
  nmspars <- names(pars)
  idvar <- grep("^var|P0\\d{1,2}$", nmspars, value = FALSE)
  
  if (gradient) {
    d1 <- rep(NA, p)
    names(d1) <- nmspars
  } else d1 <- NULL
  if (hessian) {
    d2 <- matrix(0, p, p)
    rownames(d2) <- colnames(d2) <- nmspars
  } else d2 <- NULL

  res <- switch(type, "StructTS" = tpStructTS(pars, rp), "square" = tpsq(pars), 
    "exp" = tpexp(pars), "exp10sq" = tpexp10sq(pars)) #"4" = tpStructTS(pars, rp),

  idphi <- grep("^phi\\d", nmspars)

  if (length(idphi) > 0)
  {
    stopifnot(length(idphi) == 2)

    opabsx <- 1 + abs(pars[idphi])
    aux <- pars[idphi] / opabsx
    res$pars[idphi] <- c(sum(aux), -prod(aux))

    if (gradient)
    {
      res$d1["phi1"] <- (opabsx["phi1"] - 
        pars["phi1"] * sign(pars["phi1"])) / opabsx["phi1"]^2
      res$d1["phi2"] <- -aux["phi1"] * (opabsx["phi2"] - 
        pars["phi2"] * sign(pars["phi2"])) / opabsx["phi2"]^2
    }

    if (hessian)
    {
##FIXME TOO
    }
  }

  list(pars = res$pars, gradient = res$d1, hessian = res$d2)
})

setMethod(f = "transPars", signature = "stsm", 
  definition = function(x, type = NULL,
    gradient = FALSE, hessian = FALSE,     
    rp, sclrho = 1.7, sclomega = 1.7, ftrans = NULL, ...) {

  if (is.null(x@transPars))
    return(list(pars = x@pars))

  if (is.function(x@transPars)) {
    type <- NULL
    ftrans <- x@transPars
    rp <- NULL
  } else {
    #type <- match.arg(x@transPars, c("StructTS", "square", "exp", "4", "nk"))
    type <- match.arg(x@transPars, eval(formals(stsm.class::transPars)$type))

    rp <- switch(type, 
      "StructTS" = var(x@y, na.rm = TRUE) / 100,
      NULL) # default
  }
    
  transPars(x@pars, type = type, gradient = gradient, hessian = hessian,
    rp = rp, sclrho = sclrho, sclomega = sclomega, ftrans = ftrans, ...)
})
