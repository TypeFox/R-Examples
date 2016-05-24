lnre.gigp <- function (gamma=-.5, B=.01, C=.01, param=list())
{
  if (! is.list(param)) stop("'param' argument must be a list of parameter=value pairs")
  pnames <- names(param)

  ## explicitly specified parameters override the values in 'param'
  if (!missing(gamma) || !("gamma" %in% pnames)) param$gamma <- gamma
  if (!missing(B) || !("B" %in% pnames)) param$B <- B
  if (!missing(C) || !("C" %in% pnames)) param$C <- C

  ## initialize lnre.fzm model object
  self <- list(type="gigp", name="Generalized Inverse Gauss-Poisson (GIGP)",
               param=list(), param2=list(),
               util=list(update=lnre.gigp.update, transform=lnre.gigp.transform, print=lnre.gigp.print))
  class(self) <- c("lnre.gigp", "lnre", class(self))

  ## update model parameters to specified values & compute secondary parameters
  self <- lnre.gigp.update(self, param)

  self
}


lnre.gigp.update <- function (self, param=list(), transformed=FALSE)
{
  if (! is.list(param)) stop("'param' argument must be a list")
  if (! inherits(self, "lnre.gigp")) stop("this function can only be applied to 'lnre.gigp' objects")

  if (transformed) param <- lnre.gigp.transform(param, inverse=TRUE)
  
  pnames <- names(param)
  unused <- !(pnames %in% c("gamma", "B", "C"))
  if (any(unused)) warning("unused parameters in 'param' argument: ", pnames[unused])

  if ("gamma" %in% pnames) {
    gamma <- param$gamma
    if (gamma <= -1 || gamma >= 0) stop("parameter gamma = ",gamma," out of range (-1,0)")
    self$param$gamma <- gamma
  }
  else {
    gamma <- self$param$gamma
  }
    
  if ("B" %in% pnames) {
    B <- param$B
    if (B < 0) stop("parameter B = ",B," out of range [0, Inf)")
    self$param$B <- B
  }
  else {
    B <- self$param$B
  }
  
  if ("C" %in% pnames) {
    C <- param$C
    if (C < 0) stop("parameter C = ",C," out of range [0, Inf)")
    self$param$C <- C
  }
  else {
    C <- self$param$C
  }

  Z <- 1 / C
  S <- (2 * besselK(B, gamma)) / (besselK(B, gamma+1) * B * C)
  self$param2$Z <- Z
  self$S <- S

  self
}

lnre.gigp.transform <- function (param, inverse=FALSE)
{
  gamma <- param$gamma
  B <- param$B
  C <- param$C
  new.param <- list()

  if (inverse) {
    ## gamma = -inv.logit(gamma*) = -1 / (1 + exp(-gamma*))
    if (! is.null(gamma)) new.param$gamma <- -1 / (1 + exp(-gamma))
    ## B = exp(B* - 5)
    if (! is.null(B)) new.param$B <- exp(B - 5)
    ## C = exp(C* - 5)
    if (! is.null(C)) new.param$C <- exp(C - 5)
  }
  else {
    ## gamma* = logit(-gamma) = log(-gamma / (1+gamma))
    if (! is.null(gamma)) new.param$gamma <- log(-gamma / (1+gamma))
    ## B* = log(B) + 5 [shifted so that B* = 0 -> B = .0067 is a useful init value]
    if (! is.null(B)) new.param$B <- log(B) + 5
    ## C* = log(C) + 5 [shifted so that C* = 0 -> C = .0067 is a useful init value]
    if (! is.null(C)) new.param$C <- log(C) + 5
  }

  new.param
}

lnre.gigp.print <- function (self)
{
  cat("Generalized Inverse Gauss-Poisson (GIGP) LNRE model.\n")
  cat("Parameters:\n")
  cat("   Shape:          gamma =", self$param$gamma, "\n")
  cat("   Lower decay:        B =", self$param$B, "\n")
  cat("   Upper decay:        C =", self$param$C, "\n")
  cat(" [ Zipf size:          Z =", self$param2$Z, "]\n")    
}
