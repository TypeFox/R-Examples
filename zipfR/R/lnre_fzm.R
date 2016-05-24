lnre.fzm <- function (alpha=.8, A=1e-9, B=.01, param=list())
{
  if (! is.list(param)) stop("'param' argument must be a list of parameter=value pairs")
  pnames <- names(param)

  ## explicitly specified parameters override the values in 'param'
  if (!missing(alpha) || !("alpha" %in% pnames)) param$alpha <- alpha
  if (!missing(A) || !("A" %in% pnames)) param$A <- A
  if (!missing(B) || !("B" %in% pnames)) param$B <- B

  ## initialize lnre.fzm model object
  self <- list(type="fzm", name="finite Zipf-Mandelbrot",
               param=list(), param2=list(),
               util=list(update=lnre.fzm.update, transform=lnre.fzm.transform, print=lnre.fzm.print))
  class(self) <- c("lnre.fzm", "lnre", class(self))

  ## update model parameters to specified values & compute secondary parameters
  self <- lnre.fzm.update(self, param)

  self
}


lnre.fzm.update <- function (self, param=list(), transformed=FALSE)
{
  if (! is.list(param)) stop("'param' argument must be a list")
  if (! inherits(self, "lnre.fzm")) stop("this function can only be applied to 'lnre.fzm' objects")

  if (transformed) param <- lnre.fzm.transform(param, inverse=TRUE)
  
  pnames <- names(param)
  unused <- !(pnames %in% c("alpha", "A", "B"))
  if (any(unused)) warning("unused parameters in 'param' argument: ", pnames[unused])

  if ("alpha" %in% pnames) {
    alpha <- param$alpha
    if (alpha <= 0 || alpha >= 1) stop("parameter alpha = ",alpha," out of range (0,1)")
    self$param$alpha <- alpha
  }
  else {
    alpha <- self$param$alpha
  }
    
  if ("B" %in% pnames) {
    B <- param$B
    if (B <= 0) stop("parameter B = ",B," out of range (0, Inf)")
    self$param$B <- B
  }
  else {
    B <- self$param$B
  }
  
  if ("A" %in% pnames) {
    A <- param$A
    if (A <= 0 || A >= B) stop("parameter A = ",A," out of range (0, B) = (0, ",B,")")
    self$param$A <- A
  }
  else {
    A <- self$param$A
  }

  
  C <- (1 - alpha) / ( B ^ (1 - alpha) - A ^ (1 - alpha) ) # Evert (2004), p. ~~ TODO ~~
  S <- (C / alpha) * ((A ^ -alpha) - (B ^ -alpha))
  self$param2$C <- C
  self$S <- S

  self
}

lnre.fzm.transform <- function (param, inverse=FALSE)
{
  alpha <- param$alpha
  A <- param$A
  B <- param$B
  new.param <- list()

  if (inverse) {
    ## alpha = inv.logit(alpha*) = 1 / (1 + exp(-alpha*))
    if (! is.null(alpha)) new.param$alpha <- 1 / (1 + exp(-alpha))
    ## A = exp(A* + log(1e-9))
    if (! is.null(A)) new.param$A <- exp(A + log(1e-9))
    ## B = exp(B*)
    if (! is.null(B)) new.param$B <- exp(B)
  }
  else {
    ## alpha* = logit(alpha) = log(alpha / (1-alpha))
    if (! is.null(alpha)) new.param$alpha <- log(alpha / (1-alpha))
    ## A* = log(A) - log(1e-9) [shifted so that A* = 0 is a sensible init value]
    if (! is.null(A)) new.param$A <- log(A) - log(1e-9)
    ## B* = log(B)
    if (! is.null(B)) new.param$B <- log(B)
  }

  new.param
}

lnre.fzm.print <- function (self)
{
  cat("finite Zipf-Mandelbrot LNRE model.\n")
  cat("Parameters:\n")
  cat("   Shape:          alpha =", self$param$alpha, "\n")
  cat("   Lower cutoff:       A =", self$param$A, "\n")
  cat("   Upper cutoff:       B =", self$param$B, "\n")
  cat(" [ Normalization:      C =", self$param2$C, "]\n")    
}
