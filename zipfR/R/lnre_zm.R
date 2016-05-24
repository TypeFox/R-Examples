lnre.zm <- function (alpha=.8, B=.01, param=list())
{
  if (! is.list(param)) stop("'param' argument must be a list of parameter=value pairs")
  pnames <- names(param)

  ## explicitly specified parameters override the values in 'param'
  if (!missing(alpha) || !("alpha" %in% pnames)) param$alpha <- alpha
  if (!missing(B) || !("B" %in% pnames)) param$B <- B

  ## initialize lnre.zm model object
  self <- list(type="zm", name="Zipf-Mandelbrot",
               param=list(), param2=list(),
               util=list(update=lnre.zm.update, transform=lnre.zm.transform, print=lnre.zm.print))
  class(self) <- c("lnre.zm", "lnre", class(self))

  ## update model parameters to specified values & compute secondary parameters
  self <- lnre.zm.update(self, param)

  self
}


lnre.zm.update <- function (self, param=list(), transformed=FALSE)
{
  if (! is.list(param)) stop("'param' argument must be a list")
  if (! inherits(self, "lnre.zm")) stop("this function can only be applied to 'lnre.zm' objects")

  if (transformed) param <- lnre.zm.transform(param, inverse=TRUE)
  
  pnames <- names(param)
  unused <- !(pnames %in% c("alpha", "B"))
  if (any(unused)) warning("unused parameters in 'param' argument: ", pnames[unused])

  if ("alpha" %in% pnames) {
    alpha <- param$alpha
    if (alpha <= 0 || alpha >= 1) stop("parameter alpha = ", alpha, " out of range (0,1)")
    self$param$alpha <- alpha
  }
  else {
    alpha <- self$param$alpha
  }
    
  if ("B" %in% pnames) {
    B <- param$B
    if (B <= 0) stop("parameter B = ", B, " out of range (0, Inf)")
    self$param$B <- param$B
  }
  else {
    B <- self$param$B
  }

  C <- (1 - alpha) / B^(1 - alpha)      # Evert (2004), p. 126
  self$param2$C <- C
  self$S <- Inf

  self
}

lnre.zm.transform <- function (param, inverse=FALSE)
{
  alpha <- param$alpha
  B <- param$B
  new.param <- list()

  if (inverse) {
    ## alpha = inv.logit(alpha*) = 1 / (1 + exp(-alpha*))
    if (! is.null(alpha)) new.param$alpha <- 1 / (1 + exp(-alpha))
    ## B = exp(B*)
    if (! is.null(B)) new.param$B <- exp(B)
  }
  else {
    ## alpha* = logit(alpha) = log(alpha / (1-alpha))
    if (! is.null(alpha)) new.param$alpha <- log(alpha / (1-alpha))
    ## B* = log(B)
    if (! is.null(B)) new.param$B <- log(B)
  }

  new.param
}

lnre.zm.print <- function (self)
{
  cat("Zipf-Mandelbrot LNRE model.\n")
  cat("Parameters:\n")
  cat("   Shape:          alpha =", self$param$alpha, "\n")
  cat("   Upper cutoff:       B =", self$param$B, "\n")
  cat(" [ Normalization:      C =", self$param2$C, "]\n")    
}
