
maxlik.fd.optim <- function(m, 
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999, 
  method = c("BFGS", "L-BFGS-B", "Nelder-Mead", "CG", "SANN"), 
  gr = c("analytical", "numerical"), optim.control = list())
{
  #hessian <- TRUE
  mcall <- match.call()
  method <- match.arg(method)

##FIXME pass xreg already differenced to optim to reduce computations

  bargs <- list(type = "1", mu = 0)
  if (barrier$mu > 0) #&& !is.null(barrier$mu)
  {
    nbargs0 <- names(bargs)
    bargs[nbargs <- names(barrier)] <- barrier
    if (length(nobargs <- nbargs[!nbargs %in% nbargs0]))
      warning("unknown names in 'barrier': ", paste(nobargs, collapse = ", "))
  }

  gr0 <- match.arg(gr)[1]
  if(bargs$mu > 0 && gr0 == "analytical")
    warning("barrier term was ignored.")

  if (is.null(m@cpar))
  {
    fnc <- mloglik.fd

    if(gr0 == "analytical") {
      gr <- function(x, ...) { mloglik.fd.grad(x = x, ...) }
    } else gr <- NULL
  } else {
    fnc <- mcloglik.fd
##FIXME TODO mcloglik with mcloglik.fd

    if(gr0 == "analytical") {
      gr <- function(x, ...) { mcloglik.fd.grad(x = x, ...) }
    } else gr <- NULL
  }

  if (!is.null(m@xreg))
  {
    ##NOTE
    # this overwrites the values m@pars["xreg"], so they 
    # are not used as initial values

    xregnms <- m@ss$xreg
    id.isna <- is.na(match(xregnms, names(c(m@pars, m@nopars))))

    if (all(id.isna))
    {
      dxreg <- m@fdiff(xreg, frequency(m@y))
      fitxreg <- lm(m@diffy ~ dxreg - 1, na.action = na.omit)
      #coef(lm.fit(x= dxreg, y = m@diffy))

      m@pars[xregnms] <- coef(fitxreg)

      if (gr0 == "analytical")
      {
##FIXME TODO this for the other "else" statements

        # element "xreg" is used by "mloglik.fd" and the other 
        # elements in the list are used by "mloglik.fd.deriv",
        # they are defined here to avoid computing these constant 
        # every time the gradient function is called
        # "dxreg" is obtained above whithin this "if" statement
        #dxreg <- model@fdiff(xreg, frequency(m@y))
        xreg <- list(dxreg = dxreg, fft.dxreg = apply(dxreg, 2, fft))
      }

    } else 
    if (all(!id.isna))
    { 
      # required checks wer already made by function "stsm.model"
      # when the input model "m" was created
    } else 
    if (any(id.isna))
    { # some of the xreg coefficients seem to be defined because their names
      # match the column names of "xreg", but for some others there is no match
      # do not make a guess about the input, stop for safety
      stop("some, but not all, names of the parameters in object ", sQuote("m"), 
        " match the column names of ", sQuote("xreg"), 
        ". Either all or neither of the coefficients should be defined in the ", 
        sQuote("stsm"), "input model ", sQuote("m"), ".")
    }

  } else
    xreg <- NULL

  if (method == "L-BFGS-B")
  {
    res <- optim(par = m@pars, fn = fnc, gr = gr, 
      model = m, xreg = xreg, barrier = bargs, inf = inf, method = method, 
      lower = m@lower, upper = m@upper, 
      control = optim.control, hessian = TRUE)
  } else
  {
    res <- optim(par = m@pars, fn = fnc, gr = gr, 
      model = m, xreg = xreg, barrier = bargs, inf = inf, method = method, 
      control = optim.control, hessian = TRUE)
  }

  pars0 <- m@pars
  m@pars <- res$par

  if (!is.null(m@cpar))
  {
##FIXME TODO if xreg!=NULL
    if (is.null(m@ssd)) {
      pg <- Mod(fft(m@diffy))^2 / (2 * pi * length(m@diffy))
    } else
      pg <- m@ssd

    n <- length(m@diffy)
    nh <- n / 2
    pi2 <- 2 * pi
    sgf <- stsm.sgf(m, FALSE, FALSE, FALSE)$sgf
    m@cpar[] <- (pi2 / n) * sum(pg / sgf)

    res$value <- -nh * (log(pi2) + 1) - nh * log(m@cpar) - 0.5 * sum(log(sgf))

    if (barrier$mu != 0) {
      bar <- barrier.eval(m, barrier$type, barrier$mu, 
        FALSE, FALSE)$barrier
      res$value <- res$value - bar
    }
  }

  #if (method == "AB-NM") {
  if (is.null(res$hessian)) {
    std.errors <- NULL
  } else
    std.errors <- sqrt(diag(solve(res$hessian)))

  res2 <- c(list(call = mcall, model = m, 
    init = pars0, pars = m@pars, xreg = xreg, loglik = -res$value,
    convergence = res$convergence, iter = res$counts, message = res$message,
    Mpars = NULL, steps = NULL), list(ls.iter = NULL, ls.counts = NULL),
    list(hessian = res$hessian, std.errors = sqrt(diag(solve(res$hessian))), 
      vcov.type = "optimHessian"))
  class(res2) <- "stsmFit"
  res2
}
