
maxlik.td.optim <- function(m, 
  KF.version = eval(formals(KFKSDS::KalmanFilter)$KF.version),
  KF.args = list(), check.KF.args = TRUE,
  barrier = list(type = c("1", "2"), mu = 0), inf = 99999, 
  method = c("BFGS", "L-BFGS-B", "Nelder-Mead", "CG", "SANN", "AB-NM"),
  gr = c("numerical", "analytical"), optim.control = list())
{
  #hessian <- TRUE
  mcall <- match.call()
  method <- match.arg(method)

  bargs <- list(type = "1", mu = 0)
  if (barrier$mu > 0) #&& !is.null(barrier$mu)
  {
    nbargs0 <- names(bargs)
    bargs[nbargs <- names(barrier)] <- barrier
    if (length(nobargs <- nbargs[!nbargs %in% nbargs0]))
      warning("unknown names in 'barrier': ", paste(nobargs, collapse = ", "))
  }

  KF.version <- match.arg(KF.version)[1]
  gr0 <- match.arg(gr)[1]
  if (gr0 == "analytical" && KF.version != "KFKSDS")
  {
    warning("analytical gradient is not available for 'KF.version = ", 
      KF.version, "'. Changed to 'KF.version = KFKSDS'.")
    KF.version <- "KFKSDS"
  }

  if (check.KF.args)
    KF.args <- make.KF.args(char2numeric(m), KF.version, KF.args) 

  if (gr0 == "analytical")
  {
    if (!is.null(m@cpar)) {
      warning("numerical gradient was used (analytical gradient is not available for ", 
        sQuote("KFconvar"), ").")
      gr <- NULL
    } else gr <- function(x, ...) { mloglik.td.grad(x = x, ...) }
  } else gr <- NULL

  if (!is.null(m@xreg))
  {
    ##NOTE
    # this overwrites the values m@pars["xreg"], so they 
    # are not used as initial values

    xregnms <- m@ss$xreg
    id.isna <- is.na(match(xregnms, names(c(m@pars, m@nopars))))

    if (all(id.isna))
    {
      dxreg <- m@fdiff(m@xreg, frequency(m@y))
      fitxreg <- lm(m@diffy ~ dxreg - 1, na.action = na.omit)
      #coef(lm.fit(x= dxreg, y = m@diffy))

      m@pars[xregnms] <- coef(fitxreg)
    } else 
    if (all(!id.isna))
    {
      # required checks were already made by function "stsm.model"
      # when the input model "m" was created
    } else 
    if (any(id.isna))
    {
      stop("some, but not all, names of the parameters in object ", sQuote("m"), 
        " match the column names of ", sQuote("xreg"), 
        ". Either all or neither of the coefficients should be defined in the ", 
        sQuote("stsm"), "input model ", sQuote("m"), ".")
    }
  }

  if (method == "L-BFGS-B")
  {
    res <- optim(par = m@pars, fn = mloglik.td, gr = gr, 
      model = m, KF.version = KF.version, KF.args = KF.args, 
      barrier = bargs, inf = inf, #load.package = FALSE, 
      method = method, lower = m@lower, upper = m@upper, 
      control = optim.control, hessian = TRUE)
  } else 
  if (method == "BFGS") {
    res <- optim(par = m@pars, fn = mloglik.td, gr = gr, 
      model = m, KF.version = KF.version, KF.args = KF.args, 
      barrier = bargs, inf = inf, #load.package = FALSE, 
      method = method, 
      control = optim.control, hessian = TRUE)
  } else
  if (method == "AB-NM")
  {
    res <- constrOptim(theta = m@pars, f = mloglik.td, grad = NULL,
      model = m, KF.version = KF.version, KF.args = KF.args,
      ui = diag(length(m@pars)), ci = m@lower, mu = 1e-04, 
      method = "Nelder-Mead")
  }

  pars0 <- c(m@pars, m@cpar)
  m@pars <- res$par

  if (!is.null(m@cpar))
  {
##FIXME TODO if xreg!=NULL
    tmp <- KFconvar(m, P0cov = KF.args$P0cov, debug = TRUE)
    m@cpar[] <- tmp$cpar
    res$value <- -tmp$mll

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
    init = pars0, pars = m@pars, loglik = -res$value,
    convergence = res$convergence, iter = res$counts, message = res$message,
    Mpars = NULL, steps = NULL), list(ls.iter = NULL, ls.counts = NULL),
    list(hessian = res$hessian, std.errors = std.errors, vcov.type = "optimHessian"))
  class(res2) <- "stsmFit"
  res2
}
