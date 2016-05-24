
maxlik.td.scoring <- function(m, step = NULL, 
  KF.args = list(), check.KF.args = TRUE,
  ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1),
  control = list(maxit = 100, tol = 0.001, trace = FALSE, silent = FALSE), 
  debug = FALSE)
{
  if (!is.null(step) && (step <= 0 || !is.numeric(step)))
    stop("'step' must be a numeric positive value.")

  barrier <- list(type = "1", mu = 0)
  mcall <- match.call()

  mll.step <- function(x, m, pd)
  {
    mloglik.td(x = m@pars + x * pd, model = m, 
      KF.version = "KFKSDS", KF.args = KF.args, check.KF.args = FALSE, 
      barrier = list(mu = 0), inf = 99999)
  }

  if (check.KF.args)
    KF.args <- make.KF.args(char2numeric(m), "KFKSDS", KF.args) 

  cargs <- list(maxit = 100, tol = 0.001, trace = FALSE, silent = FALSE)
  ncargs0 <- names(cargs)
  cargs[ncargs <- names(control)] <- control
  if (length(nocargs <- ncargs[!ncargs %in% ncargs0]))
    warning("unknown names in 'control': ", paste(nocargs, collapse = ", "))
  control <- cargs

  # information matrix is used by default
  #info <- "expected"
  #use.IM <- TRUE

  lsres <- list(ls.iter = NULL, ls.counts = c("fnc" = 0, "grd" = 0))
  lsintv <- c(0, ls$cap)

  convergence <- FALSE
  conv.e1 <- FALSE
  iter <- 0

  # regressor variables

  if (!is.null(m@xreg))
  {
    ##NOTE
    # this overwrites the values m@pars["xreg"], so they 
    # are not used as initial values
    dxreg <- m@fdiff(m@xreg, frequency(m@y))
    fitxreg <- lm(m@diffy ~ dxreg - 1, na.action = na.omit)
    m@pars[m@ss$xreg] <- coef(fitxreg)
  } #else xregcoefs <- NULL

  pars0 <- m@pars

  # storage matrices for tracing information

  if (control$trace) {
    Mpars <-  rbind(pars0, 
      matrix(nrow = control$maxit + 1, ncol = length(pars0))) 
  } else Mpars <- NULL
  steps <- if (control$trace) rep(NA, control$maxit + 2) else NULL

  # begin iterative process

j1 <- 1

  while (!(convergence || iter > control$maxit))
  {
    tmp <- mloglik.td.deriv(model = m, gradient = TRUE, 
      infomat = TRUE, KF.args = KF.args, version = "1", kfres = NULL)
##FIXME add option argument convergence = c(0.001, length(model@y)

    # gradient

    G <- -tmp$gradient

    # information matrix

    IM <- tmp$infomat

if (debug)
{
  eg <- eigen(IM, only.values = TRUE)$values
  if (min(eg) <= 0) {
    warning("IM is not a positive definite matrix.")
    print(eg)
  }
}

    # direction vector

    pd <- drop(solve(IM) %*% G)

    # step size (choose the step that maximizes the increase in 
    # the log-likelihood function given the direction vector 'pd')

    if (is.null(step))
    {
      lsintv[2] <- step.maxsize(m@pars, m@lower, m@upper, pd, ls$cap)

      lsout <- switch(ls$type,

      "optimize" = optimize(f = mll.step, interval = lsintv, 
        maximum = FALSE, tol = ls$tol, m = m, pd = pd),

      "brent.fmin" = Brent.fmin(a = 0, b = lsintv[2], 
        fcn = mll.step, tol = ls$tol, m = m, pd = pd))

      lambda <- lsout$minimum
      lsres$ls.iter <- c(lsres$ls.iter, lsout$iter)
      lsres$ls.counts <- lsres$ls.counts + lsout$counts

    } else
    if (is.numeric(step))
      lambda <- step

    # update (new guess)

    pars.old <- m@pars
    pars.new <- pars.old + lambda * pd
    m <- set.pars(m, pars.new)

    # stopping criteria

    if (sqrt(sum((pars.old - m@pars)^2)) < control$tol)
    {
      convergence <- TRUE
    } else
    if (sum(abs(pars.old - m@pars) > control$tol) == 1)
      j1 <- j1 + 1
    if (j1 == 10) {
      convergence <- TRUE
      conv.e1 <- TRUE
      if (!control$silent)
        warning(paste("Possible convergence problem.\n",
        "Over 10 iterations, failure to concergence was caused by just 1 parameter",
        "(not necessarily the same parameter)."))
    }

    # trace

    iter <- iter + 1

    if (control$trace)
    {
      Mpars[iter+1,] <- pars.new
      steps[iter] <- lambda
    }

if (debug)
{
  val <- logLik(object = m, domain = "time", 
    td.args = list(P0cov = FALSE, t0 = 1, KF.version = "KFKSDS"), 
    barrier = list(mu = 0), inf = 99999)
  cat(paste("\niter =", iter, "logLik =", round(val, 4), "\n"))
  print(get.pars(m))
}

if (debug && !is.null(m@lower) && !is.null(m@upper))
{
  check.bounds(m)
}
  }

  if (!control$silent && !convergence)
    warning(paste("Possible convergence problem.",
    "Maximum number of iterations reached."))

  if (control$trace)
  {
    Mpars <- na.omit(Mpars)
    attr(Mpars, "na.action") <- NULL
    steps <- na.omit(steps)
    attr(steps, "na.action") <- NULL
  }

  # output

  val <- -mloglik.td(#x = NULL, 
    model = m, KF.version = "KFKSDS", 
    KF.args = KF.args, check.KF.args = FALSE, 
    barrier = barrier, inf = 99999)

  if (convergence) {
    convergence <- "yes"
  } else
    convergence <- "maximum number of iterations was reached"
  if (conv.e1)
    convergence <- "possible convergence problem"
    #"Over 10 iterations, failure to concergence was caused by just 1 parameter"
    #"(not necessarily the same parameter)"

  IM <- mloglik.td.deriv(model = m, gradient = FALSE, 
      infomat = TRUE, KF.args = KF.args, version = "1", kfres = NULL)$infomat

  res <- c(list(call = mcall, model = m, 
    init = pars0, pars = m@pars, loglik = val,
    convergence = convergence, iter = iter, message = "",
    Mpars = Mpars, steps = steps), lsres, 
    list(infomat = IM, std.errors = sqrt(diag(solve(IM))), 
      vcov.type = "information matrix"))
  class(res) <- "stsmFit"
  res
}
