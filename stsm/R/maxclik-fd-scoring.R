
maxclik.fd.scoring <- function(m, step = NULL, 
  information = c("expected", "observed"),
  ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1),
  barrier = list(type = c("1", "2"), mu = 0), 
  control = list(maxit = 100, tol = 0.001, trace = FALSE, silent = FALSE))
{
  if (!is.null(step) && (step <= 0 || !is.numeric(step)))
    stop("'step' must be a numeric positive value.")

  mcall <- match.call()

  mfcn.step <- function(x, m, pd, barrier)
  {
    m <- set.pars(m, get.pars(m) + x * pd)
    mcloglik.fd(x = NULL, model = m, barrier = barrier)
  }

  mgr.step <- function(x, m, pd, barrier)
  {
    m@pars <- m@pars + x * pd

    gr <- mcloglik.fd.deriv(m, TRUE, FALSE, FALSE)$gradient
    gr <- drop(gr %*% pd)
    gr
  }

  cargs <- list(maxit = 100, tol = 0.001, trace = FALSE, silent = TRUE)
  ncargs0 <- names(cargs)
  cargs[ncargs <- names(control)] <- control
  if (length(nocargs <- ncargs[!ncargs %in% ncargs0]))
    warning("unknown names in 'control': ", paste(nocargs, collapse = ", "))
  control <- cargs

  info <- match.arg(information)[1]
  use.gcov <- info == "expected"
  use.hess <- info == "observed"

  lsres <- list(ls.iter = NULL, ls.counts = c("fcn" = 0, "grd" = 0))
  lsintv <- c(0, ls$cap)

  pars0 <- c(m@pars, m@cpar)
  convergence <- FALSE
  iter <- 0

  if (is.null(m@ssd)) {
    pg <- Mod(fft(m@diffy))^2 / (2 * pi * length(m@diffy))
  } else
    pg <- m@ssd

  # spectral generating function (constants)
##FIXME add if is.null(m@sgfc) then define

  Mpars <- if (control$trace) rbind(c(m@cpar, m@pars), 
    matrix(nrow = control$maxit + 1, ncol = length(m@pars) + 1)) else NULL
  steps <- if (control$trace) rep(NA, control$maxit + 2) else NULL

  while (!(convergence || iter > control$maxit))
  {
    tmp <- mcloglik.fd.deriv(m, gradient = TRUE, hessian = use.hess, infomat = use.gcov)

    G <- -tmp$gradient

    M <- switch(info, "expected" = tmp$infomat, 
      "observed" = tmp$hessian)

    if (info == "observed")
    {
      M <- force.defpos(M, 0.001, FALSE)
    }

    pd <- drop(solve(M) %*% G)

    if (is.null(step))
    {
      lsintv[2] <- step.maxsize(get.pars(m), m@lower, m@upper, pd, ls$cap)

      lsout <- switch(ls$type,

      "optimize" = optimize(f = mfcn.step, interval = lsintv, 
        maximum = FALSE, tol = ls$tol, m = m, pd = pd, 
        barrier = barrier),

      "brent.fmin" = Brent.fmin(a = 0, b = lsintv[2], 
        fcn = mfcn.step, tol = ls$tol, m = m, pd = pd, 
        barrier = barrier),

      "wolfe" = linesearch(b = lsintv[2], 
        fcn = mfcn.step, grd = mgr.step, 
        ftol = ls$ftol, gtol = ls$gtol, m = m, pd = pd, 
        barrier = barrier))

      lambda <- lsout$minimum
      lsres$ls.iter <- c(lsres$ls.iter, lsout$iter)
      lsres$ls.counts <- lsres$ls.counts + lsout$counts
    } else
    if (is.numeric(step))
      lambda <- step

    pars.old <- m@pars
    m@pars <- pars.old + lambda * pd

    if (sqrt(sum((pars.old - m@pars)^2)) < control$tol)
    {
      convergence <- TRUE
    }

    iter <- iter + 1

    if (control$trace)
    {
      sg <- stsm.sgf(m, FALSE, FALSE, FALSE)$sgf
      cv <- (2 * pi / length(m@diffy)) * sum(pg / sg)
      Mpars[iter+1,] <- c(cv, m@pars * cv)
      steps[iter] <- lambda
    }
  }

  sgf <- stsm.sgf(m, FALSE, FALSE, FALSE)$sgf
  m@cpar[] <- (2 * pi / length(m@diffy)) * sum(pg / sgf)

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

  val <- -mcloglik.fd(x = NULL, model = m, barrier = barrier)

  vcov.type <- switch(info, "expected" = "information matrix", 
    "observed" = "Hessian")

  res <- c(list(call = mcall,
    init = pars0, pars = m@pars, model = m, loglik = val,
    convergence = convergence, iter = iter, message = "",
    Mpars = Mpars, steps = steps), lsres,
    list(infomat = M, std.errors = sqrt(diag(solve(M))), 
      vcov.type = vcov.type))
  class(res) <- "stsmFit"
  res
}
