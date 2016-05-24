
##NOTE
#conv.e1 is not considered here, it is tracked in maxlik.td.scoring()

maxlik.fd.scoring <- function(m, step = NULL, 
  information = c("expected", "observed", "mix"),
  ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1),
  barrier = list(type = c("1", "2"), mu = 0), 
  control = list(maxit = 100, tol = 0.001, trace = FALSE, silent = FALSE), 
  debug = FALSE)
{
  if (!is.null(step) && (step <= 0 || !is.numeric(step)))
    stop("'step' must be a numeric positive value.")

  mcall <- match.call()

  mll.step <- function(x, m, pd, barrier)
  {
    mloglik.fd(x = m@pars + x * pd, model = m, 
      barrier = barrier) # by default: inf = 99999
  }

  dmll.step <- function(x, m, pd, xreg, barrier)
  {
    m@pars <- m@pars + x * pd
    gr <- mloglik.fd.deriv(model = m, xreg = xreg, gradient = TRUE,
      hessian = FALSE, infomat = FALSE, modcovgrad = FALSE,
      barrier = barrier, version = "2")$gradient
    drop(gr %*% pd)
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
  use.mix <- info == "mix"

  lsres <- list(ls.iter = NULL, ls.counts = c("fnc" = 0, "grd" = 0))
  lsintv <- c(0, ls$cap)

  convergence <- FALSE
  iter <- 0

  # periodogram

  if (is.null(m@ssd)) {
    pg <- Mod(fft(m@diffy))^2 / (2 * pi * length(m@diffy))
    #pi2.pg <- pi2 * pg
  } else
    pg <- m@ssd

# spectral generating function (constants)
##FIXME add if is.null(m@sgfc) then define sgfc (see repository previous version)
#this would save some computations if m@sgfc is not defined in the input model

  # regressor variables

  if (!is.null(m@xreg))
  {
    ##NOTE
    # this overwrites the values m@pars["xreg"], so they 
    # are not used as initial values
    dxreg <- m@fdiff(m@xreg, frequency(m@y))
    fitxreg <- lm(m@diffy ~ dxreg - 1, na.action = na.omit)
    m@pars[m@ss$xreg] <- coef(fitxreg)

    # constant terms to be passed to mloglik.fd.deriv() in the loop below
    xreg <- list(dxreg = dxreg, fft.dxreg = fft.dxreg <- apply(dxreg, 2, fft))
  } else
    xreg <- NULL

  pars0 <- m@pars

  # storage matrices for tracing information

  if (control$trace) {
    Mpars <-  rbind(pars0, 
      matrix(nrow = control$maxit + 1, ncol = length(pars0))) 
  } else Mpars <- NULL
  steps <- if (control$trace) rep(NA, control$maxit + 2) else NULL

  # begin iterative process

  while (!(convergence || iter > control$maxit))
  {
    tmp <- mloglik.fd.deriv(m, xreg = xreg, gradient = TRUE, 
      hessian = use.hess, infomat = use.gcov, modcovgrad = use.mix,
      barrier = barrier, version = "2")

    # gradient

    G <- -tmp$gradient

    # information matrix

    M <- switch(info, "expected" = tmp$infomat, 
      "observed" = tmp$hessian, "mix" = tmp$modcovgrad)

    if (info == "observed" || info == "mix")
    {
      M <- force.defpos(M, 0.001, FALSE)
    }

    # direction vector

    pd <- drop(solve(M) %*% G)

    # step size (choose the step that maximizes the increase in 
    # the log-likelihood function given the direction vector 'pd')

    if (is.null(step))
    {
      lsintv[2] <- step.maxsize(m@pars, m@lower, m@upper, pd, ls$cap)

      lsout <- switch(ls$type,

      "optimize" = optimize(f = mll.step, interval = lsintv, 
        maximum = FALSE, tol = ls$tol, m = m, pd = pd, 
        barrier = barrier),

      "brent.fmin" = Brent.fmin(a = 0, b = lsintv[2], 
        fcn = mll.step, tol = ls$tol, m = m, pd = pd, 
        barrier = barrier),

      "wolfe" = linesearch(b = lsintv[2], 
        fcn = mll.step, grd = dmll.step, 
        ftol = ls$ftol, gtol = ls$gtol, m = m, pd = pd, xreg = xreg,
        barrier = barrier))

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

    if (sqrt(sum((pars.old - pars.new)^2)) < control$tol)
    {
      convergence <- TRUE
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
  val <- logLik(object = m, domain = "frequency", barrier = barrier)
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

  val <- -mloglik.fd(#x = NULL, 
    model = m, barrier = barrier, inf = 99999)

  if (convergence) {
    convergence <- "yes"
  } else
    convergence <- "maximum number of iterations was reached"

  vcov.type <- switch(info, "expected" = "information matrix", 
    "observed" = "Hessian", "mix" = "modified outer product of the gradient")

##FIXME see rename "Dmat"

  res <- c(list(call = mcall, model = m, 
    init = pars0, pars = m@pars, xreg = xreg, loglik = val, 
    convergence = convergence, iter = iter, message = "",
    Mpars = Mpars, steps = steps), lsres, 
    list(Dmat = M, std.errors = sqrt(diag(solve(M))), vcov.type = vcov.type))
  class(res) <- "stsmFit"
  res
}
