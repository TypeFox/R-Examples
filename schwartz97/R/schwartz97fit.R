## <---------------------------------------------------------------------->
fit.schwartz2f <- function(data, ttm, deltat = 1 / 260,
                           s0 = data[1,1],
                           delta0 = 0,
                           mu = 0.1, sigmaS = 0.3,
                           kappa = 1, alpha = 0, sigmaE = 0.3,
                           rho = 0.7, lambda = 0,
                           meas.sd = rep(0.1, ncol(data)),
                           opt.pars = c(s0 = FALSE, delta0 = FALSE, mu = TRUE,
                             sigmaS = TRUE, kappa = TRUE, alpha = TRUE,
                             sigmaE = TRUE, rho = TRUE, lambda = FALSE),
                           opt.meas.sd = c("scalar", "all", "none"),
                           r = 0.03, silent = FALSE, ...)
{
  call <- match.call()

  opt.meas.sd <- match.arg(opt.meas.sd)

  ## ---------------------------------------------------------------------------
  ## Internal function to compute the log-likelihood
  log.likelihood.2f <- function(thetaOpt, thetaConst, thetaNames,
                                data, ttm, deltat, r, d, n,
                                meas.sd, opt.meas.sd, silent)
    {
      ## Transformations
      if(!is.na(thetaOpt["rho"])){
        thetaOpt["rho"] <- 2 * atan(thetaOpt["rho"]) / pi
      }
      for(name in setdiff(names(thetaOpt), c("delta0", "mu", "alpha", "lambda", "rho"))){
        thetaOpt[name] <- exp(thetaOpt[name])
      }

      theta <- c(thetaOpt, thetaConst)
      theta <- theta[thetaNames]

      if(opt.meas.sd == "scalar"){
        gg <- theta["meas.sd1"] * meas.sd
      }else if(opt.meas.sd == "all"){
        gg <- theta[10:(10 + d - 1)]
      }else{                            #opt.pars.sd == "none"
        gg <- meas.sd
      }

      ## Build State Space Elements
      stateSpace <- .state.space.2f(y = data, ttm = ttm,
                                    deltat = deltat,
                                    x0 = log(theta["s0"]), delta0 = theta["delta0"],
                                    kappa = theta["kappa"], mu = theta["mu"],
                                    alpha = theta["alpha"], lambda = theta["lambda"],
                                    sigmaS = theta["sigmaS"], sigmaE = theta["sigmaE"],
                                    rho = theta["rho"],
                                    gg = gg, r = r, d = d, n = n)

      logLikelihood <- fkf(a0 = stateSpace$a0,
                           P0 = stateSpace$P0,
                           Tt = stateSpace$Tt,
                           dt = stateSpace$dt,
                           HHt = stateSpace$HHt,
                           yt = stateSpace$yt,
                           Zt = stateSpace$Zt,
                           ct = stateSpace$ct,
                           GGt = stateSpace$GGt)$logLik

      n.iter <- nrow(theta.backup)
      rel.tol <- abs((logLikelihood - theta.backup[n.iter, "logLik"])/
                     theta.backup[n.iter, "logLik"])
      abs.tol <- abs(logLikelihood - theta.backup[n.iter, "logLik"])

      if(!silent){
        cat("--> i: ", n.iter + 1,
            "; logL: ",sprintf("%.4E", logLikelihood), 
            "; rel.tol: ",sprintf("%.2E", rel.tol), 
            "; abs.tol: ",sprintf("%.2E", abs.tol), "; ",
            paste(names(thetaOpt), ": ", sprintf("%.2E", thetaOpt), "; ", collapse = "", sep = ""),
            "\n", sep = "")
      }

      ## if(!is.na(thetaOpt["rho"]))
      ##   {
      ##     thetaOpt["rho"] <- tan(thetaOpt["rho"] * pi / 2)
      ##   }

      theta.backup <<- rbind(theta.backup, c(logLikelihood, rel.tol, abs.tol, thetaOpt))


      return(-logLikelihood)

    }
  ## ---------------------------------------------------------------------------

  log.data <- log(data)

  d <- ncol(data)                      # Dimension of the observations
  n <- nrow(data)                      # Number of observations

  ## Check whether arguments are feasible
  if(any(data < 0, na.rm = TRUE)){
    stop("All elements of 'data' must be positive!")
  }
  if(any(ttm < 0, na.rm = TRUE)){
    stop("All elements of 'ttm' must be positive!")
  }

  if(!(is.logical(opt.pars) & (length(opt.pars) == 9))){
    stop("'opt.pars' must be of class 'logical' and of length 9!\n")
  }
  if(!is(data, "matrix")){
    stop("'data' must be a matrix!")
  }
  if(!is(ttm, "matrix")){
    stop("'ttm' must be a matrix!")
  }
  ## if(any(!is.finite(ttm))){
  ##   stop("'ttm' contains non-finite values!")
  ## }
  if(any(meas.sd < 0)){
    stop("Elements of 'meas.sd' must not be smaller than 0!")
  }
  if(length(meas.sd) != d){
    stop("length(meas.sd) must be of the same dimension as 'data'!")
  }
  if(any(c(sigmaS, sigmaE, s0) <= 0)){
    stop("'sigmaS', 'sigmaE', and 's0' must be greater than 0!")
  }
  if(rho < -1 | rho > 1){
    stop("'rho' must be in [-1, 1]!")
  }

  ## Initialization
  thetaNames <- c("s0", "delta0", "mu", "sigmaS", "kappa",
                  "alpha", "sigmaE", "rho", "lambda",
                  paste("meas.sd", 1:d, sep = ""))

  ## Initial values must be scalars. Check it:
  .check.lengths(s0 = s0, delta0 = delta0, mu = mu, sigmaS = sigmaS, kappa = kappa,
                 alpha = alpha, sigmaE = sigmaE, rho = rho, lambda = lambda, r = r)

  theta <- c(s0, delta0, mu, sigmaS, kappa,
             alpha, sigmaE, rho, lambda, meas.sd)
  names(theta) <- thetaNames

  if(opt.meas.sd == "scalar"){
    opt.pars <- c(opt.pars, TRUE, rep(FALSE, d - 1))
    theta["meas.sd1"] <- 1
  }else if(opt.meas.sd == "all"){
    opt.pars <- c(opt.pars, rep(TRUE, d))
  }else{                                #opt.pars.sd == "none"
    opt.pars <- c(opt.pars, rep(FALSE, d))
  }

  names(opt.pars) <- thetaNames

  thetaOpt <- theta[opt.pars]
  thetaConst <- theta[!opt.pars]

  theta.backup <- rbind(c(NA, NA, NA, thetaOpt))
  colnames(theta.backup) <- c("logLik", "rel.tol", "abs.tol", names(thetaOpt))

  ## Transformations
  if(opt.pars["rho"]){
    thetaOpt["rho"] <- tan(thetaOpt["rho"] * pi / 2)
  }
  for(name in setdiff(names(thetaOpt), c("delta0", "mu", "alpha", "lambda", "rho"))){
    thetaOpt[name] <- log(thetaOpt[name])
  }
  mle <- try(optim(thetaOpt, fn = log.likelihood.2f,
                   thetaConst = thetaConst, thetaNames = thetaNames,
                   data = log.data, ttm = ttm, deltat = deltat,
                   r = r, d = d, n = n, meas.sd = meas.sd,
                   opt.meas.sd = opt.meas.sd,
                   silent = silent, ...))

  if(class(mle) == "try-error"){
    convergence <- -1
    n.iter <- nrow(theta.backup)
    message <- as.character(mle)
    thetaOpt <- theta.backup[nrow(theta.backup), -(1:3)]
  }else{
    convergence <- mle$convergence
    n.iter <- mle$counts[1]
    message <- ""
    thetaOpt <- mle$par
  }

  if(opt.pars["rho"]){
    thetaOpt["rho"] <- 2 * atan(thetaOpt["rho"]) / pi
  }
  for(name in setdiff(names(thetaOpt), c("delta0", "mu", "alpha", "lambda", "rho"))){
    thetaOpt[name] <- exp(thetaOpt[name])
  }

  theta <- c(thetaOpt, thetaConst)
  theta <- theta[thetaNames]

  if(opt.meas.sd == "scalar"){
    gg <- theta["meas.sd1"] * meas.sd
  }else if(opt.meas.sd == "all"){
    gg <- theta[10:(10 + d - 1)]
  }else{                                #opt.pars.sd == "none"
    gg <- meas.sd
  }

  stateSpace <- .state.space.2f(y = data, ttm = ttm,
                                deltat = deltat,
                                x0 = log(theta["s0"]), delta0 = theta["delta0"],
                                kappa = theta["kappa"], mu = theta["mu"],
                                alpha = theta["alpha"], lambda = theta["lambda"],
                                sigmaS = theta["sigmaS"], sigmaE = theta["sigmaE"],
                                rho = theta["rho"],
                                gg = gg,
                                r = r, d = d, n = n)

  filtered.ts <- fkf(a0 = stateSpace$a0,
                     P0 = stateSpace$P0,
                     Tt = stateSpace$Tt,
                     dt = stateSpace$dt,
                     HHt = stateSpace$HHt,
                     yt = stateSpace$yt,
                     Zt = stateSpace$Zt,
                     ct = stateSpace$ct,
                     GGt = stateSpace$GGt)

  ##   state <- cbind(S = exp(filtered.ts$att[1,]), delta = filtered.ts$att[2,])

  return(new("schwartz2f.fit",
             call = call,
             s0 = unname(theta["s0"]),
             delta0 = unname(theta["delta0"]),
             mu = unname(theta["mu"]),
             sigmaS = unname(theta["sigmaS"]),
             kappaE = unname(theta["kappa"]),
             alpha = unname(theta["alpha"]),
             sigmaE = unname(theta["sigmaE"]),
             rhoSE = unname(theta["rho"]),
             n.iter = unname(n.iter),
             llh = unname(filtered.ts$logLik),
             converged = convergence == 0,
             error.code= unname(convergence),
             error.message = message,
             fitted.params = opt.pars,
             trace.pars = theta.backup,
             r = unname(r),
             alphaT = unname(theta["alpha"] - theta["lambda"] / theta["kappa"]),
             lambdaE = unname(theta["lambda"]),
             meas.sd = gg,
             deltat = deltat))
}
