endorse <- function(Y,
                    data,
                    data.village = NA,
                    village = NA,
                    treat = NA,
                    na.strings = 99,
                    identical.lambda = TRUE,
                    covariates = FALSE,
                    formula.indiv = NA,
                    hierarchical = FALSE,
                    formula.village = NA,
                    x.start = 0,
                    s.start = 0,
                    beta.start = 1,
                    tau.start = NA,
                    lambda.start = 0,
                    omega2.start = .1,
                    theta.start = 0,
                    phi2.start = .1,
                    kappa.start = 0,
                    psi2.start = 1,
                    delta.start = 0,
                    zeta.start = 0,
                    rho2.start = 1,
                    mu.beta = 0,
                    mu.x = 0,
                    mu.theta = 0,
                    mu.kappa = 0,
                    mu.delta = 0,
                    mu.zeta = 0,
                    precision.beta = 0.04,
                    precision.x = 1,
                    precision.theta = 0.04,
                    precision.kappa = 0.04,
                    precision.delta = 0.04,
                    precision.zeta = 0.04,
                    s0.omega2 = 1,
                    nu0.omega2 = 10,
                    s0.phi2 = 1,
                    nu0.phi2 = 10,
                    s0.psi2 = 1,
                    nu0.psi2 = 10,
                    s0.sig2 = 1,
                    nu0.sig2 = 400,
                    s0.rho2 = 1,
                    nu0.rho2 = 10,
                    MCMC = 20000,
                    burn = 1000,
                    thin = 1,
                    mh = TRUE,
                    prop = 0.001,
                    x.sd = TRUE,
                    tau.out = FALSE,
                    s.out = FALSE,
                    omega2.out = TRUE,
                    phi2.out = TRUE,
                    psi2.out = TRUE, 
                    verbose = TRUE,
                    seed.store = FALSE,
                    update = FALSE,
                    update.start = NULL) {

  mda = FALSE

  ## if (!identical.lambda & hierarchical) {
  ##   stop("Options are not consistent. If 'identical.lambda = TRUE', 'hierarchical' must be set at
  ##        'FALSE.'")
  ## }

  if (covariates == FALSE) {
    formula.indiv <- ~ 1
  }

  cov.mat <- model.matrix(formula.indiv, data)
  M <- ncol(cov.mat)

  var.names.indiv <- colnames(cov.mat)

  #############################################
  ## NEED TO MODIFY
  #############################################
  data <- data[complete.cases(cov.mat),] ## do NOT work. NOT DROP NA obs
  cov.mat <- cov.mat[complete.cases(cov.mat), ]
  
  N <- nrow(data)
  #############################################
  J <- length(Y)

  if (!is.na(treat[1])) {
    if (class(treat) != "matrix" | dim(treat)[1] != N | dim(treat)[2] != J) {
      stop(paste("'treat' must be a ", N, "-by-", J, " matrix.", sep = ""))
    }
  }

  if (hierarchical) {
    village <- as.integer(as.factor(eval(parse(text = paste("data$", village,
                                                 sep = ""))))) - 1
    G <- length(unique(village))
    if (!is.na(data.village[1, 1])) {
      cov.village.mat <- model.matrix(formula.village, data.village)
      P <- ncol(cov.village.mat)
      var.names.village <- colnames(cov.village.mat)
    }
  } else {
    G <- 1
    P <- 1
    cov.village.mat <- double(1)
    village <- rep(0, times = N)
  }
  
  response <- matrix(NA, nrow = N, ncol = J)
  temp.Qs <- paste("Y$Q", 1:J, sep ="")

  if (is.na(treat[1])) {
    endorse <- matrix(0, nrow = N, ncol = J)
  } else {
    endorse <- treat
  }

  K <- rep(NA, times = J)
  
  if (is.na(na.strings[1])) {
    for (i in 1:J) {
      temp.group <- eval(parse(text = paste("length(", temp.Qs[i], ")", sep ="")))
      K[i] <- temp.group
      
      for (j in 1:temp.group) {
        varname <- eval(parse(text = paste(temp.Qs[i], "[j]", sep = "")))
        response[, i] <- eval(parse(text = paste("ifelse(!is.na(data$", varname,
                                      "), data$", varname, ", response[, i])",
                                      sep = "")))

        if (is.na(treat[1])) {
          endorse[, i] <- eval(parse(text = paste("ifelse(!is.na(data$", varname,
                                       "), j - 1, endorse[, i])",
                                       sep = "")))
        }
      }
    }
  } else {
    for (i in 1:J) {
      temp.group <- eval(parse(text = paste("length(", temp.Qs[i], ")", sep ="")))
      K[i] <- temp.group

      for (j in 1:temp.group) {
        varname <- eval(parse(text = paste(temp.Qs[i], "[j]", sep = "")))
        response[, i] <- eval(parse(text = paste("ifelse(!is.na(data$", varname,
                                      ") & !(data$", varname,
                                      " %in% na.strings), data$", varname,
                                      ", response[, i])", sep = "")))

        if (is.na(treat[1])) {
          endorse[, i] <- eval(parse(text = paste("ifelse(!is.na(data$", varname,
                                       "), j - 1, endorse[, i])",
                                       sep = "")))
        }
      }
    }
  }

  for (i in 1:J){
    response[, i] <- as.integer(as.ordered(response[, i]))
  }

  L <- apply(response, 2, max, na.rm = TRUE)
  max.L <- max(L)
  
  response <- response - 1

  for (i in 1:J) {
    response[, i] <- ifelse(is.na(response[, i]), -1, response[, i])
  }
  
  if (is.na(treat[1])) {
    K <- max(K) - 1
  } else {
    K <- max(treat)
  }

  if (update) {
    beta.start <- update.start$beta.start
    
    tau.start <- matrix(-99, nrow = J, ncol = max.L)
    for (j in 1:J) {
      tau.start[j, 1:(max.L - 1)] <- update.start$tau.start[((max.L - 1) *
                                                             (j - 1) + 1):((max.L -
                                                                            1) * j)]
      tau.start[j, L[j]] <- max(tau.start[j, ]) + 1000
    }

    x.start <- update.start$x.start
    s.start <- update.start$s.start
    lambda.start <- update.start$lambda.start
    omega2.start <- update.start$omega2.start
    theta.start <- update.start$theta.start
    phi2.start <- update.start$phi2.start
    if (covariates) delta.start <- update.start$delta.start
    if (hierarchical) {
      kappa.start <- update.start$kappa.start
      psi2.start <- update.start$psi2.start
      zeta.start <- update.start$zeta.start
      rho2.start <- update.start$rho2.start
    }

    .Random.seed <- update.start$seed
    
  } else {
    
    if (is.na(tau.start[1])) {
      tau.start <- matrix(-99, nrow = J, ncol = max.L)

      for (j in 1:J){
        temp.tau <- seq(from = 0, to = .5 * (L[j] - 2), by = .5)
        for (i in 1:(L[j] - 1)) {
          tau.start[j, i] <- temp.tau[i]
        }
        tau.start[j, L[j]] <- max(temp.tau) + 1000
      }
    } else if (class(tau.start) != "matrix" | dim(tau.start)[1] != J | dim(tau.start)[2] != max.L) {
      stop(paste("Incorrect input for tau.start. It should be a ",
                 J, "-by-", max.L, " matrix.", sep = ""))
    }
    
    if (length(x.start) == 1) {
      x.start <- rep(x.start, times = N)
    } else if (length(x.start) != N) {
      stop(paste("The length of 'x.start' should be", N))
    }

    if (length(s.start) == 1) {
      s.start <- as.numeric(matrix(s.start, nrow = N, ncol = J))
      s.start[as.integer(endorse) == 0] <- 0
    } else if (class(s.start) != "matrix" | dim(s.start)[1] != N | dim(s.start)[2] != J) {
      stop(paste("'s.start' should be a ", N, "-by-", J, " matrix.", sep = ""))
    } else {
      s.start <- as.double(s.start)
      s.start[as.integer(endorse) == 0] <- 0
    }

    if (length(beta.start) == 1) {
      beta.start <- matrix(beta.start, nrow = J, ncol = 2)
    } else if (class(beta.start) != "matrix" | dim(beta.start)[1] != J | dim(beta.start)[2] != 2) {
      stop(paste("Incorrect input for 'beta.start'. It should be a ",
                 J, "-by-", 2, " matrix.", sep = ""))
    }

    if (length(lambda.start) == 1) {
      lambda.start <- matrix(lambda.start, nrow = G + M, ncol = J * K)
    } else if (covariates & !identical.lambda & !hierarchical &
               (class(lambda.start) != "matrix" | dim(lambda.start)[1] != M |
                dim(lambda.start)[2]  != J * K)) {
      stop(paste("Incorrect input for lambda.start. It should be a ",
                 M, "-by-", J * K, " matrix.", sep = ""))
    } else if (covariates & !identical.lambda & hierarchical &
               (class(lambda.start) != "matrix" | dim(lambda.start)[1] != G + M - 1 |
                dim(lambda.start)[2]  != J * K)) {
      stop(paste("Incorrect input for lambda.start. It should be a ",
                 M - 1, "-by-", J * K, " matrix.", sep = ""))
    } else if (!covariates & !identical.lambda & (class(lambda.start) != "matrix" |
                                                  dim(lambda.start)[1] != G |
                                                  dim(lambda.start)[2] != J * K)) {
      stop(paste("Incorrect input for lambda.start. It should be a ", G, "-by-",
                 J * K, " matrix.", sep = ""))
    } else if (covariates & identical.lambda & hierarchical &
               (class(lambda.start) != "matrix" | dim(lambda.start)[1] != G + M - 1 |
                dim(lambda.start)[2] != K)) {
      stop(paste("Incorrect input for lambda.start. It should be a ",
                 G + M - 1, "-by-", K, " matrix.", sep = ""))
    } else if (!covariates & identical.lambda & hierarchical & length(lambda.start) != G * K) {
      stop(paste("Incorrect input for lambda.start. Its length should be ", G * K,
                 sep = ""))
    }

    if (length(omega2.start) == 1) {
      omega2.start <- rep(omega2.start, times = J * K)
    } else if (identical.lambda & length(omega2.start) != K) {
      stop(paste("Incorrect input for omega2.start. Its length should be ", K,
                 sep = ""))
    } else if (!identical.lambda & length(omega2.start) != J * K) {
      stop(paste("Incorrect input for omega2.start. Its length should be ", J * K,
                 sep = ""))
    }

    ## check the dim of starting values for theta
    if (length(theta.start) == 1) {
      if (hierarchical) {
        theta.start <- matrix(theta.start, nrow = M - 1 + P, ncol = K)
      } else {
        theta.start <- matrix(theta.start, nrow = M, ncol = K)
      }
    } else if (hierarchical & !identical.lambda &
               (class(theta.start) != "matrix" | dim(theta.start)[1] != M - 1 + P |
                dim(theta.start)[2] != K)) {
      stop(paste("'theta.start' should be a ", M - 1 + P, "-by-", K, " matrix.", sep = ""))
    } else if (!hierarchical & !identical.lambda &
               (class(theta.start) != "matrix" | dim(theta.start)[1] != M |
                dim(theta.start)[2] != K)) {
      stop(paste("'theta.start' should be a ", M, "-by-", K, " matrix.", sep = ""))
    }


    ## check the dim of starting values for phi2
    if (length(phi2.start) == 1) {
      if (hierarchical) {
        phi2.start <- rep(phi2.start, times = K * (M - 1 + P))
      } else {
        phi2.start <- rep(phi2.start, times = K * M)
      }
    } else if (hierarchical & !identical.lambda &
               (class(phi2.start) != "matrix" | dim(phi2.start)[1] != M - 1 + P |
                dim(phi2.start)[2] != K)) {
      stop(paste("'phi2.start' should be a ", M - 1 + P, "-by-", K, " matrix.", sep = ""))
    } else if (!hierarchical & !identical.lambda &
               (class(phi2.start) != "matrix" | dim(phi2.start)[1] != M |
                dim(phi2.start)[2] != K)) {
      stop(paste("'phi2.start' should be a ", M, "-by-", K, " matrix.", sep = ""))
    }


    ## check the dim of starting values for kappa, psi2, zeta
    if (hierarchical) {
      ## kappa
      if (length(kappa.start) == 1) {
        if (identical.lambda) {
          kappa.start <- rep(kappa.start, times = K * P)
        } else {
          kappa.start <- rep(kappa.start, times = J * K * P)
        }
      } else if (identical.lambda & (class(kappa.start) != "matrix" | dim(kappa.start)[1] != P |
                                     dim(kappa.start)[2] != K)) {
        stop(paste("Incorrect input for kappa.start. It should be a ", P, "-by", K, " matrix.",
                   sep = ""))
      } else if (!identical.lambda & (class(kappa.start) != "matrix" | dim(kappa.start)[1] != P |
                                      dim(kappa.start)[2] != J * K)) {
        stop(paste("Incorrect input for kappa.start. It should be a ", P, "-by", J * K, " matrix.",
                   sep = ""))
      }

      ## psi2
      if (length(psi2.start == 1)) {
        if (identical.lambda) {
          psi2.start <- rep(psi2.start, times = K)
        } else {
          psi2.start <- rep(psi2.start, times = J*K)
        }
      } else if (identical.lambda & length(psi2.start) != K){
        stop(paste("'psi2.start' should be a vector of length ", K, ".", sep = ""))
      } else if (!identical.lambda & length(psi2.start) != J * K){
        stop(paste("'psi2.start' should be a vector of length ", J * K, ".", sep = ""))
      }

      ## zeta
      if (length(zeta.start) == 1) {
        zeta.start <- rep(zeta.start, times = P)
      } else if (length(zeta.start) != P) {
        stop(paste("'zeta.start' should be a vector of length ", P, ".", sep = ""))
      }
    }

    if (length(delta.start) == 1) {
      if (hierarchical) {
        delta.start <- rep(delta.start, times = G + M - 1)
      } else {
        delta.start <- rep(delta.start, times = M)
      }
    } else {
      if (hierarchical & length(delta.start) != (G + M - 1)) {
        stop(paste("Incorrect input for delta.start. Its length should be ", G + M - 1,
                   sep = ""))
      } else if (!hierarchical & length(delta.start) != M) {
        stop(paste("Incorrect input for delta.start. Its length should be ", M,
                   sep = ""))
      }
    }
  }

###
### check the dimension of prior parameters
###

  ## prior mean of beta
  if (length(mu.beta) == 1) {
    mu.beta <- matrix(mu.beta, nrow = J, ncol = 2)
  } else if (class(mu.beta) != "matrix" | dim(mu.beta)[1] != J |
             dim(mu.beta)[2] != 2) {
    stop(paste("'mu.beta' should be a ", J, "-by-", 2, " matrix.", sep = ""))
  }

  ## prior mean of theta  (ADD 'mu.lambda' in the next update)
  if (length(mu.theta) == 1) {
    mu.theta <- rep(mu.theta, times = (M - 1 + P) * K)
  } else if (identical.lambda & hierarchical & (class(mu.theta) != "matrix" |
                                                dim(mu.theta)[1] != (M - 1) |
                                                dim(mu.theta)[2] != K)) {
    stop(paste("'mu.theta' (prior mean of lambda) has incorrect dimensions. It should be a ",
               M - 1, "-by-", K, " matrix.",
               sep = ""))
  } else if (!identical.lambda & hierarchical & (class(mu.theta) != "matrix" |
                                                 dim(mu.theta)[1] != M - 1 + P |
                                                 dim(mu.theta)[2] != K)) {
    stop(paste("'mu.theta' (prior mean of theta) has incorrect dimensions. It should be a ",
               M - 1 + P, "-by-", K, " matrix.",
               sep = ""))
  } else if (!identical.lambda & !hierarchical & (class(mu.theta) != "matrix" |
                                                  dim(mu.theta)[1] != M |
                                                  dim(mu.theta)[2] != K)) {
    stop(paste("'mu.theta' (prior mean of theta) has incorrect dimensions. It should be a ",
               M, "-by-", K, " matrix.",
               sep = ""))
  }

  ## prior mean of delta
  if (covariates) {
    if (length(mu.delta) == 1) {
      if (hierarchical) {
        mu.delta <- rep(mu.delta, times = M - 1)
      } else {
        mu.delta <- rep(mu.delta, times = M)
      }
    } else if (length(mu.delta) != M - 1 + G){
      stop(paste("'mu.delta' must be a vector of length ", M - 1 + G, ".", sep = ""))
    }
  } else {
    if (length(mu.x) == 1) {
      mu.delta <- rep(mu.x, times = N)
    } else if (length(mu.x) == N) {
      mu.delta <- mu.x
    } else if (length(mu.x) != N) {
      stop(paste("'mu.x' must be a vector of length ", N, ".", sep = ""))
    }
  }
      

  if (hierarchical) {
    ## prior mean of kappa
    if (length(mu.kappa) == 1) {
      mu.kappa <- rep(mu.kappa, times = K * P)
    } else if (identical.lambda & (class(mu.kappa) != "matrix" | dim(mu.kappa)[1] != P |
                                   dim(mu.kappa)[2] != K)){
      stop(paste("Incorrect input for mu.kappa. It should be a ", P, "-by-", K, " matrix.",
                 sep = ""))
    }

    ## prior mean of zeta
    if (length(mu.zeta) == 1) {
      mu.zeta <- rep(mu.zeta, times = P)
    } else if (length(mu.zeta) != P){
      stop(paste("Incorrect input for mu.zeta. Its length should be ", P,
                 sep = ""))
    }
  }

  ## prior precision of beta
  if (length(precision.beta) == 1) {
    precision.beta <- diag(precision.beta, nrow = 2, ncol = 2)
  } else if (class(precision.beta) != "matrix" | dim(precision.beta)[1] != 2 |
             dim(precision.beta)[2] != 2) {
    stop("'precision.beta' must be a 2-by-2 matrix.")
  }

  ## prior precision of delta
  if (hierarchical) {
    if (length(precision.delta) == 1) {
      precision.delta <- diag(precision.delta, nrow = M - 1, ncol = M - 1)
    } else if (class(precision.delta) != "matrix" | dim(precision.delta)[1] != M - 1 |
               dim(precision.delta)[2] != M - 1) {
      stop(paste("'precision.delta' must be a ", M - 1, "-by", M - 1, " matrix.", sep = ""))
    }
  } else {
    if (length(precision.delta) == 1) {
      precision.delta <- diag(precision.delta, nrow = M, ncol = M)
    } else if (class(precision.delta) != "matrix" | dim(precision.delta)[1] != M |
               dim(precision.delta)[2] != M) {
      stop(paste("'precision.delta' must be a ", M, "-by", M, " matrix.", sep = ""))
    }
  }

  ## prior precision of theta (lambda)
  if (length(precision.theta) == 1) {
    precision.theta <- rep(precision.theta, times = M + P)
  } else if (hierarchical & identical.lambda & length(precision.theta) != M - 1) {
    stop(paste("'precision.theta' (prior precision of lambda) must be a vector of length ",
               M - 1, ".", sep = ""))
  } else if (hierarchical & !identical.lambda & length(precision.theta) != M - 1 + P) {
    stop(paste("'precision.theta' must be a vector of length ", M - 1 + P, ".", sep = ""))
  } else if (!hierarchical & identical.lambda & length(precision.theta) != M) {
    stop(paste("'precision.theta' (prior precision of lambda) must be a vector of length ",
               M, ".", sep = ""))
  } else if (!hierarchical & !identical.lambda & length(precision.theta) != M) {
    stop(paste("'precision.theta' must be a vector of length ", M, ".", sep = ""))
  }


  if (hierarchical) {
    ## prior precision of kappa
    if (length(precision.kappa) == 1) {
      precision.kappa <- diag(precision.kappa, nrow = P, ncol = P)
    } else if (hierarchical & identical.lambda & (class(precision.kappa) != "matrix" |
                                                  dim(precision.kappa)[1] != P |
                                                  dim(precision.kappa)[2] != P)) {
      stop(paste("'precision.kappa' must be a ", P, "-by-", P, " matrix.", sep = ""))
    }

    ## prior precision of zeta
    if (length(precision.zeta) == 1) {
      precision.zeta <- diag(precision.zeta, nrow = P, ncol = P)
    } else if (hierarchical & (class(precision.zeta) != "matrix" | dim(precision.zeta)[1] != P |
                               dim(precision.zeta)[2] != P)) {
      stop(paste("'precision.zeta' must be a ", P, "-by-", P, " matrix.", sep = ""))
    }
    
  }

  if (length(prop) != J) {
    prop <- rep(prop, times = J)
  }

  printout <- floor( (MCMC - burn) / thin )

  temp <- .C("R2endorse",
             as.integer(response),
             as.integer(endorse),
             as.double(cov.mat),
             as.double(cov.village.mat),
             as.integer(village),
             as.integer(c(N, J, M, K, max.L)),
             as.integer(L),
             as.integer(c(G, P)),
             as.double(x.start),
             as.double(s.start),
             as.double(beta.start),
             as.double(tau.start),
             as.double(lambda.start),
             as.double(omega2.start),
             as.double(theta.start),
             as.double(phi2.start),
             as.double(kappa.start),
             as.double(psi2.start),
             as.double(delta.start),
             as.double(zeta.start),
             as.double(rho2.start),
             as.double(mu.beta),
             as.double(precision.beta),
             as.double(precision.x),
             as.double(mu.theta),
             as.double(precision.theta),
             as.double(mu.kappa),
             as.double(precision.kappa),
             as.double(mu.delta),
             as.double(precision.delta),
             as.double(mu.zeta),
             as.double(precision.zeta),
             as.double(c(s0.omega2, nu0.omega2, s0.phi2, nu0.phi2, s0.psi2, nu0.psi2, s0.sig2,
                         nu0.sig2, s0.rho2, nu0.rho2)),
             as.integer(c(mda, mh, x.sd, tau.out, s.out, omega2.out, phi2.out,
                          psi2.out, verbose, seed.store, covariates,
                          identical.lambda, hierarchical, MCMC, burn, thin)),
             as.double(prop),
             betaStore = double(printout * 2 * J),
             tauStore = if (tau.out) double(printout * (max.L - 1) * J) else double(1),
             xStore = if (x.sd) double(printout) else double(printout * N),
             sStore = if (s.out) double(printout * N * J) else double(1),
             lambdaStore = if (!hierarchical & !identical.lambda) double(printout * J * K * M) else if (hierarchical & identical.lambda) double(printout * K * (G + M - 1)) else if (!hierarchical & identical.lambda) double(printout * K * M) else if (hierarchical & !identical.lambda) double(printout * J * K * (G + M - 1)),
             thetaStore = if (hierarchical & !identical.lambda) double(printout * K * (M - 1 + P)) else if (!hierarchical & !identical.lambda) double(printout * K * M) else double(1),
             kappaStore = if (hierarchical & !identical.lambda) double(printout * J * K * P) else if (hierarchical & identical.lambda) double(printout * K * P) else double(1),
             deltaStore = if (hierarchical) double(printout * (G + M - 1)) else if (covariates) double(printout * M) else double(1),
             zetaStore = if (hierarchical) double(printout * P) else double(1),
             omega2Store = if (omega2.out & !identical.lambda) double(printout * J * K) else if (omega2.out) double(printout * K) else double(1),
             phi2Store = if (phi2.out & !identical.lambda & hierarchical) double(printout * K * (M - 1 + P)) else if (phi2.out & !identical.lambda & !hierarchical) double(printout * K * M) else double(1),
             psi2Store = if (psi2.out & hierarchical & identical.lambda) double(printout * K) else if (psi2.out & hierarchical & !identical.lambda) double(printout * J * K) else double(1),
	     sig2Store = if (covariates | hierarchical) double(printout) else double(1),
             rho2Store = if (hierarchical) double(printout) else double(1), 
             betaLast = if (seed.store) double(2 * J) else double(1),
             tauLast = if (seed.store) double((max.L - 1) * J) else double(1),
             xLast = if (seed.store) double(N) else double(1),
             sLast = if (seed.store) double(N * J) else double(1),
             lambdaLast = if (seed.store & !identical.lambda & !hierarchical) double(J * K * M) else if (seed.store & !identical.lambda & hierarchical) double(K * J * (M + G - 1)) else if (seed.store & identical.lambda & !hierarchical) double(K * M) else if (seed.store & identical.lambda & hierarchical) double(K * (M + G - 1)) else double(1),
             thetaLast = if (seed.store & !identical.lambda & !hierarchical) double(K * M) else if (seed.store & !identical.lambda & hierarchical) double(K * (M - 1 + P)) else double(1),
             kappaLast = if (seed.store & identical.lambda & hierarchical) double(K * P) else if (seed.store & !identical.lambda & hierarchical) double(J * K * P) else double(1),
             deltaLast = if (seed.store & hierarchical) double(G + M - 1) else if (seed.store & !hierarchical) double (M) else double(1),
             zetaLast = if (seed.store & hierarchical) double(P) else double(1),
             omega2Last = if (seed.store & !identical.lambda) double(J * K) else if (seed.store) double(K) else double(1),
             phi2Last = if (seed.store) double(K * M) else double(1),
             psi2Last = if (seed.store & identical.lambda & hierarchical) double(K) else if (seed.store & !identical.lambda & hierarchical) double(J * K) else double(1),
	     sig2Last = double(1),
             rho2Last = double(1),
             accept.ratio = double(J),
             package = "endorse")

  seedStore <- .Random.seed

  res <- list(beta = matrix(as.double(temp$betaStore), byrow = TRUE, ncol = 2 * J,
                nrow = printout),
              tau = if (tau.out) matrix(as.double(temp$tauStore), byrow = TRUE,
                ncol = (max.L - 1) * J, nrow = printout) else NULL,
              x = if (x.sd) matrix(as.double(temp$xStore)[1:printout], ncol = 1,
                nrow = printout) else matrix(as.double(temp$xStore), byrow = TRUE,
                  ncol = N, nrow = printout),
              s = if (s.out) matrix(as.double(temp$sStore), byrow = TRUE,
                ncol = N * J, nrow = printout) else NULL,
              lambda = if (identical.lambda & hierarchical) matrix(as.double(temp$lambdaStore), byrow = TRUE, ncol = (G + M - 1) * K, nrow = printout) else if (identical.lambda & !hierarchical) matrix(as.double(temp$lambdaStore), byrow = TRUE, ncol = M * K, nrow = printout) else if (!identical.lambda & hierarchical) matrix(as.double(temp$lambdaStore), byrow = TRUE, ncol = (G + M - 1) * J * K, nrow = printout) else if (!identical.lambda & !hierarchical) matrix(as.double(temp$lambdaStore), byrow = TRUE, ncol = J * M * K, nrow = printout),
              theta = if (identical.lambda) NULL else if (hierarchical) matrix(as.double(temp$thetaStore), byrow = TRUE, ncol = K * (M - 1 + P), nrow = printout) else matrix(as.double(temp$thetaStore), byrow = TRUE, ncol = K * M, nrow = printout),
              kappa = if (hierarchical & identical.lambda) matrix(as.double(temp$kappaStore), byrow = TRUE, ncol = P * K, nrow = printout) else if (hierarchical) matrix(as.double(temp$kappaStore), byrow = TRUE, ncol = P * K * J, nrow = printout) else NULL,
              delta = if (hierarchical) matrix(as.double(temp$deltaStore), byrow = TRUE, ncol = G + M - 1, nrow = printout) else if (covariates) matrix(as.double(temp$deltaStore), byrow = TRUE, ncol = M, nrow = printout) else NULL,
              zeta = if (hierarchical) matrix(as.double(temp$zetaStore), byrow = TRUE, ncol = P, nrow = printout) else NULL,
              omega2 = if (identical.lambda & omega2.out) matrix(as.double(temp$omega2Store), byrow = TRUE, ncol = K, nrow = printout) else if (omega2.out) matrix(as.double(temp$omega2Store), byrow = TRUE, ncol = J * K, nrow = printout) else NULL,
              phi2 = if (!identical.lambda & hierarchical & phi2.out) matrix(as.double(temp$phi2Store), byrow = TRUE, ncol = K * (M - 1 + P), nrow = printout) else if (!identical.lambda & !hierarchical & phi2.out) matrix(as.double(temp$phi2Store), byrow = TRUE, ncol = K * M, nrow = printout) else NULL,
              psi2 = if (hierarchical & identical.lambda & psi2.out) matrix(as.double(temp$psi2Store), byrow = TRUE, ncol = K, nrow = printout) else if (hierarchical & !identical.lambda & psi2.out) matrix(as.double(temp$psi2Store), byrow = TRUE, ncol = K * J, nrow = printout) else NULL,
	      sig2 = if (hierarchical | covariates) matrix(as.double(temp$sig2Store), ncol = 1, nrow = printout) else NULL,
              rho2 = if (hierarchical) matrix(as.double(temp$rho2Store), byrow = TRUE, ncol = 1, nrow = printout) else NULL,
              accept.ratio = if (mh) as.double(temp$accept.ratio) else NULL,
              seed = if (seed.store) seedStore else NULL,
              beta.start = if (seed.store) matrix(as.double(temp$betaLast), nrow = J, ncol = 2, byrow = TRUE) else NULL,
              tau.start = if (seed.store) as.double(temp$tauLast) else NULL,
              x.start = if (seed.store) as.double(temp$xLast) else NULL,
              s.start = if (seed.store) matrix(as.double(temp$sLast), nrow = N, ncol = J, byrow = TRUE) else NULL,
              lambda.start = if (seed.store) as.double(temp$lambdaLast) else NULL,
              theta.start = if (seed.store & !identical.lambda) as.double(temp$thetaLast) else NULL,
              kappa.start = if (seed.store & hierarchical) as.double(temp$kappaLast) else NULL,
              delta.start = if (seed.store & (covariates | hierarchical)) as.double(temp$deltaLast) else NULL,
              zeta.start = if (seed.store & hierarchical) as.double(temp$zetaLast) else NULL,
              omega2.start = if (seed.store) as.double(temp$omega2Last) else NULL,
              phi2.start = if (seed.store & !identical.lambda) as.double(temp$phi2Last) else NULL,
              psi2.start = if (seed.store & hierarchical) as.double(temp$psi2Last) else NULL,
	      sig2.start = if (seed.store & (covariates | hierarchical)) as.double(temp$sig2.start) else NULL,
              rho2.start = if (seed.store & hierarchical) as.double(temp$rho2Last) else NULL,
              village.indicator = village + 1,
              model.matrix.indiv = cov.mat,
              formula.indiv = formula.indiv,
              model.matrix.village = if (hierarchical) cov.village.mat else NULL,
              formula.village = if (hierarchical) formula.village else NULL,
              hierarchical = hierarchical,
              identical.lambda = identical.lambda,
              num.act = K,
              num.policy = J,
              treat = endorse)

  colnames(res$beta) <- paste(rep(c("alpha", "beta"), times = J),
                              rep(1:J, each = 2), sep = ".")
  res$beta <- mcmc(res$beta, start = burn + 1, end = MCMC, thin = thin)

  if (tau.out) {
    temp.names <- paste("tau", 1:J, sep = "")
    colnames(res$tau) <- paste(rep(temp.names, each = (max.L - 1)),
                               rep(1:(max.L - 1), times = J), sep = ".")
    res$tau <- mcmc(res$tau, start = burn + 1, end = MCMC, thin = thin)
  }

  if (x.sd) {
    colnames(res$x) <- "sd.x"
    res$x <- mcmc(res$x, start = burn + 1, end = MCMC, thin = thin)
  } else {
    colnames(res$x) <- paste("x", 1:N, sep = ".")
    res$x <- mcmc(res$x, start = burn + 1, end = MCMC, thin = thin)
  }

  if (s.out) {
    colnames(res$s) <- paste("s", rep(1:nrow(data), each = J),
                             rep(1:J, times = nrow(data)), sep = ".")
  }

  
  if (identical.lambda) {
    if (hierarchical) {
      if (covariates) {
        colnames(res$lambda) <- paste(rep(c(paste("group", 1:G, sep = ""), var.names.indiv[2:M]), times = K),
                                      rep(1:K, each = G + M - 1), sep = ".")
      } else {
        colnames(res$lambda) <- paste(rep(paste("group", 1:G, sep = ""), times = K),
                                      rep(1:K, each = G), sep = ".")
      }
    } else {
      colnames(res$lambda) <- paste(rep(var.names.indiv, times = K),
                                    rep(1:K, each = M), sep =".")
    }
  } else {
    if (hierarchical) {
      temp.village.label <- paste("group", 1:G, sep = "")
      if (covariates) {
        temp.col.names <- paste(c(temp.village.label, var.names.indiv[2:M]),
                                rep(1:J, each = (M + G - 1) * K), sep = ".")
        colnames(res$lambda) <- paste(temp.col.names, rep(rep(1:K, each = G + M - 1), times = J),
                                      sep = ".")
      } else {
        temp.col.names <- paste(temp.village.label, rep(1:J, each = G * K), sep = ".")
        colnames(res$lambda) <- paste(temp.col.names, rep(rep(1:K, each = G), times = J),
                                      sep = ".")
      }
    } else {
      temp.col.names <- paste(var.names.indiv, rep(1:J, each = M * K), sep = ".")
      colnames(res$lambda) <- paste(temp.col.names, rep(rep(1:K, each = M), times = J), sep = ".")
    }
  }

  res$lambda <- mcmc(res$lambda, start = burn + 1, end = MCMC, thin = thin)


  if (identical.lambda) {
    colnames(res$omega2) <- paste("omega2", 1:K, sep = ".")
  } else {
    colnames(res$omega2) <- paste("omega2", rep(1:J, each = K), rep(1:K,
                                                       times = J), sep = ".")
  }
  res$omega2 <- mcmc(res$omega2, start = burn + 1, end = MCMC, thin = thin)  


  if (!identical.lambda) {
    if (hierarchical) {
      temp.col.names <- c(rep("indiv", times = (M - 1)), rep("group", times = P))
      if (M > 1) {
        temp.col.names <- paste(temp.col.names, c(var.names.indiv[2:M], var.names.village),
                                sep = ".")
        temp.col.names <- paste(temp.col.names, rep(1:K, each = (M - 1) + P), sep = ".")
      } else {
        temp.col.names <- paste(temp.col.names, rep(var.names.village, times = K),
                                rep(1:K, each = P), sep = ".")
      }
      colnames(res$theta) <- temp.col.names
    } else {
      if (covariates) {
        temp.col.names <- paste(var.names.indiv, rep(1:K, each = M), sep = ".")
        colnames(res$theta) <- temp.col.names
      } else {
        temp.col.names <- paste("theta", 1:K, sep = "")
        colnames(res$theta) <- temp.col.names
      }

      temp.col.names <- paste("phi2", temp.col.names, sep = ".")
      colnames(res$phi2) <- temp.col.names
    }
    res$theta <- mcmc(res$theta, start = burn + 1, end = MCMC, thin = thin)
    res$phi2 <- mcmc(res$phi2, start = burn + 1, end = MCMC, thin = thin)
  }

  if (hierarchical) {
    if (identical.lambda) {
      colnames(res$kappa) <- paste(rep(var.names.village, times = K),
                                   rep(1:K, each = P), sep = ".")
      res$kappa <- mcmc(res$kappa, start = burn + 1, end = MCMC, thin = thin)
      if (psi2.out) {
        colnames(res$psi2) <- paste("psi2", 1:K, sep = ".")
        res$psi2 <- mcmc(res$psi2, start = burn + 1, end = MCMC, thin = thin)
      }
    } else {
      colnames(res$kappa) <- paste(rep(var.names.village, times = J * K),
                                   rep(rep(1:J, each = P), times = K),
                                   rep(1:K, each = J * P), sep = ".")
      res$kappa <- mcmc(res$kappa, start = burn + 1, end = MCMC, thin = thin)
      if (psi2.out) {
        colnames(res$psi2) <- paste("psi2",
                                    rep(1:J, times = K),
                                    rep(1:K, each = J), sep = ".")
        res$psi2 <- mcmc(res$psi2, start = burn + 1, end = MCMC, thin = thin)
      }
    }
  }

  if (hierarchical) {
    if (covariates) {
      colnames(res$delta) <- c(paste("village", 1:G, sep = "."), var.names.indiv[2:M])
    } else {
      colnames(res$delta) <- paste("village", 1:G, sep = ".")
    }
    res$delta <- mcmc(res$delta, start = burn + 1, end = MCMC, thin = thin)  
  } else if (covariates) {
    colnames(res$delta) <- var.names.indiv
    res$delta <- mcmc(res$delta, start = burn + 1, end = MCMC, thin = thin)  
  }

  if (hierarchical) {
    colnames(res$sig2) <- "sig2"
    res$sig2 <- mcmc(res$sig2, start = burn + 1, end = MCMC, thin = thin)
  }

  if (hierarchical) {
    colnames(res$zeta) <- var.names.village
    res$zeta <- mcmc(res$zeta, start = burn + 1, end = MCMC, thin = thin)
    colnames(res$rho2) <- "rho2"
    res$rho2 <- mcmc(res$rho2, start = burn + 1, end = MCMC, thin = thin)
  }

  names(res$accept.ratio) <- paste("Q", 1:J, sep = "")

  return(res)
}  
