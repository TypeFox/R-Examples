
coverage <- function(gbp.object, A.or.r, reg.coef, mean.PriorDist, nsim = 100) {

  # Rao-Blackwellized criterion
  coverageRB <- matrix(NA, nrow = length(gbp.object$se), ncol = nsim)

  # 1-0 criterion that is 1 if interval includes true parameter, 0 if not
  coverageS <- matrix(NA, nrow = length(gbp.object$se), ncol = nsim)

  if (missing(A.or.r) & missing(reg.coef) & missing(mean.PriorDist)) {
    only.gbp.result <- TRUE
  } else {
    only.gbp.result <- FALSE
  }

  betas <- NA

  if (gbp.object$model == "pr" & is.na(gbp.object$prior.mean)) {
    print("Model is Poisson and the prior mean is unknown. We currently do not provide Frequency Method Checking for this model.")
    stop()
  }
######
  if (gbp.object$model == "br") {
    if (sum(is.na(gbp.object$weight)) == 1 & sum(is.na(gbp.object$p)) == 1) {
      AR <- FALSE
    } else {
      AR <- TRUE
    }
  }

  # if model=BRIMM	
  if (gbp.object$model == "br") {
    if (only.gbp.result) {

      # 1. initial values
      if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
        temp.x <- as.matrix(rep(1, length(gbp.object$se)))
        betas <- as.vector(gbp.object$beta.new)
        p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
      } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
        temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
        betas <- as.vector(gbp.object$beta.new)
        p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
        X <- gbp.object$X
      } else if (!is.na(gbp.object$prior.mean)) {
        p0 <- gbp.object$prior.mean
      }
  
      n <- gbp.object$se
      r <- as.numeric(exp(-gbp.object$a.new))
      priormeanused <- p0

      # 2. generate p matrix
      sim.p <- matrix(rbeta(length(n) * nsim, r * p0, r * (1 - p0)), nrow = length(n))

      # 3. generate z (data) matrix
      sim.z <- matrix(rbinom(nrow(sim.p) * nsim, n, sim.p), nrow = length(n))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
          if (AR == 0) {
            out <- if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
                     gbp(sim.z[, i], n, model = "binomial", confidence.lvl = gbp.object$confidence.lvl)
                   } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
                     gbp(sim.z[, i], n, X, model = "binomial", confidence.lvl = gbp.object$confidence.lvl)
                   } else if (!is.na(gbp.object$prior.mean)) {
                     gbp(sim.z[, i], n, mean.PriorDist = p0, model = "binomial", 
                         confidence.lvl = gbp.object$confidence.lvl)
                   }
          } else {
            out <- if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
                     gbp(sim.z[, i], n, model = "binomial", confidence.lvl = gbp.object$confidence.lvl,
                         n.AR = gbp.object$n.AR, n.AR.factor = gbp.object$n.AR.factor,
                         trial.scale = gbp.object$trial.scale, save.result = FALSE, 
                         t = gbp.object$c, u = gbp.object$u)
                   } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
                     gbp(sim.z[, i], n, X, model = "binomial", confidence.lvl = gbp.object$confidence.lvl,
                         n.AR = gbp.object$n.AR, n.AR.factor = gbp.object$n.AR.factor, 
                         trial.scale = gbp.object$trial.scale, save.result = FALSE,
                         t = gbp.object$c, u = gbp.object$u)
                   } else if (!is.na(gbp.object$prior.mean)) {
                     gbp(sim.z[, i], n, mean.PriorDist = p0, model = "binomial", 
                         confidence.lvl = gbp.object$confidence.lvl, n.AR = gbp.object$n.AR, 
                         n.AR.factor = gbp.object$n.AR.factor,
                         trial.scale = gbp.object$trial.scale, save.result = FALSE,
                         t = gbp.object$c, u = gbp.object$u)
                   }
          }

          a1 <- r * p0 + sim.z[, i]
          a0 <- r * (1 - p0) + n - sim.z[, i]
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pbeta(upp, a1, a0) - pbeta(low, a1, a0)
          coverageS[, i] <- ifelse(low <= sim.p[, i] & sim.p[, i] <= upp, 1, 0)

        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	

    } else if (!only.gbp.result) {

      # 1. initial values
      if (missing(A.or.r) & missing(reg.coef) & !missing(mean.PriorDist)) {
        p0 <- mean.PriorDist
        r <- as.numeric(exp(-gbp.object$a.new))
      } else if (missing(A.or.r) & !missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          print("reg.coef cannot be designated because there is no covariate.")
          stop()
        } else if (identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(rep(1, length(gbp.object$se)))
          betas <- as.vector(reg.coef)
        } else {
          temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
          betas <- as.vector(reg.coef)
        }
        p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
        r <- as.numeric(exp(-gbp.object$a.new))
      } else if (!missing(A.or.r) & missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          p0 <- gbp.object$prior.mean
        } else if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(rep(1, length(gbp.object$se)))
          betas <- as.vector(gbp.object$beta.new)
          p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
        } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
          betas <- as.vector(gbp.object$beta.new)
          p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
        }
        r <- A.or.r
      } else if (missing(A.or.r) & !missing(reg.coef) & !missing(mean.PriorDist)) {
        print("reg.coef and mean.PriorDist cannot be designated at the same time because once we know mean.PriorDist we do not need to estimate reg.coef.")
        stop()
      } else if (!missing(A.or.r) & missing(reg.coef) & !missing(mean.PriorDist)) {
        p0 <- mean.PriorDist
        r <- A.or.r
      } else if (!missing(A.or.r) & !missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          print("reg.coef cannot be designated because second-level mean is known in the gbp object to begin with.")
          stop()
        } else if (identical(gbp.object$prior.mean, NA)) {
          if (identical(gbp.object$X, NA)) {
            temp.x <- as.matrix(rep(1, length(gbp.object$se)))
            betas <- as.vector(reg.coef)
          } else {
            temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
            betas <- as.vector(reg.coef)
          }
          p0 <- exp(temp.x %*% betas) / (1 + exp(temp.x %*% betas))
          r <- A.or.r
        }
      } else if (!missing(A.or.r) & !missing(reg.coef) & !missing(mean.PriorDist)) {
        print("reg.coef and mean.PriorDist cannot be designated at the same time because once we know mean.PriorDist we do not need to estimate reg.coef.")
        stop()
      }

      n <- gbp.object$se
      priormeanused <- p0
 
      # 2. generate p matrix
      sim.p <- matrix(rbeta(length(n) * nsim, r * p0, r * (1 - p0)), nrow = length(n))

      # 3. generate z (data) matrix
      sim.z <- matrix(rbinom(nrow(sim.p) * nsim, n, sim.p), nrow = length(n))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({

          if (AR == 0) {
            out <- if (!missing(mean.PriorDist)) {
                     gbp(sim.z[, i], n, mean.PriorDist = mean.PriorDist, model = "binomial", 
                         confidence.lvl = gbp.object$confidence.lvl)
                   } else if (!missing(A.or.r) & missing(reg.coef)) { 
                     if (!identical(gbp.object$prior.mean, NA)) {
                       gbp(sim.z[, i], n, mean.PriorDist = gbp.object$prior.mean, model = "binomial", 
                           confidence.lvl = gbp.object$confidence.lvl)
                     } else if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, model = "binomial", confidence.lvl = gbp.object$confidence.lvl)
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, gbp.object$X, model = "binomial", confidence.lvl = gbp.object$confidence.lvl)
                     }
                   } else if (!missing(reg.coef)) {
                     if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, model = "binomial", confidence.lvl = gbp.object$confidence.lvl)
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, gbp.object$X, model = "binomial", confidence.lvl = gbp.object$confidence.lvl)
                     }
                   }
          
          } else {
            out <- if (!missing(mean.PriorDist)) {
                     gbp(sim.z[, i], n, mean.PriorDist = mean.PriorDist, model = "binomial", 
                         confidence.lvl = gbp.object$confidence.lvl, n.AR = gbp.object$n.AR,
                         n.AR.factor = gbp.object$n.AR.factor,
                         trial.scale = gbp.object$trial.scale, save.result = FALSE,
                         t = gbp.object$c, u = gbp.object$u)
                   } else if (!missing(A.or.r) & missing(reg.coef)) { 
                     if (!identical(gbp.object$prior.mean, NA)) {
                       gbp(sim.z[, i], n, mean.PriorDist = gbp.object$prior.mean, model = "binomial", 
                           confidence.lvl = gbp.object$confidence.lvl, n.AR = gbp.object$n.AR,
                           n.AR.factor = gbp.object$n.AR.factor,
                           trial.scale = gbp.object$trial.scale, save.result = FALSE,
                           t = gbp.object$c, u = gbp.object$u)
                     } else if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, model = "binomial", confidence.lvl = gbp.object$confidence.lvl,
                           n.AR = gbp.object$n.AR, trial.scale = gbp.object$trial.scale, 
                           n.AR.factor = gbp.object$n.AR.factor, save.result = FALSE,
                           t = gbp.object$c, u = gbp.object$u)
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, gbp.object$X, model = "binomial", confidence.lvl = gbp.object$confidence.lvl,
                           n.AR = gbp.object$n.AR, trial.scale = gbp.object$trial.scale, 
                           n.AR.factor = gbp.object$n.AR.factor, save.result = FALSE,
                           t = gbp.object$c, u = gbp.object$u)
                     }
                   } else if (!missing(reg.coef)) {
                     if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, model = "binomial", confidence.lvl = gbp.object$confidence.lvl,
                           n.AR = gbp.object$n.AR, trial.scale = gbp.object$trial.scale, 
                           n.AR.factor = gbp.object$n.AR.factor, save.result = FALSE,
                           t = gbp.object$c, u = gbp.object$u)
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, gbp.object$X, model = "binomial", confidence.lvl = gbp.object$confidence.lvl,
                           n.AR = gbp.object$n.AR, trial.scale = gbp.object$trial.scale, 
                           n.AR.factor = gbp.object$n.AR.factor, save.result = FALSE,
                           t = gbp.object$c, u = gbp.object$u)
                     }
                   }

          }

          a1 <- r * p0 + sim.z[, i]
          a0 <- r * (1 - p0) + n - sim.z[, i]
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pbeta(upp, a1, a0) - pbeta(low, a1, a0)
          coverageS[, i] <- ifelse(low <= sim.p[, i] & sim.p[, i] <= upp, 1, 0)

        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	
    }

  } else if(gbp.object$model == "pr") {

    if (only.gbp.result) {

      # 1. initial values
      if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
        temp.x <- as.matrix(rep(1, length(gbp.object$se)))
        betas <- as.vector(gbp.object$beta.new)
        lambda0 <- exp(temp.x %*% betas)
      } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
        temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
        betas <- as.vector(gbp.object$beta.new)
        lambda0 <- exp(temp.x %*% betas)
        X <- gbp.object$X
      } else if (!is.na(gbp.object$prior.mean)) {
        lambda0 <- gbp.object$prior.mean
      }
  
      n <- gbp.object$se
      r <- exp(-gbp.object$a.new)
      priormeanused <- lambda0

      # 2. generate lambda matrix
      sim.lambda <- matrix(rgamma(length(n) * nsim, r * lambda0, r), nrow = length(n))

      # 3. generate z (data) matrix
      sim.z <- matrix(rpois(nrow(sim.lambda) * nsim, n * sim.lambda), nrow = length(n))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
          out <- if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
                   gbp(sim.z[, i], n, model = "poisson", confidence.lvl = gbp.object$confidence.lvl)
                 } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
                   gbp(sim.z[, i], n, X, model = "poisson", confidence.lvl = gbp.object$confidence.lvl)
                 } else if (!is.na(gbp.object$prior.mean)) {
                   gbp(sim.z[, i], n, mean.PriorDist = lambda0, model = "poisson", confidence.lvl = gbp.object$confidence.lvl)
                 }
          
          sh <- r * lambda0 + sim.z[, i]
          rt <- n + r
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pgamma(upp, sh, rt) - pgamma(low, sh, rt)
          coverageS[, i] <- ifelse(low <= sim.lambda[, i] & sim.lambda[, i] <= upp, 1, 0)

        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	

    } else if (!only.gbp.result) {

      # 1. initial values
      if (missing(A.or.r) & missing(reg.coef) & !missing(mean.PriorDist)) {
        lambda0 <- mean.PriorDist
        r <- exp(-gbp.object$a.new)
      } else if (missing(A.or.r) & !missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          print("reg.coef cannot be designated because there is no covariate.")
          stop()
        } else if (identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(rep(1, length(gbp.object$se)))
          betas <- as.vector(reg.coef)
        } else {
          temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
          betas <- as.vector(reg.coef)
        }
        lambda0 <- exp(temp.x %*% betas)
        r <- exp(-gbp.object$a.new)
      } else if (!missing(A.or.r) & missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          lambda0 <- gbp.object$prior.mean
        } else if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(rep(1, length(gbp.object$se)))
          betas <- as.vector(gbp.object$beta.new)
          lambda0 <- exp(temp.x %*% betas)
        } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
          betas <- as.vector(gbp.object$beta.new)
          lambda0 <- exp(temp.x %*% betas)
        }
        r <- A.or.r
      } else if (missing(A.or.r) & !missing(reg.coef) & !missing(mean.PriorDist)) {
        print("reg.coef and mean.PriorDist cannot be designated at the same time because once we know mean.PriorDist we do not need to estimate reg.coef.")
        stop()
      } else if (!missing(A.or.r) & missing(reg.coef) & !missing(mean.PriorDist)) {
        lambda0 <- mean.PriorDist
        r <- A.or.r
      } else if (!missing(A.or.r) & !missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          print("reg.coef cannot be designated because second-level mean is known in the gbp object to begin with.")
          stop()
        } else if (identical(gbp.object$prior.mean, NA)) {
          if (identical(gbp.object$X, NA)) {
            temp.x <- as.matrix(rep(1, length(gbp.object$se)))
            betas <- as.vector(reg.coef)
          } else {
            temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
            betas <- as.vector(reg.coef)
          }
          lambda0 <- exp(temp.x %*% betas)
          r <- A.or.r
        }
      } else if (!missing(A.or.r) & !missing(reg.coef) & !missing(mean.PriorDist)) {
        print("reg.coef and mean.PriorDist cannot be designated at the same time because once we know mean.PriorDist we do not need to estimate reg.coef.")
        stop()
      }

      n <- gbp.object$se
      priormeanused <- lambda0
 
      # 2. generate lambda matrix
      sim.lambda <- matrix(rgamma(length(n) * nsim, r * lambda0, r), nrow = length(n))

      # 3. generate z (data) matrix
      sim.z <- matrix(rpois(nrow(sim.lambda) * nsim, n * sim.lambda), nrow = length(n))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
            out <- if (!missing(mean.PriorDist)) {
                     gbp(sim.z[, i], n, mean.PriorDist = mean.PriorDist, model = "poisson", 
                         confidence.lvl = gbp.object$confidence.lvl)
                   } else if (!missing(A.or.r) & missing(reg.coef)) { 
                     if (!identical(gbp.object$prior.mean, NA)) {
                       gbp(sim.z[, i], n, mean.PriorDist = gbp.object$prior.mean, model = "poisson", 
                           confidence.lvl = gbp.object$confidence.lvl)
                     } else if (identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, model = "poisson", confidence.lvl = gbp.object$confidence.lvl)
                     } else if (!identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, gbp.object$X, model = "poisson", confidence.lvl = gbp.object$confidence.lvl)
                     }
                   } else if (!missing(reg.coef)) {
                     if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, model = "poisson", confidence.lvl = gbp.object$confidence.lvl)
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.z[, i], n, gbp.object$X, model = "poisson", confidence.lvl = gbp.object$confidence.lvl)
                     }
                   }
          
          sh <- r * lambda0 + sim.z[, i]
          rt <- n + r
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pgamma(upp, sh, rt) - pgamma(low, sh, rt)
          coverageS[, i] <- ifelse(low <= sim.lambda[, i] & sim.lambda[, i] <= upp, 1, 0)

        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	
    }

  } else if (gbp.object$model == "gr") {	

     if (only.gbp.result) {
      # 1. initial values
      if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
        temp.x <- as.matrix(rep(1, length(gbp.object$se)))
        betas <- as.vector(gbp.object$beta.new)
        mu0 <- temp.x %*% betas
      } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
        temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
        betas <- as.vector(gbp.object$beta.new)
        mu0 <- temp.x %*% betas
        X <- gbp.object$X
      } else if (!is.na(gbp.object$prior.mean)) {
        mu0 <- gbp.object$prior.mean
      }
  
      se <- gbp.object$se
      A <- exp(gbp.object$a.new)
      priormeanused <- mu0

      # 2. generate mu matrix
      sim.mu <- matrix(rnorm(length(se) * nsim, mu0, sqrt(A)), nrow = length(se))

      # 3. generate y (data) matrix
      sim.y <- matrix(rnorm(nrow(sim.mu) * nsim, sim.mu, se), nrow = length(se))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
          out <- if (is.na(gbp.object$prior.mean) & identical(gbp.object$X, NA)) {
                   gbp(sim.y[, i], se, confidence.lvl = gbp.object$confidence.lvl)
                 } else if (is.na(gbp.object$prior.mean) & !identical(gbp.object$X, NA)) {
                   gbp(sim.y[, i], se, X, confidence.lvl = gbp.object$confidence.lvl)
                 } else if (!is.na(gbp.object$prior.mean)) {
                   gbp(sim.y[, i], se, mean.PriorDist = mu0, confidence.lvl = gbp.object$confidence.lvl)
                 }
          postmean <- mu0 * (se^2 / (se^2 + A)) + sim.y[, i] * (A / (se^2 + A))
          postsd <- sqrt(se^2 * (A / (se^2 + A)))
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pnorm(upp, postmean, postsd) - pnorm(low, postmean, postsd)
          coverageS[, i] <- ifelse(low <= sim.mu[, i] & sim.mu[, i] <= upp, 1, 0)

        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	

    } else if (!only.gbp.result) {

      # 1. initial values
      if (missing(A.or.r) & missing(reg.coef) & !missing(mean.PriorDist)) {
        mu0 <- mean.PriorDist
        A <- exp(gbp.object$a.new)
      } else if (missing(A.or.r) & !missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          print("reg.coef cannot be designated because there is no covariate.")
          stop()
        } else if (identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(rep(1, length(gbp.object$se)))
          betas <- as.vector(reg.coef)
        } else {
          temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
          betas <- as.vector(reg.coef)
        }
        mu0 <- temp.x %*% betas
        A <- exp(gbp.object$a.new)
      } else if (!missing(A.or.r) & missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          mu0 <- gbp.object$prior.mean
        } else if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(rep(1, length(gbp.object$se)))
          betas <- as.vector(gbp.object$beta.new)
          mu0 <- temp.x %*% betas
        } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)) {
          temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
          betas <- as.vector(gbp.object$beta.new)
          mu0 <- temp.x %*% betas
        }
        A <- A.or.r
      } else if (missing(A.or.r) & !missing(reg.coef) & !missing(mean.PriorDist)) {
        print("reg.coef and mean.PriorDist cannot be designated at the same time because once we know mean.PriorDist we do not need to estimate reg.coef.")
        stop()
      } else if (!missing(A.or.r) & missing(reg.coef) & !missing(mean.PriorDist)) {
        mu0 <- mean.PriorDist
        A <- A.or.r
      } else if (!missing(A.or.r) & !missing(reg.coef) & missing(mean.PriorDist)) {
        if (!identical(gbp.object$prior.mean, NA)) {
          print("reg.coef cannot be designated because second-level mean is known in the gbp object to begin with.")
          stop()
        } else if (identical(gbp.object$prior.mean, NA)) {
          if (identical(gbp.object$X, NA)) {
            temp.x <- as.matrix(rep(1, length(gbp.object$se)))
            betas <- as.vector(reg.coef)
          } else {
            temp.x <- as.matrix(cbind(rep(1, length(gbp.object$se)), gbp.object$X))
            betas <- as.vector(reg.coef)
          }
          mu0 <- temp.x %*% betas
          A <- A.or.r
        }
      } else if (!missing(A.or.r) & !missing(reg.coef) & !missing(mean.PriorDist)) {
        print("reg.coef and mean.PriorDist cannot be designated at the same time because once we know mean.PriorDist we do not need to estimate reg.coef.")
        stop()
      }

      se <- gbp.object$se
      priormeanused <- mu0
 
      # 2. generate mu matrix
      sim.mu <- matrix(rnorm(length(se) * nsim, mu0, sqrt(A)), nrow = length(se))

      # 3. generate y (data) matrix
      sim.y <- matrix(rnorm(nrow(sim.mu) * nsim, sim.mu, se), nrow = length(se))

      # 4. simulation
      for (i in 1 : nsim) {
        tryCatch({
            out <- if (!missing(mean.PriorDist)) {
                     gbp(sim.y[, i], se, mean.PriorDist = mean.PriorDist, model = "gaussian", 
                         confidence.lvl = gbp.object$confidence.lvl)
                   } else if (!missing(A.or.r) & missing(reg.coef)) { 
                     if (!identical(gbp.object$prior.mean, NA)) {
                       gbp(sim.y[, i], se, mean.PriorDist = gbp.object$prior.mean, model = "gaussian", 
                           confidence.lvl = gbp.object$confidence.lvl)
                     } else if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.y[, i], se, model = "gaussian", confidence.lvl = gbp.object$confidence.lvl)
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.y[, i], se, gbp.object$X, model = "gaussian", confidence.lvl = gbp.object$confidence.lvl)
                     }
                   } else if (!missing(reg.coef)) {
                     if (identical(gbp.object$prior.mean, NA) & identical(gbp.object$X, NA)){
                       gbp(sim.y[, i], se, model = "gaussian", confidence.lvl = gbp.object$confidence.lvl)
                     } else if (identical(gbp.object$prior.mean, NA) & !identical(gbp.object$X, NA)){
                       gbp(sim.y[, i], se, gbp.object$X, model = "gaussian", confidence.lvl = gbp.object$confidence.lvl)
                     }
                   }

          postmean <- mu0 * (se^2 / (se^2 + A)) + sim.y[, i] * (A / (se^2 + A))
          postsd <- sqrt(se^2 * (A / (se^2 + A)))
          low <- out$post.intv.low
          upp <- out$post.intv.upp
          coverageRB[, i] <- pnorm(upp, postmean, postsd) - pnorm(low, postmean, postsd)
          coverageS[, i] <- ifelse(low <= sim.mu[, i] & sim.mu[, i] <= upp, 1, 0)

        }, error = function(x) {
                     print(c(i,"error"))
                   }, warning = function(x) {
                                  print(c(i, "warning"))
                                })
      }	
    }
  }

  # average coverage probability
  result <- round(rowMeans(coverageRB, na.rm = TRUE), 3)
  avr.cov <- round(mean(result), 3)
  se.cov <- round(sqrt(apply(coverageRB, 1, var, na.rm = TRUE) / 
                            sum(!is.na(coverageRB[1, ]))), 4)
  avr.se.cov <- round(sqrt( sum(apply(coverageRB, 1, var, na.rm = TRUE) / 
                            sum(!is.na(coverageRB[1, ]))) / length(result)^2 ), 4)
  result2 <- round(rowMeans(coverageS, na.rm = TRUE), 3)
  avr.cov2 <- round(mean(result2), 3)
  se.cov2 <- round(sqrt(apply(coverageS, 1, var, na.rm = TRUE) / 
                        sum(!is.na(coverageS[1, ]))), 4)
  effective.n <- sum(!is.na(coverageS[1, ]))

  # plotting coverage graph
  par(xaxs = "r", yaxs = "r", mai = c(1, 0.9, 1, 0.3), las = 1)
  n.units <- length(gbp.object$se)
  plot(1 : length(gbp.object$se), result, ylim = c(0.5, 1), col = 2,
       ylab = "Coverage rate estimate",
       xlab = paste("Group_ j", ", ", "j = 1, ...,", n.units), 
       main = "Estimated coverage rate for each group",
       lwd = 3, lty = 1, yaxt = "n", xaxt = "n")
  axis(2, at = seq(0.5, 1, by = 0.05), labels = TRUE)

  if (n.units <= 10) {
    axis(1, at = seq(1, n.units, by = 1), labels = TRUE)
  } else if (n.units <= 30) {
    axis(1, at = seq(1, n.units, by = 2), labels = TRUE)
  } else if (n.units <= 50) {
    axis(1, at = seq(1, n.units, by = 5), labels = TRUE)
  } else {
    axis(1, at = seq(1, n.units, by = as.integer(n.units / 10)), labels = TRUE)
  }

  abline(h = gbp.object$confidence.lvl)

  if (is.na(gbp.object$prior.mean) & missing(mean.PriorDist)) {
    if (gbp.object$model == "gr") {
      legend("bottomleft", c(paste("Model: Normal-Normal"), 
                             "Red circles: RB coverage estimates",
                             paste("# of simulations per group: ", effective.n),
                             paste("A for data generation: ", round(A, 2)), 
                             paste("beta", 0 : (length(betas) - 1), " for data generation: ", round(betas, 3), 
                                   sep = ""), 
                             paste("Overall coverage estimate: ", avr.cov),
                             paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
      case <- 1
    } else {
      modelspec <- ifelse(gbp.object$model == "br", "Binomial-Beta", "Poisson-Gamma")
      legend("bottomleft", c(paste("Model: ", modelspec), 
                             "Red circles: RB coverage estimates",
                             paste("# of simulations per group: ", effective.n),
                             paste("r for data generation: ", round(r, 2)), 
                             paste("beta", 0 : (length(betas) - 1), " for data generation: ", round(betas, 3), 
                                   sep = ""), 
                             paste("Overall coverage estimate: ", avr.cov),
                             paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
      case <- 2   
    }
  } else if (is.na(gbp.object$prior.mean) & !missing(mean.PriorDist)) {
    if (gbp.object$model == "gr") {
      legend("bottomleft", c(paste("Model: Normal-Normal"),
                             "Red circles: RB coverage estimates",
                             paste("# of simulations per group: ", effective.n),
                             paste("A for data generation: ", round(A, 2)), 
                             paste("Known prior mean: ", round(priormeanused, 2)), 
                             paste("Overall coverage estimate: ", avr.cov),
                             paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
      case <- 3
    } else {
      modelspec <- ifelse(gbp.object$model == "br", "Binomial-Beta", "Poisson-Gamma")
      legend("bottomleft", c(paste("Model: ", modelspec), 
                             "Red circles: RB coverage estimates",
                             paste("# of simulations per group: ", effective.n),
                             paste("r for data generation: ", round(r, 2)), 
                             paste("Known prior mean: ", round(priormeanused, 2)), 
                             paste("Overall coverage estimate: ", avr.cov),
                             paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
      case <- 4
    }

  } else if (!is.na(gbp.object$prior.mean) & !missing(mean.PriorDist)) {  # if prior mean is assigned
    if (gbp.object$model == "gr") {
      legend("bottomleft", c(paste("Model: Normal-Normal"),
                             "Red circles: RB coverage estimates",
                             paste("# of simulations per group: ", effective.n),
                             paste("A for data generation: ", round(A, 2)), 
                             paste("Known prior mean: ", round(priormeanused, 2)), 
                             paste("Overall coverage estimate: ", avr.cov),
                             paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
      case <- 5
    } else {
      modelspec <- ifelse(gbp.object$model == "br", "Binomial-Beta", "Poisson-Gamma")
      legend("bottomleft", c(paste("Model: ", modelspec), 
                             "Red circles: RB coverage estimates",
                             paste("# of simulations per group: ", effective.n),
                             paste("r for data generation: ", round(r, 2)), 
                             paste("Known prior mean: ", round(priormeanused, 2)), 
                             paste("Overall coverage estimate: ", avr.cov),
                             paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
      case <- 6
    }
  } else if (!is.na(gbp.object$prior.mean) & missing(mean.PriorDist)) {  # if prior mean is assigned
    if (gbp.object$model == "gr") {
      legend("bottomleft", c(paste("Model: Normal-Normal"),
                             "Red circles: RB coverage estimates",
                             paste("# of simulations per group: ", effective.n),
                             paste("A for data generation: ", round(A, 2)), 
                             paste("Known prior mean: ", round(priormeanused, 2)), 
                             paste("Overall coverage estimate: ", avr.cov),
                             paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
      case <- 7
    } else {
      modelspec <- ifelse(gbp.object$model == "br", "Binomial-Beta", "Poisson-Gamma")
      legend("bottomleft", c(paste("Model: ", modelspec), 
                             "Red circles: RB coverage estimates",
                             paste("# of simulations per group: ", effective.n),
                             paste("r for data generation: ", round(r, 2)), 
                             paste("Known prior mean: ", round(priormeanused, 2)), 
                             paste("Overall coverage estimate: ", avr.cov),
                             paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
      case <- 8
    }
  }

  if (gbp.object$model == "gr") {
    A.r <- A
  } else {
    A.r <- r
  }

  # print output
  output <- list(coverageRB = result, coverageS = result2, 
                 average.coverageRB = avr.cov, average.coverageS = avr.cov2, 
                 overall.coverageRB = avr.cov, se.overal.coverageRB = avr.se.cov,
                 se.coverageRB = se.cov, se.coverageS = se.cov2, 
                 raw.resultRB = coverageRB, raw.resultS = coverageS,
                 confidence.lvl = gbp.object$confidence.lvl, effective.n = effective.n, 
                 model = gbp.object$model, case = case, betas =  betas, A.r = A.r, 
                 priormeanused = priormeanused)
  return(output)
}

coverage.plot <- function(cov) {

  n.groups <- length(cov$coverageRB)
  effective.n <- cov$effective.n
  if (cov$model == "gr") {
    A <- cov$A.r
  } else {
    r <- cov$A.r
  }
  betas <- cov$betas
  avr.cov <- round(mean(cov$coverageRB), 3)
  avr.se.cov <- round(sqrt(sum(cov$se.coverageRB^2) / n.groups^2), 4)
  priormeanused <- cov$priormeanused

  # plotting coverage graph
  par(xaxs = "r", yaxs = "r", mai = c(1, 0.9, 1, 0.3), las = 1)
  plot(1 : n.groups, cov$coverageRB, ylim = c(0.5, 1), col = 2,
       ylab = "Coverage rate estimate",
       xlab = paste("Group_ j", ", ", "j = 1, ...,", n.groups), 
       main = "Estimated coverage rate for each group",
       lwd = 3, lty = 1, yaxt = "n", xaxt = "n")
  axis(2, at = seq(0.5, 1, by = 0.05), labels = TRUE)

  if (n.groups <= 10) {
    axis(1, at = seq(1, n.groups, by = 1), labels = TRUE)
  } else if (n.groups <= 30) {
    axis(1, at = seq(1, n.groups, by = 2), labels = TRUE)
  } else if (n.groups <= 50) {
    axis(1, at = seq(1, n.groups, by = 5), labels = TRUE)
  } else {
    axis(1, at = seq(1, n.groups, by = as.integer(n.groups / 10)), labels = TRUE)
  }

  abline(h = cov$confidence.lvl)

  if (cov$case == 1) {
    legend("bottomleft", c(paste("Model: Normal-Normal"), 
                           "Red circles: RB coverage estimates",
                           paste("# of simulations per group: ", effective.n),
                           paste("A for data generation: ", round(A, 2)), 
                           paste("beta", 0 : (length(betas) - 1), " for data generation: ", round(betas, 3), 
                                 sep = ""), 
                           paste("Overall coverage estimate: ", avr.cov),
                           paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
  } else if (cov$case == 2) {
    modelspec <- ifelse(cov$model == "br", "Binomial-Beta", "Poisson-Gamma")
    legend("bottomleft", c(paste("Model: ", modelspec), 
                           "Red circles: RB coverage estimates",
                           paste("# of simulations per group: ", effective.n),
                           paste("r for data generation: ", round(r, 2)), 
                           paste("beta", 0 : (length(betas) - 1), " for data generation: ", round(betas, 3), 
                                 sep = ""), 
                           paste("Overall coverage estimate: ", avr.cov),
                           paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")   
  } else if (cov$case == 3) {
    legend("bottomleft", c(paste("Model: Normal-Normal"),
                           "Red circles: RB coverage estimates",
                           paste("# of simulations per group: ", effective.n),
                           paste("A for data generation: ", round(A, 2)), 
                           paste("Known prior mean: ", round(priormeanused, 2)), 
                           paste("Overall coverage estimate: ", avr.cov),
                           paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
  } else if (cov$case == 4) {
    modelspec <- ifelse(cov$model == "br", "Binomial-Beta", "Poisson-Gamma")
    legend("bottomleft", c(paste("Model: ", modelspec), 
                           "Red circles: RB coverage estimates",
                           paste("# of simulations per group: ", effective.n),
                           paste("r for data generation: ", round(r, 2)), 
                           paste("Known prior mean: ", round(priormeanused, 2)), 
                           paste("Overall coverage estimate: ", avr.cov),
                           paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
  } else if (cov$case == 5) {
    legend("bottomleft", c(paste("Model: Normal-Normal"),
                           "Red circles: RB coverage estimates",
                           paste("# of simulations per group: ", effective.n),
                           paste("A for data generation: ", round(A, 2)), 
                           paste("Known prior mean: ", round(priormeanused, 2)), 
                           paste("Overall coverage estimate: ", avr.cov),
                           paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
  } else if (cov$case == 6) {
    modelspec <- ifelse(cov$model == "br", "Binomial-Beta", "Poisson-Gamma")
    legend("bottomleft", c(paste("Model: ", modelspec), 
                           "Red circles: RB coverage estimates",
                           paste("# of simulations per group: ", effective.n),
                           paste("r for data generation: ", round(r, 2)), 
                           paste("Known prior mean: ", round(priormeanused, 2)), 
                           paste("Overall coverage estimate: ", avr.cov),
                           paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
  } else if (cov$case == 7) {
    legend("bottomleft", c(paste("Model: Normal-Normal"),
                           "Red circles: RB coverage estimates",
                           paste("# of simulations per group: ", effective.n),
                           paste("A for data generation: ", round(A, 2)), 
                           paste("Known prior mean: ", round(priormeanused, 2)), 
                           paste("Overall coverage estimate: ", avr.cov),
                           paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
  } else if (cov$case == 8) {
    modelspec <- ifelse(cov$model == "br", "Binomial-Beta", "Poisson-Gamma")
    legend("bottomleft", c(paste("Model: ", modelspec), 
                           "Red circles: RB coverage estimates",
                           paste("# of simulations per group: ", effective.n),
                           paste("r for data generation: ", round(r, 2)), 
                           paste("Known prior mean: ", round(priormeanused, 2)), 
                           paste("Overall coverage estimate: ", avr.cov),
                           paste("SE(overall coverage estimate): ", avr.se.cov)), bty = "n")
  }
}
