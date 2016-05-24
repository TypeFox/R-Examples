# semm.R
#
# created: Okt/20/2014, KN
# last mod: Aug/27/2015, NU

#--------------- main functions ---------------

# Calculate mu of multivariate normal distribution for joint vector of
# indicators (See equation 4 in Jedidi, Jagpal & DeSarbo, 1997).
# The order is (x, y) as opposed to the paper.
mu_semm <- function(matrices) {

    #check_filled(matrices)
    # TODO: Remove?

    mu.y <- matrices$nu.y + matrices$Lambda.y %*% solve(matrices$Beta) %*%
            (matrices$alpha + matrices$Gamma %*% matrices$tau)
    mu.x <- matrices$nu.x + matrices$Lambda.x %*% matrices$tau
    mu <- rbind(mu.x, mu.y) # vertical vector

    mu
}

# Calculate sigma of multivariate normal distribution for joint vector of
# indicators (y, x). (See equation 5 in Jedidi, Jagpal & DeSarbo, 1997)
# The order is (x, y), as opposed to the paper. Therefore the rows and cols in
# sigma are switched.
sigma_semm <- function(matrices) {

    #check_filled(matrices)

    # Lambda.y * B^-1
    Ly.Binv <- matrices$Lambda.y %*% solve(matrices$Beta)

    s22 <- Ly.Binv %*% (matrices$Gamma %*% matrices$Phi %*% t(matrices$Gamma) +
                        matrices$Psi) %*% t(Ly.Binv) + matrices$Theta.e

    s21 <- Ly.Binv %*% matrices$Gamma %*% matrices$Phi %*% t(matrices$Lambda.x)
    s12 <- t(s21)
    s11 <- matrices$Lambda.x %*% matrices$Phi %*% t(matrices$Lambda.x) + matrices$Theta.d
    sigma <- rbind(cbind(s11,s12), cbind(s21, s22))

    # TODO check if this warning is really necessary
    if (!isSymmetric(sigma)) warning("Sigma is not symmetric. This is probably due to numerical calculation.")
    tryCatch(solve(sigma), error = function(e) stop("Sigma is singular."))

    sigma
}

# Expectation step of the EM-algorithm (see Jedidi, Jagpal & DeSarbo, 1997)
estep_semm <- function(model, parameters, data) {

  model.filled <- fill_model(model=model, parameters=parameters)

  P <- NULL
  for (c in seq_len(model$info$num.classes)) {
    # class weight
    w.c <- model$info$w[c]

    p.ij <- w.c * dmvnorm(data, mean=mu_semm(model.filled$matrices[[c]]),
                          sigma=sigma_semm(model.filled$matrices[[c]]))
    if (sum(p.ij) == 0) stop("Posterior probability could not be calculated properly. Choose different starting parameters.")
    P <- cbind(P, p.ij, deparse.level=0)
  }
  P <- P / rowSums(P)
  P
}

# Negative log likelihood function which will be optimized in M-step (see below)
loglikelihood_semm <- function(parameters, matrices, data, p, w) {
    # fill matrices
    for (i in seq_along(matrices)) {
        matrix.i <- matrices[[i]]
        # number of NA's in matrix
        num.na <- length(matrix.i[is.na(matrix.i)])
        if (num.na > 0) {
            matrix.i[is.na(matrix.i)] <- parameters[1:num.na]
            parameters <- parameters[-(1:num.na)]
            matrices[[i]] <- matrix.i
        }
    }
    tryCatch({ matrices$Psi <- fill_symmetric(matrices$Psi) },
                                       error=function(e) e,
                                       warning=function(w) w)
    tryCatch({ matrices$Phi <- fill_symmetric(matrices$Phi) },
                                       error=function(e) e,
                                       warning=function(w) w)

    N <- nrow(data)
    N.c <- sum(p)
    mu <- mu_semm(matrices)
    sigma <- sigma_semm(matrices)
    T.c <- 1/N.c * Reduce('+', lapply(1:N, function(i)(
                                       p[i] * (data[i,]-mu) %*%
                                       t(data[i,]-mu))))
    res <- 1/2 * N.c * (log(det(sigma)) + sum(diag(T.c %*% solve(sigma))) -
                        2*log(w))
    res
}

# Negative log likelihood function for maximization of all classes at once
loglikelihood_semm_constraints <- function(parameters, model, data, P) {
    model.filled <- fill_model(model=model, parameters=parameters)
    N <- nrow(data)
    res <- 0

    for (c in seq_len(model$info$num.classes)) {
        w.c <- model$info$w[c]
        N.c <- sum(P[,c])
        mu.c <- mu_semm(model.filled$matrices[[c]])
        sigma.c <- sigma_semm(model.filled$matrices[[c]])
        T.c <- 1/N.c * Reduce('+', lapply(1:N, function(i)(
                                           P[i,c] * (data[i,] - mu.c) %*%
                                           t(data[i,] - mu.c))))
        res <- res + (1/2 * N.c * (log(det(sigma.c)) + sum(diag(T.c %*%
                                           solve(sigma.c))) - 2*log(w.c)))
    }

    res
}


# Maximization step of the EM-algorithm (see Jedidi, Jagpal & DeSarbo, 1997)
mstep_semm <- function(model, parameters, data, P, neg.hessian=FALSE,
                       optimizer=c("nlminb", "optim"),
                       max.mstep, control=list(), ...) {

  optimizer <- match.arg(optimizer)

  if (model$info$constraints == "direct1") {
    # Maximizing each class separately
    num.classes <- model$info$num.classes
    class.pars <- get_class_parameters(model, parameters)

    est <- lapply(seq_len(num.classes), function(c) {
            if (optimizer == "nlminb") {
              if (is.null(control$iter.max)) {
                control$iter.max <- max.mstep
              } else {
                warning("iter.max is set for nlminb. max.mstep will be ignored.")
              }
              suppress_NaN_warnings(
                res <- nlminb(start=class.pars[[c]],
                              objective=loglikelihood_semm,
                              data=data, matrices=model$matrices[[c]],
                              p=P[,c], w=model$info$w[[c]],
                              upper=model$info$bounds$upper[[c]],
                              lower=model$info$bounds$lower[[c]],
                              control=control, ...)
              )
            } else {
              if (is.null(control$maxit)){
                control$maxit <- max.mstep
              } else {
                warning("maxit is set for optim. max.mstep will be ignored.")
              }
              res <- optim(par=class.pars[[c]],
                           fn=loglikelihood_semm, data=data,
                           matrices=model$matrices[[c]],
                           p=P[,c], w=model$info$w[[c]],
                           upper=model$info$bounds$upper[[c]],
                           lower=model$info$bounds$lower[[c]],
                           method="L-BFGS-B", control=control, ...)
            }
    })
    if (optimizer == "optim") {
      for (c in seq_len(num.classes)) {
        names(est[[c]]) <- gsub("value", "objective", names(est[[c]]))
      }
    }
    res <- list(objective=0)
    for (c in seq_len(num.classes)) {
      res$par[[c]] <- est[[c]]$par
      res$objective <- res$objective + est[[c]]$objective
      res$convergence[[c]] <- est[[c]]$convergence
      res$iterations <- est[[c]]$iterations
    }
    names(res$par) <- paste0("class", seq_len(num.classes))

    if (neg.hessian == TRUE) {
      for (c in seq_len(num.classes)) {
        res$hessian[[c]] <- fdHess(pars=est[[c]]$par,
                                   fun=loglikelihood_semm,
                                   matrices=model$matrices[[c]],
                                   data=data, p=P[,c],
                                   w=model$info$w[[c]])$Hessian
      }
    names(res$hessian) <- paste0("class", seq_len(num.classes))
    }
    res

  } else {
    # Maximization of all classes together
    if (optimizer == "nlminb") {
      if (is.null(control$iter.max)) {
        control$iter.max <- max.mstep
      } else {
        warning("iter.max is set for nlminb. max.mstep will be ignored.")
      }
      suppress_NaN_warnings(
        est <- nlminb(start=parameters,
                      objective=loglikelihood_semm_constraints, data=data,
                      model=model, P=P,
                      upper=unlist(model$info$bounds$upper),
                      lower=unlist(model$info$bounds$lower),
                      control=control, ...)
      )
    } else {
      if (is.null(control$maxit)){
        control$maxit <- max.mstep
      } else {
        warning("maxit is set for optim. max.mstep will be ignored.")
      }
      est <- optim(par=parameters, fn=loglikelihood_semm_constraints,
                   model=model, data=data, P=P,
                   upper=unlist(model$info$bounds$upper),
                   lower=unlist(model$info$bounds$lower),
                   method="L-BFGS-B", control=control, ...)
      # fit est to nlminb output
      names(est) <- gsub("value", "objective", names(est))
    }
    if (neg.hessian == TRUE) {
      est$hessian <- fdHess(pars=est$par,
                            fun=loglikelihood_semm_constraints,
                            model=model, data=data, P=P)$Hessian
    }
    est
  }
}

#--------------- helper functions ---------------

# Make a list of class specific parameter vectors
get_class_parameters <- function(model, parameters) {

  mod.filled <- fill_model(model=model, parameters=parameters)

  dat.filled <- as.data.frame(mod.filled)
  dat <- as.data.frame(model)

  class.pars <- list()
  for (class in names(model$matrices)) {
    class.pars[[class]] <- dat.filled[is.na(dat[, class]), class]
  }
  class.pars
}

# Obtain two neg. Hessians for several classes for indirect and direct2
# approach
get_hessian <- function(object) {
  if (is.null(object$neg.hessian)) {
    hessian.l <- NULL
  } else {
    if (object$info$constraints == "direct2") {
      names.class1 <- names(object$coefficients$class1)
      for (class in names(object$coef)[-1]) {
        names.class <- names(coef(object)[[class]])
        par.names <- c(paste0("class1.",names.class1),
          paste0(class,".",names.class[grep("Gamma", names.class)]))
        # Is there a Beta?
        if (object$info$num.eta > 1) {
        par.names <- c(par.names, paste0(class,".",names.class[grep("Beta",
          names.class)]))
        }
        par.names <- c(par.names, paste0(class,".",names.class[grep("Psi", names.class)]),
          paste0(class,".",names.class[grep("Phi", names.class)]),
          paste0(class,".",names.class[grep("alpha", names.class)]),
          paste0(class,".",names.class[grep("tau", names.class)]))
        # Is there an Omega?
        if (object$model.class == "nsemm") {
          par.names <- c(par.names, paste0(class,".",names.class[grep("Omega", names.class)]))
        }
      }
    } else {
      names.class1 <- names(object$coefficients$class1)
      for (class in names(object$coef)[-1]) {
        names.class <- names(coef(object)[[class]])
        par.names <- c(paste0("class1.",names.class1),
          paste0(class,".",names.class[grep("Phi", names.class)]),
          paste0(class,".",names.class[grep("tau", names.class)]))
      }
    }
    dimnames(object$neg.hessian) <- list(par.names, par.names)
    hessian.l <- list()
    for (class in names(object$coef)) {
      hessian.l[[class]] <- object$neg.hessian[grep("class1",
        rownames(object$neg.hessian)), grep("class1",
        colnames(object$neg.hessian))]
      dimnames(hessian.l[[class]]) <- list(names(object$coefficients$class1),
        names(object$coefficients$class1))
    }
    for (class in names(hessian.l)[-1]) {
      par.names.c <- par.names[grep(class, par.names)]
      hessian.l[[class]][gsub(paste0("^",class,".(.*)$"), "\\1",
      par.names.c), gsub(paste0("^",class,".(.*)$"), "\\1", par.names.c)] <-
        object$neg.hessian[grep(class, rownames(object$neg.hessian)),
        grep(class, colnames(object$neg.hessian))]
    }
    }
  hessian.l
}

# Make a list of class specific parameter names
get_class_parnames <- function(model) {

  dat <- as.data.frame(model)

  class.names <- list()
  for (class in names(model$matrices)) {
    class.names[[class]] <- as.character(dat$label[is.na(dat[, class])])
  }
  class.names
}

# Suppress all warnings that contain 'NaN'
suppress_NaN_warnings <- function(expr) {
  withCallingHandlers(expr, warning=function(w) {
                      if (grepl("NaN", conditionMessage(w)))
                        invokeRestart("muffleWarning")
  })
}
