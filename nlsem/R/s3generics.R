# s3generics.R
#
# created Nov/03/2014, KN
# last mod Sep/17/2015, NU

#--------------- main functions ---------------

as.data.frame.singleClass <- as.data.frame.semm <- as.data.frame.nsemm <- function(x, ...) {
    data <- data.frame(
        label = names(unlist(x$matrices$class1)))
    for (c in seq_len(length(x$matrices))) {
        temp <- data.frame(unlist(x$matrices[[c]], use.names=FALSE))
        names(temp) <- paste0("class", c)
        data <- cbind(data, temp)
    }
    attr(data, "constraints") <- x$info$constraints
    data
}

simulate.nsemm <- function(object, nsim=1, seed=NULL, n=400, m=16,
                           parameters, ...) {

  set.seed(seed)

  mod.filled <- fill_model(model=object, parameters=parameters)

  num.classes <- mod.filled$info$num.classes
  w <- mod.filled$info$w

  # simulate n data points for each mixture as lms
  dat.sim <- lapply(seq_len(num.classes), function(c) {
                    simulate(lms_ify(object, c),
                             nsim=nsim, seed=seed, n=n, m=m,
                             parameters=get_class_parameters(object,
                             parameters)[[c]], ...)

  })

  # see simulate.singleClass for explanation
  border <- cumsum(w)
  prob <- runif(n)

  dat <- NULL
  for (i in seq_len(n)) {
    ind <- sum(prob[i] > border) + 1
    dat <- rbind(dat, dat.sim[[ind]][i,])
  }

  dat
}

simulate.semm <- function(object, nsim=1, seed=NULL, n=400, parameters,
                          ...) {

  set.seed(seed)

  mod.filled <- fill_model(model=object, parameters=parameters)

  num.classes <- mod.filled$info$num.classes
  w <- mod.filled$info$w

  # simulate n data points from each mixture distribution
  dat.sim <- lapply(1:num.classes, function(c) {
                    rmvnorm(n,
                            mean=mu_semm(matrices=mod.filled$matrices[[c]]),
                            sigma=sigma_semm(matrices=mod.filled$matrices[[c]]))
                            })

  # see simulate.singleClass for explanation
  border <- cumsum(w)
  prob <- runif(n)

  dat <- NULL
  for (i in seq_len(n)){
    ind <- sum(prob[i] > border) + 1
    dat <- rbind(dat, dat.sim[[ind]][i,])
  }

  dat
}

simulate.singleClass <- function(object, nsim=1, seed=NULL, n=400, m=16,
                                 parameters, ...) {

  if (object$info$num.eta > 1) {
    stop("Cannot simulate data for a model with more than one eta, yet.")
  }

  # set seed
  set.seed(seed)

  # Gauss-Hermite quadrature
  k <- get_k(object$matrices$class1$Omega)
  if (k != 0){
      quad <- quadrature(m=m, k=k)
      V <- quad$n
      w <- quad$w
  } else {
      V <- 0
      w <- 1
      # do not need mixtures, if I do not have interactions
  }

  parameters <- convert_parameters_singleClass(object, parameters)
  names(object$matrices$class1)[grep("Phi", names(object$matrices$class1))] <- "A"

  mod.filled <- fill_model(model=object, parameters=parameters)

  # simulate n data points from each mixture distribution
  dat.sim <- lapply(seq_along(w), function(i){
                      rmvnorm(n, 
                      mean=mu_lms(model=mod.filled, z=V[i,]),
                      sigma=sigma_lms(model=mod.filled, z=V[i,]))
                      })
  dat <- dat.sim[[1]]

  # decide which data points from each mixture should be included in
  # simulated data set: weights give intervall borders between 0 and 1;
  # we draw random numbers from a uniform distribution and check in what
  # intervall they lie: the ith element from that distribution will be
  # taken and put into the data frame
  if (k != 0){
    border <- cumsum(w)
    prob <- runif(n)

    dat <- NULL
    for (i in seq_len(n)){
      ind <- sum(prob[i] > border) + 1
      dat <- rbind(dat, dat.sim[[ind]][i,])
    }
  }
  dat
}

print.singleClass <- print.semm <- print.nsemm <- function(x, ...) {
  cat("Model of class", class(x), "\n\n")

  cat("Number of latent endogenous variables:", x$info$num.eta, "(with", x$info$num.y, "indicators)\n")
  cat("Number of latent exogenous variables:", x$info$num.xi, "(with", x$info$num.x, "indicators)\n")
  if (class(x) != "singleClass") {
    cat("Number of latent classes:", x$info$num.classes, "\n")
  }
  cat("\nStructural model:\n")
  for (class in names(x$matrices)) {
    cat("------------------------------------------------------------\n")
    if (length(names(x$matrices)) > 1) {
      cat(paste0("Specifications for ", class, ":"), "\n")
    }
    # rel.lat
    rel.lat <- get_rel.lat(x$matrices[[class]][["Gamma"]],
      x$matrices[[class]][["Beta"]])
    for (r in rel.lat[[1]]) {
      cat(r, "\n")
    }
    for (r in rel.lat[[2]]) {
      cat(r, "\n")
    }
    # interaction
    if (class(x) != "semm") {
      interaction <- get_interaction(x$matrices[[class]][["Omega"]])
      for (i in interaction) {
        cat(i, "\n")
      }
    }
  }
  cat("------------------------------------------------------------\n")

  cat("\nMeasurement model:\n")
  for (class in names(x$matrices)) {
    cat("------------------------------------------------------------\n")
    if (length(names(x$matrices)) > 1) {
      cat(paste0("Specifications for ", class, ":"), "\n")
    }
    xi <- get_xi(x$matrices[[class]][["Lambda.x"]])
    for (i in xi) {
      cat(i, "\n")
    }
    eta <- get_eta(x$matrices[[class]][["Lambda.y"]])
    for (e in eta) {
      cat(e, "\n")
    }
  }
  cat("------------------------------------------------------------\n")

  if (class(x) != "singleClass") {
    constraints <- x$info$constraints
    switch(EXPR = constraints,
      indirect = {
        cat(paste0("\nModel constraints set to ", constraints, ":\n"))
        cat("Parameters of all classes will be set to equal except for taus and Phi.\n\n")
      },
      direct1 = {
        cat(paste0("\nModel constraints set to ", constraints, ":\n"))
        cat("All parameters of all classes will be estimated separately.\n\n")
      },
      direct2 = {
        cat(paste0("\nModel constraints set to ", constraints, ":\n"))
        cat("Parameters of measurement model will be equal for all classes.\n")

      }
    )     # end switch
  }

}

print.emEst <- function(x, digits=3, ...) {
  cat("Fitted model of class", x$model.class, "with", x$info$num.xi +
    x$info$num.eta, "latent variables and", x$info$num.x + x$info$num.y,
    "indicators.\n\n")

  if (x$model.class != "singleClass") {
    constraints <- x$info$constraints
    switch(EXPR = constraints,
      indirect = {
        cat(paste0("Model constraints set to ", constraints, ":\n"))
        cat("Parameters of all classes equal except for taus and Phi.\n\n")
      },
      direct1 = {
        cat(paste0("Model constraints set to ", constraints, ":\n"))
        cat("All parameters of all classes different.\n\n")
      },
      direct2 = {
        cat(paste0("Model constraints set to ", constraints, ":\n"))
        cat("Parameters of measurement model equal for all classes.\n\n")

      }
    )     # end switch

    cat("Class weights:", round(x$info$w, digits=digits), "\n\n")
  }
  cat("Estimated parameters:\n")
  print(round(as.data.frame(x$coef), digits=digits))


}

summary.emEst <- function(object, print.likelihoods = FALSE, ...) {

  # estimates
  est <- object$coefficients

  # standard errors
  if (is.numeric(est)) {
    s.error <- calc_standard_error(object$neg.hessian)
    tvalue <- est / s.error
    pvalue <- 2 * pnorm(-abs(tvalue))
    est.table <- cbind(est, s.error, tvalue, pvalue)
    dimnames(est.table)  <- list(names(est), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  } else {
    if (object$info$constraints == "direct1") {
      neg.hessian <- object$neg.hessian
    } else {
      neg.hessian <- get_hessian(object)
    }
    est.table <- Reduce('rbind', lapply(seq_along(est), function(c) {
      s.error <- calc_standard_error(neg.hessian[[c]])
      tvalue <- est[[c]] / s.error
      pvalue <- 2 * pnorm(-abs(tvalue))
      est.table <- cbind(est[[c]], s.error, tvalue, pvalue)
      dimnames(est.table)  <- list(paste0("class", c, ".", names(est[[c]])),
        c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
      est.table
    }))
  }

  # loglikelihoods
  iterations   <- length(object$loglikelihoods) 
  abs.change   <- c(0, diff(object$loglikelihoods))
  rel.change   <- rel_change(object$loglikelihoods)
  logLik.table <- cbind(object$loglikelihoods, abs.change, rel.change)
  dimnames(logLik.table) <- list(1:iterations, c("loglikelihood", "difference", "relative change"))

  ans <- list(model=object$model.class,
              estimates=est.table,
              iterations=iterations,
              finallogLik=object$objective,
              logLikelihoods=logLik.table)

  if (object$model.class == "semm" || object$model.class == "nsemm") {
    ans$class.weights <- object$info$w
  }

  ans$print.likelihoods <- print.likelihoods

  class(ans) <- "summary.emEst"

  ans
}

print.summary.emEst <- function(x, digits=max(3, getOption("digits") - 3),
                                cs.ind=2:3, ...) {
    
  cat("\nSummary for model of class", x$model, "\n")
  cat("\nEstimates:\n")
  printCoefmat(x$estimates, digits=digits, cs.ind=cs.ind, ...)
  cat("\nNumber of iterations:", x$iterations, "\nFinal loglikelihood:",
    round(x$finallogLik, 3), "\n") 
  if (x$model == "semm" || x$model == "nsemm"){
    cat("\nClass weights:", round(x$class.weights, digits), "\n\n")
  }
  if (x$print.likelihoods) {
    cat("\n", "\nLikelihoods:\n")
    printCoefmat(x$logLikelihoods, digits=6, cs.ind=2:3, ...)
  }

}

print.qmlEst <- function(x, digits=3, ...) {
  cat("Fitted model of class", x$model.class, "with", x$info$num.xi +
    x$info$num.eta, "latent variables and", x$info$num.x + x$info$num.y,
    "indicators.\n\n")

  cat("Estimated parameters:\n")
  print(round(as.data.frame(x$coef), digits=digits))
}

summary.qmlEst <- function(object, ...) {

  # estimates
  est <- object$coefficients

  # standard errors
  s.error <- calc_standard_error(object$neg.hessian)
  tvalue <- est / s.error
  pvalue <- 2 * pnorm(-abs(tvalue))
  est.table <- cbind(est, s.error, tvalue, pvalue)
  dimnames(est.table)  <- list(names(est), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

  ans <- list(model=object$model.class,
              estimates=est.table,
              iterations=object$iterations,
              finallogLik=object$objective)

  class(ans) <- "summary.qmlEst"

  ans
}

print.summary.qmlEst <- function(x, digits=max(3, getOption("digits") - 3),
                                 cs.ind=2:3, ...) {
    
  cat("\nSummary for model of class", x$model, "estimated with QML\n")
  cat("\nEstimates:\n")
  printCoefmat(x$estimates, digits=digits, cs.ind=cs.ind, ...)
  cat("\nNumber of iterations:", x$iterations,
      "\nFinal loglikelihood:", round(x$finallogLik, 3), "\n") 
}

logLik.emEst <- logLik.qmlEst <- function(object, ...){
  if(length(list(...)))
      warning("extra arguments discarded")

  out <- object$objective
  attr(out, "df") <- length(unlist(object$coef))
  class(out) <- "logLik"
  out
}

anova.emEst <- anova.qmlEst <- function(object, ..., test=c("Chisq", "none")) {
  # Adapted from anova.polr by Brian Ripley
  
  test <- match.arg(test)
  dots <- list(...)
  if (length(dots) == 0)
      stop('anova is not implemented for a single "emEst" object')

  mlist <- list(object, ...)
  if (any(!sapply(mlist, function(x) x$model.class == "singleClass"))) {
      stop('Likelihood Ratio Test only meaningful for models of class "singleClass".')
  }

  nlist <- sapply(mlist, function(x) x$info$n)
  if (!all(nlist == object$info$n)) {
      stop("SEM have not all been fitted to the same data set.")
  }

  names(mlist) <- c(deparse(substitute(object)),
              as.character(substitute(...[]))[2:length(mlist)])
  if (any(!sapply(mlist, inherits, "emEst")))
      stop('not all objects are of class "emEst"')
  nt <- length(mlist)

  dflist <- sapply(mlist, function(x) length(unlist(x$coef)))

  s <- order(dflist)
  mlist <- mlist[s]
  dflist <- dflist[s]

  lls <- sapply(mlist, function(x) -2*x$objective)
  tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
  df <- c(NA, diff(dflist))
  x2 <- c(NA, -diff(lls))
  pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
  out <- data.frame(Model=names(mlist), Resid.df=dflist, Deviance=lls,
                    Test=tss, Df=df, LRtest=x2, Prob=pr)
  names(out) <- c("Model", "Numb. coef", "-2logLik", "Test",
                  "   Df", "LR stat.", "Pr(>Chi)")
  rownames(out) <- 1:nt
  if (test == "none") out <- out[, 1:6]
  class(out) <- c("Anova", "data.frame")
  attr(out, "heading") <- "Chi Square test statistic for (nonlinear) SEM\n"
  out
}

AIC.emEst <- AIC.qmlEst <- function(object, ..., k=2) {

  dots <- list(...)
  if (length(dots) == 0){
      out <- as.numeric(-2*logLik(object) + k*length(object$coef))
  } else {
    mlist <- list(object, ...)
    names(mlist) <- c(deparse(substitute(object)),
                      as.character(substitute(...[]))[2:length(mlist)])
    nt <- length(mlist)

    dflist <- sapply(mlist, function(x) length(unlist(x$coef)))

    aic <- sapply(mlist, function(x) -2*logLik(x) + k*length(unlist(x$coef)))
    s <- order(aic, decreasing=TRUE)
    aic <- aic[s]
    dflist <- dflist[s]
    out <- data.frame(df=dflist, AIC=aic)
    rownames(out) <- names(mlist)
  }
  out
}

BIC.emEst <- BIC.qmlEst <- function(object, ...) {

  dots <- list(...)
  if (length(dots) == 0){
      out <- as.numeric(-2*logLik(object) +
             log(object$info$n)*length(unlist(object$coef)))
  } else {
    mlist <- list(object, ...)
    names(mlist) <- c(deparse(substitute(object)),
                      as.character(substitute(...[]))[2:length(mlist)])
    nt <- length(mlist)

    dflist <- sapply(mlist, function(x) length(unlist(x$coef)))

    bic <- sapply(mlist, function(x) -2*logLik(x) + log(object$info$n)*length(unlist(x$coef)))
    s <- order(bic, decreasing=TRUE)
    bic <- bic[s]
    dflist <- dflist[s]
    out <- data.frame(df=dflist, BIC=bic)
    rownames(out) <- names(mlist)
  }
  out
}

plot.emEst <- function(x, y, ...) {

  plot(x$loglikelihoods, type="l", xlab="Number of iterations", 
       ylab="log likelihood", axes=F, ...)
  axis(1, at=1:length(x$loglikelihoods))
  axis(2)
  box()
}

#--------------- helper functions ---------------

# calculates relative change defined as absolute difference divided by
# maximum absolute value
rel_change <- function(x) {
    
  if (length(x) == 1) {
      rel.change <- 0
  } else {
    rel.change <- numeric(length(x))
    for (i in 2:length(rel.change)){
      rel.change[i] <- abs(x[i-1]-x[i])/max(abs(x[i-1]), abs(x[i]))
    }
  }
  rel.change
}

calc_standard_error <- function(neg.hessian) {
  s.error <- tryCatch({
      sqrt(diag(solve(neg.hessian)))
    }, error=function(e) {
      NA
    }, warning=function(w) {
       if (grepl("NaN", conditionMessage(w))) {
         suppressWarnings(sqrt(diag(solve(neg.hessian))))
      } else{
         sqrt(diag(solve(neg.hessian)))
      }
    })
    if (all(is.na(s.error))) warning("Standard errors could not be computed, because negative Hessian was either not available or singular.")
    if (any(is.nan(s.error))) warning("Standard errors for some coefficients could not be computed.") 
  s.error
}

# Extracting information from matrices to print them in model summary:
# print.singleClass/semm/nsemm
get_xi <- function(Lambda.x) {

  ind <- list()
  for (i in seq_len(ncol(Lambda.x))) {
    ind[[i]] <- seq_len(nrow(Lambda.x))[-which(Lambda.x[,i] == 0)]
    if (length(ind[[i]]) == 0) {
      ind[[i]] <- 1
    }
  }
  xs <- sapply(ind, function(x) paste0("x", x, collapse=" + "))
  res <- paste(paste0("xi", seq_len(ncol(Lambda.x)), " ="), xs)
  res

}

get_eta <- function(Lambda.y) {

  ind <- list()
  for (i in seq_len(ncol(Lambda.y))) {
    ind[[i]] <- seq_len(nrow(Lambda.y))[-which(Lambda.y[,i] == 0)]
    if (length(ind[[i]]) == 0) {
      ind[[i]] <- 1
    }
  }
  xs <- sapply(ind, function(x) paste0("y", x, collapse=" + "))
  res <- paste(paste0("eta", seq_len(ncol(Lambda.y)), " ="), xs)
  res


}

get_rel.lat <- function(Gamma, Beta) {

  # eta and xi
  ind <- which(is.na(Gamma), arr.ind=T)
  g <- NULL
    if (nrow(ind) > 0) {
    ind.fill <- matrix(c(paste0("eta", ind[,"row"]), paste0("xi",
      ind[,"col"])), nrow=dim(ind)[1], ncol=dim(ind)[2])
    for (i in unique(ind.fill[,1])) {
      g <- c(g, paste(i, "=", paste(ind.fill[ind.fill[,1] == i, 2],
        collapse=" + ")))
    }
  }

  # eta and eta
  ind <- which(is.na(Beta), arr.ind=T)
  b <- NULL
  if (nrow(ind) > 0) {
    ind.fill <- matrix(c(paste0("eta", ind[,"row"]), paste0("eta",
      ind[,"col"])), nrow=dim(ind)[1], ncol=dim(ind)[2])
    for (i in unique(ind.fill[,1])) {
      b <- c(b, paste(i, "=", paste(ind.fill[ind.fill[,1] == i, 2],
        collapse=" + ")))
    }
  }
  list(g, b)
}

get_interaction <- function(Omega) {

  if (anyNA(Omega)) {
    if (is.matrix(Omega) || (length(dim(Omega)) == 3 & dim(Omega)[3] == 1)) {
      Omega <- matrix(Omega, nrow(Omega))
      ind <- which(is.na(Omega), arr.ind=T)
      ind.fill <- matrix(c(paste0("xi", ind[,"row"]), paste0("xi",
        ind[,"col"])), nrow=dim(ind)[1], ncol=2)
      int <- apply(ind.fill, 1, function(x) paste(x, collapse=":"))
      o <- paste("eta1 =", paste(int, collapse=" + "))
    } else {
      ind <- which(is.na(Omega), arr.ind=T)
      ind.fill <- matrix(c(paste0("xi", ind[,"dim1"]), paste0("xi",
        ind[,"dim2"])), nrow=dim(ind)[1], ncol=2)
      int <- apply(ind.fill, 1, function(x) paste(x, collapse=":"))
      o <- NULL
      for (i in unique(ind[, "dim3"])) {
        o <- c(o, paste0("eta", i, " = ", paste(int[ind[, "dim3"] == i],
          collapse=" + ")))
      }
    }
  } else {
  o <- NULL
  }
  o
}

