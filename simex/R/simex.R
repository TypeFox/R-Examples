# Aenderungen gekennzeichnet über #**Heidi**#
# Nenne momentan Funktion "testsimex" zum testen, ob sie genau so gut funktioniert, wie simex aus bisherigem Package
# Soll spaeter aber wieder "simex" heissen

##' @param model the naive model
##' @param SIMEXvariable character or vector of characters containing the names ofthe variables with measurement error
##' @param measurement.error matrix of standard deviations of the known measurement errors
##' @param lambda vector of lambdas for which the simulation step should be done (without 0)
##' @param B number of iterations for each lambda
##' @param fitting.method fitting method \code{linear}, \code{quadratic}, \code{nonlinear} are implemented. (first 4 letters are enough)
##' @param jackknife.estimation specifying the extrapolation method for jackknife variance estimation. Can be set to \code{FALSE} if it should not be performed
##' @param asymptotic logical, indicating if asymptotic variance estimation should be done, in the naive model the option \code{x = TRUE} has to be set
##' @param x object of class 'simex'
##' @param digits number of digits to be printed
##' @param object object of class 'simex'
##' @param xlab optional name for the X-Axis
##' @param ylab vector containing the names for the Y-Axis
##' @param ask logical. If \code{TRUE}, the user is asked for input, before a new figure is drawn
##' @param show vector of logicals indicating for wich variables a plot should be produced
##' @param newdata optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used
##' @param \dots arguments passed to other functions

simex <-
  function (
    model,                          # the naive model
    SIMEXvariable,                  # vector of names for the SIMEXvariables
    measurement.error,              # vector of standard deviations of the measurement errors
    lambda = c(0.5, 1, 1.5, 2),     # vector of lambda which contains the values for the ...
    B = 100,                        # numeric: number of simulations to be made in each step
    fitting.method = "quadratic",   # fitting method for the extrapolation step
    jackknife.estimation = "quadratic", # extrapolation function for the variance estimation
    asymptotic = TRUE)              # logical: do asymptotic variance estimation
  {
    fitting.method <- substr(fitting.method, 1, 4)
    if (!any(fitting.method == c("quad", "line", "nonl"))) {
      warning("Fitting Method not implemented. Using: quadratic", call. = FALSE)
      fitting.method <- "quad"
    }
    if (jackknife.estimation != FALSE)
      jackknife.estimation <- substr(jackknife.estimation, 1, 4)
    if (!any(jackknife.estimation == c("quad", "line", "nonl", FALSE))) {
      warning("Fitting Method (jackknife) not implemented. Using: quadratic",
              call. = FALSE)
      jackknife.estimation <- "quad"
    }
    if (!is.character(SIMEXvariable))
      stop("SIMEXvariable must be character", call. = FALSE)
    if (any(lambda <= 0)) {
      warning("lambda should not contain 0 or negative values. 0 or negative values will be ignored",
              call. = FALSE)
      lambda <- lambda[lambda >= 0]
    }
    if (!any(names(model) == "x") && asymptotic)
      stop("The option x must be enabled in the naive model for asymptotic variance estimation",
           call. = FALSE)
#**Heidi**#
    measurement.error <- as.matrix(measurement.error)
    SIMEXdata <- model$model #das war vorher nach "for (j in 1:B)"
    if (NROW(measurement.error) != NROW(SIMEXdata) && NROW(measurement.error) == 1){
      measurement.error <- matrix(measurement.error, nrow = NROW(SIMEXdata), ncol = NCOL(measurement.error), byrow = TRUE)
    }
    if (NROW(measurement.error) != NROW(SIMEXdata) && NROW(measurement.error) != 1){
      stop("NROW(measurement.error) must be either 1 or must take the number of rows of the data used.",
           call. = FALSE)
    }
#    if (length(measurement.error) != length(SIMEXvariable))
#      stop("SIMEXvariable and measurement.error must have the same length",
#           call. = FALSE)
#********#
#**Heidi**# 
    # vorher konnte es auch nicht 0 sein, jetzt darf es das
    if (any(measurement.error < 0))
      stop("measurement.error is negative", call. = FALSE)
    # es darf aber nicht konstant 0 sein
    any0 <- ( apply(measurement.error, 2, all.equal, current = rep(0, times = NROW(measurement.error))) == TRUE )
    if ( sum(any0) > 0 ) 
      stop("measurement.error is constant 0 in column(s) ", which(any0), call. = FALSE)
    
    # Abbruch falls das Modell nicht lm oder glm bzw. heterosk. me und asymptotic = TRUE
    if ( asymptotic == TRUE & ((class(model)[1] != "glm" & class(model)[1] != "lm") | dim(unique(measurement.error))[1] != 1) )
      stop("Asymptotic is only implemented for naive models of class lm or glm with homoscedastic measurement error.")
#*********#  
    cl <- match.call()
    # defining the vector for the solutions of the simulations
    ncoef <- length(model$coefficients)
    ndes <- dim(model$model)[1]
    p.names <- names(coef(model))
    nlambda <- length(lambda)
    estimates <- matrix(data = NA, nlambda + 1, ncoef) # +1 because "0" will be added
    theta <- matrix(data = NA, B, ncoef)
    colnames(theta)<- p.names
    theta.all <-  vector(mode = "list", nlambda)
    if (jackknife.estimation != FALSE) {
      var.exp <- list()
#**Heidi**#
# simex::: kann spaeter wieder weggelassen werden
      var.exp[[1]] <- simex:::extract.covmat(model)
#*********#      
    }
    if (asymptotic) {
      psi <- matrix(rep(0, ndes * ncoef), ncol = ncoef, nrow = ndes)
      psi <- resid(model, type = "response") * model$x
      PSI <- psi
      am <- list()
      a <- list()
      xi <- model$x
      dh <- rep(1, ndes)
      if (class(model)[1] == "glm")
        dh <- model$family$mu.eta(model$linear.predictors)
      for (k in 1:ndes)
        a[[k]] <- dh[k] * xi[k, ] %*% t(xi[k, ])
      a.mat <- matrix(unlist(a), nrow = length(a), byrow =TRUE)
      ab <- matrix(colSums(a.mat), nrow = NROW(a[[1]]), byrow = FALSE)
      am[[1]] <- -ab / ndes
      a <- list()
    }
    # assigning the naive estimator
    estimates[1, ] <- model$coefficients
    # The Simulation step
    # outer loop doing the simulations for each lambda
    for (i in 1:length(lambda)) {
      if (jackknife.estimation != FALSE)
        variance.est <- matrix(0, ncol = ncoef, nrow = ncoef)
      if (asymptotic) {
        psi <- matrix(0, ncol = ncoef, nrow = ndes)
        a <- list()
        for(k in 1:ndes)
          a[[k]] <- matrix(0, nrow = ncoef, ncol = ncoef)
      }
      # inner loop, doing the simulations
      for (j in 1:B) {
        SIMEXdata <- model$model              
#**Heidi**#        
#        epsilon <- rnorm(n = NROW(SIMEXdata))
        epsilon <- matrix(rnorm(n=NROW(SIMEXdata)*length(SIMEXvariable))
                          ,ncol= length(SIMEXvariable)
                          ,nrow =  NROW(SIMEXdata))
        # and adding the random error
#        SIMEXdata[, SIMEXvariable] <- SIMEXdata[,SIMEXvariable] +
#          sapply((sqrt(lambda[i]) * measurement.error), "*", epsilon)
        SIMEXdata[, SIMEXvariable] <- SIMEXdata[, SIMEXvariable] + 
          ( sqrt(lambda[i]) *  epsilon * measurement.error  )
#*********#
        # updating the model and calculating the estimate
        model.SIMEX <- update(model, data = data.frame(SIMEXdata))
        theta[j, ] <- model.SIMEX$coefficients
        if (jackknife.estimation != FALSE) {
#**Heidi**#
# simex::: kann spaeter wieder weggelassen werden
          variance.est <- variance.est + simex:::extract.covmat(model.SIMEX)
#*********#
        }
        if (asymptotic) {
          xi <- model.SIMEX$x
          psi <- psi + (resid(model.SIMEX, type = "response") * xi)
          dh <- rep(1, ndes)
          if (class(model)[1] == "glm")
            dh <- model$family$mu.eta(model.SIMEX$linear.predictors)
          for (k in 1:ndes)
            a[[k]] <- a[[k]] - dh[k] * xi[k, ] %*% t(xi[k, ])
        }
      }
      # taking the mean of the estimate -> SIMEX estimate
      estimates[i + 1, ] <- colMeans(theta)
      theta.all[[i]] <- theta
      # Variance estimation via the Jackknife
      if (jackknife.estimation != FALSE) {
        variance.est <- variance.est / B
        s2 <- cov(theta)
        var.exp[[i + 1]]<- variance.est - s2
      }
      if (asymptotic) {
        xiB <- psi / B
        PSI <- cbind(PSI, xiB)
        a.mat <- matrix(unlist(a), nrow = length(a), byrow =TRUE)
        ab <- matrix(colSums(a.mat), nrow = NROW(a[[1]]), byrow = FALSE)
        am[[i + 1]] <- ab / (B * ndes)
      }
    }
    # extrapolation step
    SIMEX.estimate <- vector(mode = "numeric", length = ncoef)
    colnames(estimates) <- p.names
    lambda <- c(0, lambda)
    # fitting the extrapolation function
    switch(fitting.method,
           "quad" = extrapolation <- lm(estimates ~ lambda + I(lambda^2)),
           "line "= extrapolation <- lm(estimates ~ lambda),
           "nonl" = extrapolation <- fit.nls(lambda, p.names, estimates)
    )
    # predicting the SIMEX estimate
    if (fitting.method == "nonl") {
      for (i in 1:length(p.names))
        SIMEX.estimate[i] <- predict(extrapolation[[p.names[i]]],
                                     newdata = data.frame(lambda = -1))
    } else {
      SIMEX.estimate <- predict(extrapolation, newdata = data.frame(lambda = -1))
    }
    ######  Jackknife Estimation
    if (jackknife.estimation != FALSE) {
      variance.jackknife <- matrix(unlist(var.exp), ncol = ncoef^2 , byrow = TRUE)
      switch(jackknife.estimation,
             "quad" = extrapolation.variance <-
               lm(variance.jackknife ~ lambda + I(lambda^2)),
             "line" = extrapolation.variance <- lm(variance.jackknife ~ lambda),
             "nonl" = extrapolation.variance <-
               fit.nls(lambda, 1:NCOL(variance.jackknife), variance.jackknife)
      )
      variance.jackknife2 <- vector("numeric", ncoef^2)
      switch(jackknife.estimation,
             "nonl"= for(i in 1:NCOL(variance.jackknife))
               variance.jackknife2[i] <- predict(extrapolation.variance[[i]],
                                                 newdata = data.frame(lambda = -1)),
             "quad"= variance.jackknife2 <- predict(extrapolation.variance,
                                                    newdata = data.frame(lambda = -1)),
             "line"= variance.jackknife2 <- predict(extrapolation.variance,
                                                    newdata = data.frame(lambda = -1))
      )
      variance.jackknife <- rbind(variance.jackknife2, variance.jackknife)
      variance.jackknife.lambda <- cbind(c(-1,lambda), variance.jackknife)
      variance.jackknife <- matrix(variance.jackknife[1, ], nrow = ncoef,
                                   ncol = ncoef, byrow = TRUE)
      dimnames(variance.jackknife) <- list(p.names, p.names)
    }
    # Asymptotic estimation
    if (asymptotic) {
      c11 <- cov(PSI)
      a11 <- diag.block(am)
      a11.inv <- solve(a11)
      sigma <- a11.inv %*% c11 %*% t(a11.inv)
#**Heidi**#
# simex::: kann spaeter wieder weggelassen werden
      s <- simex:::construct.s(ncoef, lambda, fitting.method, extrapolation)
      d.inv <- solve(s %*% t(s))
      sigma.gamma <- d.inv %*% s %*% sigma %*% t(s) %*% d.inv
      g <- list()
      switch(fitting.method,
             "quad" = g <- c(1, -1, 1),
             "line" = g <- c(1, -1),
             "nonl" = for(i in 1:ncoef)
               g[[i]] <- c(-1, -(coef(extrapolation[[i]])[3] - 1)^-1,
                           coef(extrapolation[[i]])[2] / (coef(extrapolation[[i]])[3] - 1)^2)
      )
      g <- diag.block(g, ncoef)
      variance.asymptotic <- (t(g) %*% sigma.gamma %*% g) / ndes
      dimnames(variance.asymptotic) <- list(p.names, p.names)
    }
    # creating class "simex"
    theta <- matrix(unlist(theta.all), nrow = B)
    theta.all <- list()
    for (i in 1:ncoef)
      theta.all[[p.names[i]]] <-
        data.frame(theta[, seq(i, ncoef * nlambda, by = ncoef)])
    z <- cbind(lambda, estimates)
    z <- rbind(c(-1, SIMEX.estimate), z) # returning the estimated values
    colnames(z) <- c("lambda", names(coef(model)))
    erg <- list(
      coefficients = z[1, -1], # SIMEX corrected coefficients
      SIMEX.estimates = z, # all thetas as a matrix
      lambda = lambda, # vector for the values for lambda
      model = model, # the naive model
      measurement.error = measurement.error, # vector of values of measurement.error
      B = B, # number of Simulations
      extrapolation = extrapolation, # model of the extrapolation
      fitting.method = fitting.method,# which fitting method was used
      SIMEXvariable = SIMEXvariable,
      theta = theta.all,
      call = cl
    )
    class(erg) <- ("simex")
    fitted.values <- predict(erg, newdata = model$model[, -1, drop = FALSE],
                             type = "response")
    erg$fitted.values <- fitted.values
    if (is.factor(model$model[, 1])) {
      erg$residuals <-
        as.numeric(levels(model$model[, 1]))[model$model[, 1]] - fitted.values
    } else {
      erg$residuals <- model$model[, 1] - fitted.values
    }
    if (jackknife.estimation != FALSE) {
      erg$extrapolation.variance <- extrapolation.variance
      erg$variance.jackknife <- variance.jackknife
      erg$variance.jackknife.lambda <- variance.jackknife.lambda
    }
    if (asymptotic) {
      erg$PSI <- PSI
      erg$c11 <- c11
      erg$a11 <- a11
      erg$sigma <- sigma
      erg$sigma.gamma <- sigma.gamma
      erg$g <- g
      erg$s <- s
      erg$variance.asymptotic <- variance.asymptotic
    }
    return(erg)
  }

