
ictregBayes <- function(formula, data = parent.frame(), treat = "treat", J, constrained.single = "full", 
                        constrained.multi = TRUE, fit.start = "lm", n.draws = 10000, burnin = 5000, thin = 0, 
                        delta.start, psi.start, Sigma.start, Phi.start, delta.mu0, psi.mu0, delta.A0, psi.A0, 
                        Sigma.df, Sigma.scale, Phi.df, Phi.scale, delta.tune, psi.tune, gamma.tune, zeta.tune, formula.mixed, group.mixed,
                        verbose = TRUE, sensitive.model = "logit", df = 5,
                        endorse.options,
                        ...){

  ictreg.call <- match.call() 

  ## removed to fix error in package compiling
  ##require(magic)
  
  # set up data frame, with support for standard and modified responses
  mf <- match.call(expand.dots = FALSE)
  
  # make all other call elements null in mf <- NULL in next line
  mf$treat <- mf$J <- mf$constrained.single <- mf$constrained.multi <- mf$verbose <- mf$n.draws <- mf$burnin <- mf$thin <- mf$delta.start <- mf$psi.start <- mf$delta.mu0 <- mf$psi.mu0 <- mf$delta.A0 <- mf$psi.A0 <- mf$delta.tune <- mf$psi.tune <- mf$fit.start <- mf$formula.mixed <- mf$group.mixed <- mf$Sigma.start <- mf$Phi.start <- mf$Sigma.df <- mf$Sigma.scale <- mf$Phi.df <- mf$Phi.scale <- mf$gamma.tune <- mf$zeta.tune <- mf$sensitive.model <- mf$df <- mf$endorse.options <- NULL

    ##mf$ceiling <- mf$floor <- NULL

  mf[[1]] <- as.name("model.frame")
  mf$na.action <- 'na.pass'
  mf <- eval.parent(mf)

  # define design, response data frames
  x.all <- model.matrix.default(attr(mf, "terms"), mf)
  y.all <- model.response(mf)

  # get mixed effects group-level predictors
  mixed <- missing("group.mixed") == FALSE

  endorse <- missing("endorse.options") == FALSE
  ##endorse <- FALSE

  if (mixed == TRUE) {

    if (missing("formula.mixed") == TRUE)
      formula.mixed <- ~ 1

    z.all <- model.matrix(formula.mixed, data)

  }
 
  # list-wise missing deletion
  na.x <- apply(is.na(x.all), 1, sum)
  na.y <- is.na(y.all)

  if (mixed == TRUE) {
    na.z <- apply(is.na(z.all), 1, sum)
    na.cond <- na.x==0 & na.y==0 & na.z==0
  } else {
    na.cond <- na.x==0 & na.y==0
  }

  ## group indicator for mixed effects regression
  if (mixed == TRUE)
    grp <- data[na.cond == TRUE, paste(group.mixed)]
  
  ## treatment indicator for subsetting the dataframe
  t <- data[na.cond == TRUE, paste(treat)]

  if(class(t) == "factor") {
    
    levels(t) <- tolower(levels(t))
    
    if (length(which(levels(t) == "control")) == 1) {
      t <- relevel(t, ref = "control")
    } else {
      warning("Note: using the first level as the control condition, but it is not labeled 'control'.")
    }
        
    condition.labels <- levels(t)
    t <- as.numeric(t) - 1
    treatment.labels <- condition.labels[2:length(condition.labels)]
    control.label <- condition.labels[1]
    
  } else {
    condition.labels <- sort(unique(t))
    treatment.labels <- condition.labels[condition.labels != 0]
    control.label <- 0
  }
  
  ## list wise delete
  y.all <- y.all[na.cond == TRUE]
  x.all <- x.all[na.cond == TRUE, , drop = FALSE]
  if (mixed == TRUE)
    z.all <- z.all[na.cond == TRUE, , drop = FALSE]
  
  ## set up data objects for y and x for each group from user input
  x.treatment <- x.all[t != 0, , drop = FALSE]
  y.treatment <- subset(y.all, t != 0)
  x.control <- x.all[t == 0 , , drop = FALSE]
  y.control <- subset(y.all, t==0)
  if (mixed == TRUE) {
    z.treatment <- z.all[t != 0, , drop = FALSE]
    z.control <- z.all[t == 0 , , drop = FALSE]
  }

  ## J is defined by user for the standard design
  if(missing("J")) {
    J <- max(y.treatment) - 1
  }

  condition.values <- sort(unique(t))
  treatment.values <- 1:length(condition.values[condition.values!=0])

  if(length(treatment.values) > 1) {
    multi <- TRUE
  } else {
    multi <- FALSE
  }
  
  n <- nrow(x.treatment) + nrow(x.control)
  if (mixed == TRUE)
    n.grp <- length(unique(grp))

  coef.names <- colnames(x.all)
  if (mixed == TRUE)
    coef.names.mixed <- colnames(z.all)

  nPar <- ncol(x.all)
  if (mixed == TRUE)
    nPar.mixed <- ncol(z.all)

  intercept.only <- ncol(x.all)==1 & sum(x.all[,1]==1) == n

  if (intercept.only == TRUE) {
    x.vars <- "1"
  } else {
    x.vars <- coef.names[-1]
  }

  if (mixed == TRUE) {
    
    intercept.only.mixed <- ncol(z.all)==1 & sum(z.all[,1]==1) == n
    
    if (intercept.only.mixed == TRUE) {
      x.vars.mixed <- "1"
    } else {
      x.vars.mixed <- coef.names.mixed[-1]
    }
    
  }
  
  logistic <- function(x) exp(x)/(1+exp(x))
  logit <- function(x) return(log(x)-log(1-x))

  if (missing("delta.mu0")) {
    if (multi == FALSE)
      delta.mu0 <- rep(0, nPar)
    else
      delta.mu0 <- rep(0, nPar + (1-constrained.multi))
  }
  if (missing("psi.mu0")) {
    if (multi == FALSE)
      psi.mu0 <- rep(0, nPar + (constrained.single == "intercept"))
    else
      psi.mu0 <- rep(0, nPar)
  }
  
  if (missing("delta.A0")) {
    if (multi == FALSE)
      delta.A0 <- .01 * diag(nPar)
    else
      delta.A0 <- .01 * diag(nPar + (1-constrained.multi))
  }
  if (missing("psi.A0")) {
    if (multi == FALSE) {
      psi.A0 <- .01 * diag(nPar + (constrained.single == "intercept"))
    } else {
      psi.A0 <- 0.01*diag(nPar)
    }
  }
  if (missing("delta.tune") & endorse == FALSE) {
    stop("The Metropolis tuning input object delta.tune is required.")

  }
  if (missing("psi.tune")) {
    stop("The Metropolis tuning input object psi.tune is required.") 

  }

  if (mixed == TRUE) {
    if (missing("Sigma.start")) {
      Sigma.start <- matrix(0.005, ncol = ncol(z.all), nrow = ncol(z.all))
      diag(Sigma.start) <- 0.01
    }
    if (missing("Phi.start")) {
      Phi.start <- matrix(0.01, ncol = ncol(z.all), nrow = ncol(z.all))
      diag(Phi.start) <- 0.025
    }
    if (missing("Sigma.df")) {
      Sigma.df <- ncol(z.all)+1
    }
    if (missing("Sigma.scale")) {
      Sigma.scale <- diag(ncol(z.all))*0.1
    }
    if (missing("Phi.df")) {
      Phi.df <- ncol(z.all)+1
    }
    if (missing("Phi.scale")) {
      Phi.scale <- diag(ncol(z.all))*0.1
    }
    if (missing("gamma.tune")) {
      gamma.tune <- rep(0.001, n.grp)
    }
    if (missing("zeta.tune")) {
      zeta.tune <- rep(0.001, n.grp)
    }
  }
    
  if (multi == FALSE) {

    mle.constrained <- constrained.single == "full"

    if ((missing(delta.start) & endorse == FALSE) | missing(psi.start))
      ictreg.fit <- ictreg(formula, data = data, treat = treat, J = J, method = fit.start, constrained = mle.constrained)

    if (missing(delta.start) & endorse == FALSE)
      delta.start <- ictreg.fit$par.treat
    if (missing(psi.start))
      psi.start <- ictreg.fit$par.control

    if (constrained.single == "intercept")
      psi.start <- c(psi.start, 0)

    constrained.pass <- ifelse(constrained.single == "full", 0, ifelse(constrained.single == "none", 1, 2))

    if (endorse == FALSE)  {
      if (mixed == FALSE) {
        
        ## do standard single sensitive item design
    
        ictregBayes.fit <- ictregBayes.fit(Y = y.all, treat = t, X = x.all, J = J, constrained = constrained.pass,
                                           n.draws = n.draws, burnin = burnin, thin = thin, verbose = verbose,
                                           delta.start = delta.start, psi.start = psi.start,
                                           delta.mu0 = delta.mu0, psi.mu0 = psi.mu0, delta.A0 = delta.A0, psi.A0 = psi.A0,
                                           delta.tune = delta.tune, psi.tune = psi.tune,
                                           ##ceiling = ceiling, floor = floor,
                                           0, 0,
                                           robit = (sensitive.model == "robit"),
                                           probit = (sensitive.model == "probit"), df = df,
                                           ...)
        
        if (constrained.single == "full" | constrained.single == "intercept")
          M <- 1
        else
          M <- 2

        if(sensitive.model != "logit")
          M <- M - 1
        
        fit <- c()
        fit$delta <- ictregBayes.fit[, 1:nPar, drop = FALSE]
        fit$psi <- ictregBayes.fit[, (nPar+1):(ncol(ictregBayes.fit)-M-1), drop = FALSE]
        if(sensitive.model == "logit")
          fit$delta.accept <- ictregBayes.fit[nrow(ictregBayes.fit), ncol(ictregBayes.fit)-M, drop = FALSE]
        fit$psi.accept <- ictregBayes.fit[nrow(ictregBayes.fit),
                                          (ncol(ictregBayes.fit)-M+
                                           as.numeric(sensitive.model == "logit")):ncol(ictregBayes.fit), drop = FALSE]
        
        rownames(fit$delta) <- rownames(fit$psi) <-
          seq(from = n.draws - burnin + 1 + thin + 1,
              to = n.draws - burnin + thin + 1 + floor((n.draws - burnin)/(thin + 1)) * (thin + 1), by = thin + 1)                      
        fit$delta <- mcmc(data = fit$delta, start = 1, thin = 1, end = nrow(fit$delta))
        fit$psi <- mcmc(data = fit$psi, start = 1, thin = 1, end = nrow(fit$psi))
        
      } else {
        
        ## do mixed effects single sensitive item design
        
        fit <- ictregBayesMixed.fit(Y = y.all, treat = t, X = x.all, Z = z.all, J = J,
                                    grp = as.numeric(grp), constrained = constrained.pass,
                                    n.draws = n.draws, burnin = burnin, thin = thin,  verbose = verbose,
                                    delta.start = delta.start, psi.start = psi.start,
                                    Sigma.start = Sigma.start, Phi.start = Phi.start,
                                    Sigma.df = Sigma.df, Sigma.scale = Sigma.scale, Phi.df = Phi.df, Phi.scale = Phi.scale,
                                    delta.mu0 = delta.mu0, psi.mu0 = psi.mu0, delta.A0 = delta.A0,
                                    psi.A0 = psi.A0, delta.tune = delta.tune, psi.tune = psi.tune,
                                    gamma.tune = gamma.tune, zeta.tune = zeta.tune,
                                    ##ceiling = ceiling, floor = floor,
                                    ...)
                
        rownames(fit$delta) <- rownames(fit$psi) <- names(fit$Sigma) <-
          names(fit$Phi) <- rownames(fit$gamma) <- rownames(fit$zeta) <-
            seq(from = n.draws - burnin + 1 + thin + 1,
                to = n.draws - burnin + thin + 1 + floor((n.draws - burnin)/(thin + 1)) * (thin + 1), by = thin + 1)
        
        end <- nrow(fit$delta)
        fit$delta <- mcmc(data = fit$delta, start = 1, thin = 1, end = end)
        fit$psi <- mcmc(data = fit$psi, start = 1, thin = 1, end = end)
        fit$Sigma <- mcmc(data = fit$Sigma, start = 1, thin = 1, end = end)
        fit$Phi <- mcmc(data = fit$Phi, start = 1, thin = 1, end = end)
        fit$gamma <- mcmc(data = fit$gamma, start = 1, thin = 1, end = end)
        fit$zeta <- mcmc(data = fit$zeta, start = 1, thin = 1, end = end)
        
      }
      
      constrained.output <- constrained.single
    } else {
      
      ## do endorse

      fit <- eval(as.call(c(list(ictregBayesEndorse,Y.list = y.all, treat.list = t, J.list = J,
                                 ##ceiling = ceiling, floor = floor, 
                                 ##constrained = constrained.multi,
                                 verbose = verbose,
                                 psi.start = psi.start,
                                 psi.mu0 = psi.mu0, 
                                 psi.A0 = psi.A0,
                                 psi.tune = psi.tune,
                                 robit = (sensitive.model == "robit"), probit = (sensitive.model == "probit"),
                                 df = df), endorse.options)))

      constrained.output <- TRUE

     }
    
  } else {
    
    ## multi code
    
    if (missing(delta.start) | missing(psi.start)) {
      if (constrained.multi == T)
        ictreg.fit <- ictreg(formula, data = data, treat = treat, J = J, method = "ml", multi.condition = "none")
      else
        ictreg.fit <- ictreg(formula, data = data, treat = treat, J = J, method = "ml", multi.condition = "level")

      par.treat <- list()
      for (m in 1:length(treatment.values))
        par.treat[[as.character(m)]] <- ictreg.fit$par.treat[[m]]
    }
     
    if (missing(delta.start))
      delta.start <- par.treat
    if (missing(psi.start))
      psi.start <- ictreg.fit$par.control

    ##if (length(ceiling) == 1)
    ##  ceiling <- rep(ceiling, length(treatment.values))
    ##if (length(floor) == 1)
    ##  floor <- rep(floor, length(treatment.values))

    if (endorse == FALSE) {
      
      if (mixed == FALSE) {
        
        ictregBayesMulti.fit <- ictregBayesMulti.fit(Y = y.all, treat = t, X = x.all, J = J,
                                                     constrained = constrained.multi,
                                                     n.draws = n.draws, burnin = burnin, thin = thin,  verbose = verbose,
                                                     delta.start = delta.start, psi.start = psi.start,
                                                     delta.mu0 = delta.mu0, psi.mu0 = psi.mu0, delta.A0 = delta.A0,
                                                     psi.A0 = psi.A0, delta.tune = delta.tune, psi.tune = psi.tune,
                                                     ##ceiling = ceiling, floor = floor,
                                                     ...)
        
       
        iter.names <-  seq(from = n.draws - burnin + 1 + thin + 1,
                           to = n.draws - burnin + thin + 1 + floor((n.draws - burnin)/(thin + 1)) * (thin + 1),
                           by = thin + 1)
        
        fit <- c()
        fit$delta <- list()
        for (m in 1:length(treatment.labels)) {
          fit$delta[[m]] <- ictregBayesMulti.fit[, ((m-1)*nPar + 1):(m*nPar), drop = FALSE]
          rownames(fit$delta[[m]]) <- iter.names
          fit$delta[[m]] <- mcmc(data = fit$delta[[m]], start = 1, thin = 1, end = nrow(fit$delta[[m]]))
        }
        fit$psi <- ictregBayesMulti.fit[, (length(treatment.labels)*nPar + 1):
                                        (ncol(ictregBayesMulti.fit)-length(treatment.labels)-1), drop = FALSE]
        if(sensitive.model == "logit")
          fit$delta.accept <- ictregBayesMulti.fit[nrow(ictregBayesMulti.fit),
                                                   (ncol(ictregBayesMulti.fit)-length(treatment.labels)):
                                                   (ncol(ictregBayesMulti.fit)-1), drop = FALSE]
        rownames(fit$psi) <- iter.names
        
        fit$psi <- mcmc(data = fit$psi, start = 1, thin = 1, end = nrow(fit$psi))
        
        fit$psi.accept <- ictregBayesMulti.fit[nrow(ictregBayesMulti.fit), ncol(ictregBayesMulti.fit), drop = FALSE]
        
      } else {
        
        ## mixed multi model
        
        fit <- ictregBayesMultiMixed.fit(Y = y.all, treat = t, X = x.all, Z = z.all, J = J,
                                         grp = as.numeric(grp), constrained = constrained.multi,
                                         n.draws = n.draws, burnin = burnin, thin = thin,  verbose = verbose,
                                         delta.start = delta.start, psi.start = psi.start,
                                         Sigma.start = Sigma.start, Phi.start = Phi.start,
                                         Sigma.df = Sigma.df, Sigma.scale = Sigma.scale, Phi.df = Phi.df,
                                         Phi.scale = Phi.scale,
                                         delta.mu0 = delta.mu0, psi.mu0 = psi.mu0, delta.A0 = delta.A0,
                                         psi.A0 = psi.A0, delta.tune = delta.tune, psi.tune = psi.tune,
                                         gamma.tune = gamma.tune, zeta.tune = zeta.tune,
                                         ##ceiling = ceiling, floor = floor,
                                         ...)
        
        rownames(fit$delta) <- rownames(fit$psi) <- names(fit$Sigma) <-
          names(fit$Phi) <- rownames(fit$gamma) <- rownames(fit$zeta) <-
            seq(from = n.draws - burnin + 1 + thin + 1,
                to = n.draws - burnin + thin + 1 + floor((n.draws - burnin)/(thin + 1)) * (thin + 1),
                by = thin + 1)
        
        delta <- list()
        for (m in 1:length(treatment.labels)) {
          delta[[m]] <- mcmc(data = fit$delta[, ((m-1)*nPar + 1):(m*nPar)],
                             start = 1, thin = 1, end = nrow(fit$delta))
        }
        
        fit$delta <- delta
        fit$psi <- mcmc(data = fit$psi, start = 1, thin = 1, end = nrow(fit$psi))
        fit$Sigma <- mcmc(data = fit$Sigma, start = 1, thin = 1, end = nrow(fit$Sigma))
        fit$Phi <- mcmc(data = fit$Phi, start = 1, thin = 1, end = length(fit$Phi))
        fit$gamma <- mcmc(data = fit$gamma, start = 1, thin = 1, end = nrow(fit$gamma))
        fit$zeta <- mcmc(data = fit$zeta, start = 1, thin = 1, end = nrow(fit$zeta))
        
      } ## end mixed model
      
      constrained.output <- constrained.multi
      
    } else {

      ## start combined model (endorse)
      
      fit <- eval(as.call(c(list(ictregBayesEndorse,Y.list = y.all, treat.list = t, J.list = J,
                                ##ceiling = ceiling, floor = floor, 
                                ##constrained = constrained.multi,
                                verbose = verbose,
                                psi.start = psi.start,
                                psi.mu0 = psi.mu0, 
                                psi.A0 = psi.A0,
                                psi.tune = psi.tune,
                                robit = (sensitive.model == "robit"), probit = (sensitive.model == "probit"),
                                df = df), endorse.options)))

      constrained.output <- TRUE
      
      
    }
  }
  
  fit$x <- x.all
  fit$multi <- multi
  fit$mixed <- mixed
  fit$constrained <- constrained.output
  if(endorse == FALSE)
    fit$delta.start <- delta.start
  fit$psi.start <- psi.start
  if (endorse == FALSE)
    fit$delta.mu0 <- delta.mu0
  fit$psi.mu0 <- psi.mu0
  if (endorse == FALSE)
    fit$delta.A0 <- delta.A0
  fit$psi.A0 <- psi.A0
  if (endorse == FALSE)
    fit$delta.tune <- delta.tune
  fit$psi.tune <- psi.tune
  fit$J <- J
  fit$treat.labels <- treatment.labels
  fit$control.label <- control.label
  fit$coef.names <- coef.names
  fit$sensitive.model <- sensitive.model
  fit$call <- match.call()
  
  class(fit) <- "ictregBayes"
  
  return(fit)
  
}

coef.ictregBayes <- function(object, ranef = FALSE, ...) {

  if (object$multi == TRUE) {
    delta.coef <- list()
    for(m in 1:length(object$delta)) {
      delta.coef[[object$treat.labels[[m]]]] <- apply(object$delta[[m]], 2, mean)
      names(delta.coef[[object$treat.labels[[m]]]]) <- object$coef.names
    }
  } else {
    delta.coef <- apply(object$delta, 2, mean)
    names(delta.coef) <- object$coef.names
  }

  psi.coef <- apply(object$psi, 2, mean)
  names(psi.coef) <- object$coef.names
 
  if (ranef == FALSE) {
    return.object <- list(delta = delta.coef, psi = psi.coef)
  } else {
    gamma.coef <- apply(object$gamma, 2, mean)
    names(gamma.coef) <- object$coef.names
    
    zeta.coef <- apply(object$zeta, 2, mean)
    names(zeta.coef) <- object$coef.names
    
    return.object <- list(delta = delta.coef, psi = psi.coef,
                          ranef.gamma = gamma.coef, ranef.zeta = zeta.coef)
  }
  
  return.object

}

coef.ictregBayes.list <- function(object, ranef = FALSE, ...) {

  if (object$multi == TRUE) {
    delta.list <- list()
    for (m in 1:length(object$treat.labels))
      delta.list[[m]] <- as.mcmc(do.call(rbind, object$delta[[m]]))
    object$delta <- delta.list
  } else {
    object$delta <- as.mcmc(do.call(rbind, as.list(object$delta)))
  }      
  
  object$psi <- as.mcmc(do.call(rbind, as.list(object$psi)))

  if (ranef == TRUE) {
    object$gamma <- as.mcmc(do.call(rbind, as.list(object$gamma)))
    object$zeta <- as.mcmc(do.call(rbind, as.list(object$zeta)))
  }
  
  class(object) <- "ictregBayes"

  coef(object, ranef = ranef, ... = ...)

}

sd.ictregBayes <- function(object, ranef = FALSE, ...) {
  
  if (object$multi == TRUE) {
    delta.coef <- list()
    for(m in 1:length(object$delta)) {
      delta.coef[[object$treat.labels[[m]]]] <- apply(object$delta[[m]], 2, sd)
      names(delta.coef[[object$treat.labels[[m]]]]) <- object$coef.names
    }
  } else {
    delta.coef <- apply(object$delta, 2, sd)
    names(delta.coef) <- object$coef.names
  }

  psi.coef <- apply(object$psi, 2, sd)
  names(psi.coef) <- object$coef.names
 
  if (ranef == FALSE) {
    return.object <- list(delta = delta.coef, psi = psi.coef)
  } else {
    gamma.coef <- apply(object$gamma, 2, sd)
    names(gamma.coef) <- object$coef.names
    
    zeta.coef <- apply(object$zeta, 2, sd)
    names(zeta.coef) <- object$coef.names
    
    return.object <- list(delta = delta.coef, psi = psi.coef,
                          ranef.gamma = gamma.coef, ranef.zeta = zeta.coef)
  }
  
  return.object


}

sd.ictregBayes.list <- function(object, ranef = FALSE, ...) {

  if (object$multi == TRUE) {
    delta.list <- list()
    for (m in 1:length(object$treat.labels))
      delta.list[[m]] <- as.mcmc(do.call(rbind, object$delta[[m]]))
    object$delta <- delta.list
  } else {
    object$delta <- as.mcmc(do.call(rbind, as.list(object$delta)))
  }      
  
  object$psi <- as.mcmc(do.call(rbind, as.list(object$psi)))

  if (ranef == TRUE) {
    object$gamma <- as.mcmc(do.call(rbind, as.list(object$gamma)))
    object$zeta <- as.mcmc(do.call(rbind, as.list(object$zeta)))
  }
  
  class(object) <- "ictregBayes"

  sd.ictregBayes(object, ranef = ranef, ... = ...)

}

vcov.ictregBayes <- function(object, ranef = FALSE, ...) {
 
  if (object$multi == TRUE)
    delta.draws <- cbind(do.call(cbind, object$delta))
  else
    delta.draws <- object$delta

  
  if (ranef == TRUE)
    cov(cbind(delta.draws, object$psi, object$gamma, object$zeta, object$Sigma, object$Phi))
  else
    cov(cbind(delta.draws, object$psi))
}

vcov.ictregBayes.list <- function(object, ranef = FALSE, ...) {
  
  if (object$multi == TRUE) {
    delta.list <- list()
    for (m in 1:length(object$treat.labels))
      delta.list[[m]] <- as.mcmc(do.call(rbind, as.list(object$delta[[m]])))
    object$delta <- delta.list
  } else {
    object$delta <- as.mcmc(do.call(rbind, as.list(object$delta)))
  }      
  
  object$psi <- as.mcmc(do.call(rbind, as.list(object$psi)))

  if (ranef == TRUE) {
    object$gamma <- as.mcmc(do.call(rbind, as.list(object$gamma)))
    object$zeta <- as.mcmc(do.call(rbind, as.list(object$zeta)))
    object$Sigma <- as.mcmc(do.call(rbind, as.list(object$Sigma)))
    object$Phi <- as.mcmc(do.call(rbind, as.list(object$Phi)))
  }
  
  class(object) <- "ictregBayes"

  vcov(object, ranef = ranef, ... = ...)
  
}

print.ictregBayes <- print.ictregBayes.list <- function(x, ...) {
  
  cat("\nItem Count Technique Bayesian Regression \n\nCall: ")
  
  dput(x$call)
  
  cat("\nCoefficient estimates\n")

  print(coef(x))

  treat.print <- c()
  for (i in 1:length(x$treat.labels)) {
    treat.print <- c(treat.print, "'", x$treat.labels[i], "'", sep = "")
    if (i != length(x$treat.labels))
      treat.print <- c(treat.print, " and ")
  }
  
  cat("\nNumber of control items J set to ", x$J, ". Treatment groups were indicated by ", sep = "")
  cat(treat.print, sep ="")
  cat(" and the control group by '", x$control.label, "'.\n\n", sep = "")
    
  invisible(x)
  
}


as.list.ictregBayes <- function(...) {
  
  x <- list(...)

  if (x[[1]]$multi == TRUE) {
    delta.list <- list()
    for (m in 1:length(x[[1]]$treat.labels)) {
      delta.indiv.list <- list()
      for (i in 1:length(x))
        delta.indiv.list[[i]] <- as.mcmc(as.matrix(x[[i]]$delta[[m]]))
      
      delta.list[[m]] <- as.mcmc.list(delta.indiv.list)
    }
  } else {
    delta.list <- list()
    for (i in 1:length(x))
      delta.list[[i]] <- x[[i]]$delta
      
    delta.list <- as.mcmc.list(delta.list)
  }
  
  psi.list <- list()
  for (i in 1:length(x))
    psi.list[[i]] <- x[[i]]$psi

  psi.list <- as.mcmc.list(psi.list)

  if (x[[1]]$mixed == TRUE) {
    gamma.list <- list()
    for (i in 1:length(x))
      gamma.list[[i]] <- x[[i]]$gamma
    
    gamma.list <- as.mcmc.list(gamma.list)
    
    zeta.list <- list()
    for (i in 1:length(x))
      zeta.list[[i]] <- x[[i]]$zeta
    
    zeta.list <- as.mcmc.list(zeta.list)

    Phi.list <- list()
    for (i in 1:length(x))
      Phi.list[[i]] <- x[[i]]$Phi
    
    Phi.list <- as.mcmc.list(Phi.list)

    Sigma.list <- list()
    for (i in 1:length(x))
      Sigma.list[[i]] <- x[[i]]$Sigma
    
    Sigma.list <- as.mcmc.list(Sigma.list)
  }
  
  return.object <- x[[1]]
  return.object$delta <- delta.list
  return.object$psi <- psi.list
  if (x[[1]]$mixed == TRUE) {
    return.object$gamma <- gamma.list
    return.object$zeta <- zeta.list
    return.object$Phi <- Phi.list
    return.object$Sigma <- Sigma.list
  }
  
  class(return.object) <- "ictregBayes.list"

  return.object
  
}

summary.ictregBayes <- function(object, ...) {
  structure(object, class = c("summary.ictregBayes", class(object)))
}

print.summary.ictregBayes <- function(x, ...) {
  
  cat("\nItem Count Technique Bayesian Regression \n\nCall: ")
  
  dput(x$call)

  if (x$multi == TRUE) {
    for (k in 1:length(x$treat.labels)) {
      cat(paste("\nSensitive item (", x$treat.labels[k], ")", "\n", sep = ""))
      print(matrix(c(round(cbind(coef(x)$delta[[k]], sd.ictregBayes(x)$delta[[k]]),5)),
                   nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
                   dimnames = list(x$coef.names, c("Est.", "S.E."))))
      if(x$sensitive.model == "logit")
        cat("\nMetropolis acceptance ratio:", round(x$delta.accept[[k]], 3), "\n")    
    }
  } else {
    
    cat("\nSensitive item \n")
    print(matrix(c(round(cbind(coef(x)$delta, sd.ictregBayes(x)$delta),5)), nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
                 dimnames = list(x$coef.names, c("Est.", "S.E."))))
    if(x$sensitive.model == "logit")
      cat("\nMetropolis acceptance ratio:", round(x$delta.accept, 3), "\n")
    
  }
    
  cat("\nControl items \n")
  print(matrix(c(round(cbind(coef(x)$psi, sd.ictregBayes(x)$psi),5)), nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
                  dimnames = list(x$coef.names, c("Est.", "S.E."))))
  cat("\nMetropolis acceptance ratio:", round(x$psi.accept, 3), "\n")

  treat.print <- c()
  for (i in 1:length(x$treat.labels)) {
    treat.print <- c(treat.print, "'", x$treat.labels[i], "'", sep = "")
    if (i != length(x$treat.labels))
      treat.print <- c(treat.print, " and ")
  }
  
  cat("\nNumber of control items J set to ", x$J, ". Treatment groups were indicated by ", sep = "")
  cat(treat.print, sep ="")
  cat(" and the control group by '", x$control.label, "'.\n\n", sep = "")
    
  invisible(x)
  
}

summary.ictregBayes.list <- function(object, ...) {
  structure(object, class = c("summary.ictregBayes.list", class(object)))
}

print.summary.ictregBayes.list <- function(x, ...) {
  
  cat("\nItem Count Technique Bayesian Regression \n\nCall: ")
  
  dput(x$call)

  cat("\nSummary from",length(x$psi),"chains\n\n")

  for (k in 1:length(x$treat.labels)) {
    cat(paste("\nSensitive item (", x$treat.labels[k], ")", "\n", sep = ""))
      print(matrix(c(round(cbind(coef(x)$delta[[k]], sd.ictregBayes.list(x)$delta[[k]]),5)),
                   nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
                   dimnames = list(x$coef.names, c("Est.", "S.E."))))
    if(x$sensitive.model == "logit")
      cat("\nMetropolis acceptance ratio:", round(x$delta.accept[[k]], 3), "\n")
    cat("\nGelman-Rubin statistics:\n")
        
    if(x$multi == TRUE){
      gelmanrubin <- round(gelman.diag(x$delta[[k]])$psrf[,1],4)
    } else {
      gelmanrubin <- round(gelman.diag(x$delta)$psrf[,1],4)
    }
    names(gelmanrubin) <- x$coef.names
    
    print(gelmanrubin) 
  }
  
  cat("\nControl items \n")
  print(matrix(c(round(cbind(coef(x)$psi, sd.ictregBayes.list(x)$psi),5)), nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
                  dimnames = list(x$coef.names, c("Est.", "S.E."))))
  cat("\nMetropolis acceptance ratio:", round(x$psi.accept, 3), "\n")

  cat("\nGelman-Rubin statistics:\n")

  gelmanrubin <- round(gelman.diag(x$psi)$psrf[,1],4)
  names(gelmanrubin) <- x$coef.names

  print(gelmanrubin)

  treat.print <- c()
  for (i in 1:length(x$treat.labels)) {
    treat.print <- c(treat.print, "'", x$treat.labels[i], "'", sep = "")
    if (i != length(x$treat.labels))
      treat.print <- c(treat.print, " and ")
  }

  cat("\nNumber of control items J set to ", x$J, ". Treatment groups were indicated by ", sep = "")
  cat(treat.print, sep ="")
  cat(" and the control group by '", x$control.label, "'.\n\n", sep = "")
    
  invisible(x)
  
}

predict.ictregBayes.list <- function(object, ...) {
  
  if (object$multi == TRUE) {
    delta.list <- list()
    for (m in 1:length(object$treat.labels))
      delta.list[[m]] <- as.mcmc(do.call(rbind, object$delta[[m]]))
    object$delta <- delta.list
  } else {
    object$delta <- as.mcmc(do.call(rbind, as.list(object$delta)))
  }      
  
  object$psi <- as.mcmc(do.call(rbind, as.list(object$psi)))

  class(object) <- "ictregBayes"

  predict(object, ... = ...)

}

predict.ictregBayes <- function(object, newdata, newdata.diff, direct.glm, se.fit = FALSE,
                                interval = c("none","confidence"), level = .95, sensitive.item, ...){

  n.draws <- nrow(object$psi)
  
  if(missing(interval)) interval <- "none"

  if(missing(sensitive.item)) {
    sensitive.item <- 1
    if(object$multi==TRUE)
      warning("Using the first sensitive item for predictions. Change with the sensitive.item option.")
  }

  nPar <- length(object$coef.names)
    
  diff <- missing(newdata.diff) == FALSE
  direct <- missing(direct.glm) == FALSE
  if (diff == TRUE & direct == TRUE)
    stop("The difference and direct comparison options cannot be used simultaneously. Try just one.")
  
  logistic <- function(object) exp(object)/(1+exp(object))
  
  if(missing(newdata)) {
    xvar <- object$x
  } else { 
    if(nrow(newdata)==0)
      stop("No data in the provided data frame.")
    xvar <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))), newdata)
  }

  if (class(object$delta) != "list")
    draws.list <- object$delta
  else
    draws.list <- object$delta[[sensitive.item]]
  
  if (direct == TRUE) {
    beta.direct <- coef(direct.glm)
    vcov.direct <- vcov(direct.glm)
    
    xvar.direct <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))),
                                direct.glm$data)

    if (ncol(xvar.direct) != nPar)
      stop("Different number of covariates in direct and list regressions.")

    draws.direct <- mvrnorm(n = n.draws, mu = beta.direct, Sigma = vcov.direct)

    pred.direct.mean <- pred.direct.diff.mean <- rep(NA, n.draws)

  }

  if (diff == TRUE) {
    xvar.diff <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))),
                              newdata.diff)
    pred.diff.mean <- pred.diff.diff.mean <- rep(NA, n.draws)
  }
  
  pred.list.mean <- rep(NA, n.draws)
  
  for (d in 1:n.draws) {
        
    pred.list <- logistic(as.matrix(xvar) %*% as.matrix(draws.list[d, ]))
    pred.list.mean[d] <- mean(pred.list)

    if (diff == TRUE) {
      pred.diff <- logistic(as.matrix(xvar.diff) %*% as.matrix(draws.list[d, ]))
      pred.diff.mean[d] <- mean(pred.list)
      pred.diff.diff.mean[d] <- pred.list.mean[d] - pred.diff.mean[d]
    }
    
    if (direct == TRUE) {
      pred.direct <- logistic(as.matrix(xvar.direct) %*% as.matrix(draws.direct[d,]))
      pred.direct.mean[d] <- mean(pred.direct)
      pred.direct.diff.mean[d] <- pred.list.mean[d] - pred.direct.mean[d]
    }
    
  }

  critical.value <- qt(1-(1-level)/2, df = nrow(xvar))
  
  est.list <- mean(pred.list.mean)
  se.list <- sd(pred.list.mean)
  ci.upper.list <- est.list + critical.value * se.list
  ci.lower.list <- est.list - critical.value * se.list

  if (direct == TRUE) {
    est.direct <- mean(pred.direct.mean)
    se.direct <- sd(pred.direct.mean)
    est.direct.diff <- mean(pred.direct.diff.mean)
    se.direct.diff <- sd(pred.direct.diff.mean)
    ci.upper.direct <- est.direct + critical.value * se.direct
    ci.lower.direct <- est.direct - critical.value * se.direct
    ci.upper.direct.diff <- est.direct.diff + critical.value * se.direct.diff
    ci.lower.direct.diff <- est.direct.diff - critical.value * se.direct.diff

    return.object <- list(fit = rbind(est.list, est.direct, est.direct.diff))
      
    if ( interval == "confidence") {
      return.object <- as.data.frame(rbind(c(est.list, ci.lower.list, ci.upper.list),
                                      c(est.direct, ci.lower.direct, ci.upper.direct),
                                      c(est.direct.diff, ci.lower.direct.diff, ci.upper.direct.diff)))
      colnames(return.object) <- c("fit","lwr","upr")
      rownames(return.object) <- c("List", "Direct", "Difference (list - direct)")
    }
      
    if (se.fit == T)
      return.object$se.fit <- c(se.list, se.direct, se.direct.diff)

    attr(return.object, "concat") <- TRUE
    
  }

  if (diff == TRUE) {
    est.diff <- mean(pred.diff.mean)
    se.diff <- sd(pred.diff.mean)
    est.diff.diff <- mean(pred.diff.diff.mean)
    se.diff.diff <- sd(pred.diff.diff.mean)
    ci.upper.diff <- est.diff + critical.value * se.diff
    ci.lower.diff <- est.diff - critical.value * se.diff
    ci.upper.diff.diff <- est.diff.diff + critical.value * se.diff.diff
    ci.lower.diff.diff <- est.diff.diff - critical.value * se.diff.diff

    return.object <- list(fit = rbind(est.list, est.diff, est.diff.diff))
      
    if ( interval == "confidence") {
      return.object <- as.data.frame(rbind(c(est.list, ci.lower.list, ci.upper.list),
                                      c(est.diff, ci.lower.diff, ci.upper.diff),
                                      c(est.diff.diff, ci.lower.diff.diff, ci.upper.diff.diff)))
      colnames(return.object) <- c("fit","lwr","upr")
      rownames(return.object) <- c("Dataset 1", "Dataset 2", "Difference (1 - 2)")
    }
      
    if (se.fit == T)
      return.object$se.fit <- c(se.list, se.diff, se.diff.diff)

    attr(return.object, "concat") <- TRUE
    
  }

  if (diff == FALSE & direct == FALSE) {

    return.object <- list(fit = est.list)
      
    if ( interval == "confidence") {
      return.object <- as.data.frame(rbind(c(est.list, ci.lower.list, ci.upper.list)))
      names(return.object) <- c("fit","lwr","upr")
    }
      
    if (se.fit == T)
      return.object$se.fit <- c(se.list)
  }

  class(return.object) <- "predict.ictreg"
  
  return.object
  
}
 
ictregBayes.fit <- function(Y, treat, X, J, constrained,
                            ceiling, floor,
                            n.draws, burnin, thin, verbose,
                            delta.start, psi.start, delta.mu0,
                            psi.mu0, delta.A0, psi.A0, delta.tune,
                            psi.tune, robit, probit, df) {

  levels.treat <- as.numeric(names(table(treat)))
  
  ## "treat" variable should be a factor variable here where the base level is control
  tmax <- length(levels.treat) - 1
  
  n <- length(Y)
  k <- ncol(X)
  n.par <- 2*(k + 1)
  keep <- thin + 1
  
  if (constrained == 1) {
    if (is.list(psi.start)) {
      psi.start <- c(psi.start$psi0, psi.start$psi1)
    } else {
      psi.start <- rep(psi.start, 2)
    }      
    if (is.list(psi.mu0)) {
      psi.mu0 <- c(psi.mu0$psi0, psi.mu0$psi1)
    } else {
      psi.mu0 <- rep(psi.mu0, 2)
    }
    if (is.list(psi.A0)) {
      psi.A0 <- c(as.double(psi.A0$psi0), as.double(psi.A0$psi1))
    } else {
      psi.A0 <- rep(as.double(psi.A0), 2)
    }
    if (is.list(psi.tune)) {
      psi.tune <- c(as.double(psi.tune$psi0), as.double(psi.tune$psi1))
    } else {
      psi.tune <- rep(as.double(psi.tune), 2)
    }
    n.par <- 3*(k + 1)
  } else if (constrained == 2) {
    n.par <- n.par + 1
  }
  
  if (probit) {
    
    ## no need for deltaCounter
    n.par <- n.par - tmax
    ## no unconstrained model for now
    if (constrained == 1)
      stop("only constrained model is allowed for probit regression")

    res <- .C("ictregBinomMultiProbit", as.integer(Y), as.integer(J),
              as.integer(n), as.integer(n.draws), as.integer(treat), as.integer(1), 
              as.double(X), as.double(delta.start),
              as.double(psi.start), as.integer(k),
              as.double(delta.mu0),
              as.double(psi.mu0), as.double(delta.A0),
              as.double(psi.A0), 
              as.double(psi.tune), 
              0, 0,
              as.integer(burnin), as.integer(keep), as.integer(verbose),
              allresults = double(n.par*floor((n.draws - burnin)/keep)),
              PACKAGE = "list")$allresults
    
  } else if (robit) {
    
    ## no need for deltaCounter
    n.par <- n.par - tmax
    ## no unconstrained model for now
    if (constrained == 1)
      stop("only constrained model is allowed for robit regression")

    res <- .C("ictregBinomMultiRobit", as.integer(Y), as.integer(J),
              as.integer(n), as.integer(n.draws), as.integer(treat), as.integer(1), 
              as.double(X), as.double(delta.start),
              as.double(psi.start), as.integer(k),
              as.double(delta.mu0),
              as.double(psi.mu0), as.double(delta.A0),
              as.double(psi.A0), 
              as.double(psi.tune), as.integer(df),
              0, 0,
              as.integer(burnin), as.integer(keep), as.integer(verbose),
              allresults = double(n.par*floor((n.draws - burnin)/keep)),
              PACKAGE = "list")$allresults
    
  } else {
    res <- .C("ictregBinom", as.integer(Y), as.integer(J), as.integer(n),
              as.integer(n.draws), as.integer(treat), as.double(X), 
              as.double(delta.start), as.double(psi.start), as.integer(k),
              as.double(delta.mu0), as.double(psi.mu0), as.double(delta.A0),
              as.double(psi.A0), as.double(delta.tune), as.double(psi.tune),
              as.integer(constrained),
              ##as.integer(ceiling), as.integer(floor),
              0, 0,
              as.integer(burnin), as.integer(keep), as.integer(verbose),
              allresults = double(n.par*floor((n.draws - burnin)/keep)),
              PACKAGE = "list")$allresults
  }
  res <- matrix(res, byrow = TRUE, ncol = n.par)
  
  class(res) <- "ictregBayes"
  return(res)
            
}

ictregBayesMulti.fit <- function(Y, treat, X, J, constrained,
                                 ##ceiling, floor,
                                 n.draws, burnin, thin, verbose,
                                 delta.start, psi.start, delta.mu0, psi.mu0, delta.A0, psi.A0,
                                 delta.tune, psi.tune, robit = FALSE, df = 5) {
  
  n <- length(Y)
  k <- ncol(X)

  levels.treat <- as.numeric(names(table(treat)))
  
  ## "treat" variable should be a factor variable here where the base level is control
  tmax <- length(levels.treat) - 1
  keep <- thin + 1

  ## starting values, prior, and tuning parameters for sensitive item
  ## parameters: either the same starting values and prior for all
  ## sensitive items or a list with the names identical to the levels
  ## of the treatment factor variable
  if (is.list(delta.start)) {
    delta.start.all <- NULL
    for (i in 1:tmax) {
      delta.start.all <- c(delta.start.all, delta.start[[levels.treat[i+1]]])
    }
  } else {
    delta.start.all <- rep(delta.start, tmax)
  }

  if (is.list(delta.mu0)) {
    delta.mu0.all <- NULL
    for (i in 1:tmax) {
      delta.mu0.all <- c(delta.mu0.all, delta.mu0[[levels.treat[i+1]]])
    }
  } else {
    delta.mu0.all <- rep(delta.mu0, tmax)
  }
  
  if (is.list(delta.A0)) {
    delta.A0.all <- NULL
    for (i in 1:tmax) {
      delta.A0.all <- c(delta.A0.all, as.double(delta.A0[[levels.treat[i+1]]]))
    }
  } else {
    delta.A0.all <- rep(as.double(delta.A0), tmax)
  }

  if (is.list(delta.tune)) {
    delta.tune.all <- NULL
    for (i in 1:tmax) {
      delta.tune.all <- c(delta.tune.all, as.double(delta.tune[[levels.treat[i+1]]]))
    }
  } else {
    delta.tune.all <- rep(as.double(delta.tune), tmax)
  }

  ## number of parameters
  if (constrained) {
    n.par <- (tmax + 1) * (k + 1) 
  } else {
    n.par <- (tmax + 1) * (k + 1) + tmax
  }

  ## converting a treatment variable to an integer variale
  ##treat <- as.integer(treat) - 1
  ## passing the stuff to the C program
  if (robit) {
    ## no need for deltaCounter
    n.par <- n.par - tmax
    ## no unconstrained model for now
    if (!constrained)
      stop("only constrained model is allowed for robit regression")
    res <- .C("ictregBinomMultiRobit", as.integer(Y), as.integer(J),
              as.integer(n), as.integer(n.draws), as.integer(treat), as.integer(tmax), 
              as.double(X), as.double(delta.start.all),
              as.double(psi.start), as.integer(k),
              as.double(delta.mu0.all),
              as.double(psi.mu0), as.double(delta.A0.all),
              as.double(psi.A0), 
              as.double(psi.tune), as.integer(df),
              0, 0,
              as.integer(burnin), as.integer(keep), as.integer(verbose),
              allresults = double(n.par*floor((n.draws - burnin)/keep)),
              PACKAGE = "list")$allresults
  } else {
    res <- .C("ictregBinomMulti", as.integer(Y), as.integer(J),
              as.integer(n), as.integer(n.draws), as.integer(treat), as.integer(tmax), 
              as.double(X), as.double(delta.start.all),
              as.double(psi.start), as.integer(k), as.double(delta.mu0.all),
              as.double(psi.mu0), as.double(delta.A0.all),
              as.double(psi.A0), as.double(delta.tune.all),
              as.double(psi.tune), as.integer(!constrained),
              ##as.integer(ceiling), as.integer(floor),
              0,0,
              as.integer(burnin), as.integer(keep), as.integer(verbose),
              allresults = double(n.par*floor((n.draws - burnin)/keep)),
              PACKAGE = "list")$allresults
  }
    
  res <- matrix(res, byrow = TRUE, ncol = n.par) 

  class(res) <- "ictregBayesMulti"
  return(res)
            
}


ictregBayesMixed.fit <- function(Y, treat, X, Z, J, grp, constrained,
		     		 ##ceiling, floor,
                                 n.draws, burnin,
		     		 thin, verbose, delta.start,
		     		 psi.start, Sigma.start, Phi.start,
		     		 delta.mu0, psi.mu0, delta.A0, psi.A0,
		     		 Sigma.df, Sigma.scale, Phi.df,
		     		 Phi.scale, delta.tune, psi.tune,
		     		 gamma.tune, zeta.tune) {

  n <- length(Y)
  k <- ncol(X)
  m <- ncol(Z)
  n.grp <- length(table(grp))
  keep <- thin + 1
  alldraws <- floor((n.draws - burnin) / keep)
  ## fixed effects, Sigma, Phi, random effects, acceptance ratios
  n.par <- 2 * (k + m*(m + 1)/2 + n.grp * m + 1 + n.grp)    
  if (constrained == 2) 
    n.par <- n.par + 1
  
  ## this code assumes the equal number of obs within each group
  res <- .C("ictregBinomMixed", as.integer(Y), as.integer(J),
            as.integer(n), as.integer(n.draws), as.integer(treat), as.double(X),
            as.double(delta.start), as.double(psi.start), as.integer(k),
            as.double(delta.mu0), as.double(psi.mu0), as.double(delta.A0),
            as.double(psi.A0), as.double(delta.tune), as.double(psi.tune),
            as.integer(constrained),
            ##as.integer(ceiling), as.integer(floor),
            0,0,
            as.integer(grp-1), as.integer(n.grp), as.integer(max(table(grp))),
            as.double(t(Z)), as.integer(m), as.double(gamma.tune),
            as.double(zeta.tune), as.double(Sigma.start), as.double(Phi.start),
            as.integer(Sigma.df), as.double(Sigma.scale), as.integer(Phi.df),
            as.double(Phi.scale), as.integer(burnin), as.integer(keep),
            as.integer(verbose), allresults = double(n.par*alldraws), PACKAGE =
            "list")$allresults

  res <- matrix(res, byrow = TRUE, ncol = n.par)

  if (constrained == 2) {
    return(list(delta = res[, 1:k, drop = FALSE], psi = res[, (k + 1):(2*k + 1) , drop = FALSE], 
                Sigma = res[, (2*k + 2):(2*k + m*(m+1)/2 + 1) , drop = FALSE],
                Phi = res[, (2*k + m*(m+1)/2 + 2):(2*k + m*(m+1) + 1) , drop = FALSE], 
                gamma = res[, (2*k + m*(m+1) + 2):(2*k + m*(m+1) + n.grp*m + 1) , drop = FALSE], 
                zeta = res[, (2*k + m*(m+1) + n.grp*m + 2):(2*k + m*(m+1) + 2*n.grp*m + 1) , drop = FALSE], 
                delta.accept = res[alldraws, n.par - 2*n.grp - 1], 
                psi.accept = res[alldraws, n.par - 2*n.grp],
                gamma.accept = res[alldraws, (n.par - 2*n.grp + 1):(n.par - n.grp)],
                zeta.accept = res[alldraws, (n.par - n.grp + 1):n.par]))
  } else {
    return(list(delta = res[, 1:k, drop = FALSE], psi = res[, (k + 1):(2 * k) , drop = FALSE], 
                Sigma = res[, (2*k + 1):(2*k + m*(m+1)/2) , drop = FALSE],
                Phi = res[, (2*k + m*(m+1)/2 + 1):(2*k + m*(m+1)) , drop = FALSE], 
                gamma = res[, (2*k + m*(m+1) + 1):(2*k + m*(m+1) + n.grp*m) , drop = FALSE], 
                zeta = res[, (2*k + m*(m+1) + n.grp*m + 1):(2*k + m*(m+1) + 2*n.grp*m) , drop = FALSE], 
                delta.accept = res[alldraws, n.par - 2*n.grp - 1], 
                psi.accept = res[alldraws, n.par - 2*n.grp],
                gamma.accept = res[alldraws, (n.par - 2*n.grp + 1):(n.par - n.grp)],
                zeta.accept = res[alldraws, (n.par - n.grp + 1):n.par]))
  }
}



ictregBayesMultiMixed.fit <- function(Y, treat, X, Z, J, grp, constrained,
                                      ##ceiling, floor,
                                      n.draws, burnin,
                                      thin, verbose, delta.start, psi.start, Sigma.start, 
                                      Phi.start, delta.mu0, psi.mu0, delta.A0, psi.A0, 
                                      Sigma.df, Sigma.scale, Phi.df, Phi.scale, delta.tune, 
                                      psi.tune, gamma.tune, zeta.tune) {

  n <- length(Y)
  k <- ncol(X)
  m <- ncol(Z)
  n.grp <- length(table(grp))

  ## "treat" variable should be a factor variable here where the base level is control
  levels.treat <- as.numeric(names(table(treat)))
  tmax <- length(levels.treat) - 1  

  ## starting values, prior, and tuning parameters for sensitive item
  ## parameters: either the same starting values and prior for all
  ## sensitive items or a list with the names identical to the levels
  ## of the treatment factor variable
  if (is.list(delta.start)) {
    delta.start.all <- NULL
    for (i in 1:tmax) {
      delta.start.all <- c(delta.start.all, delta.start[[levels.treat[i+1]]])
    }
  } else {
    delta.start.all <- rep(delta.start, tmax)
  }

  if (is.list(delta.mu0)) {
    delta.mu0.all <- NULL
    for (i in 1:tmax) {
      delta.mu0.all <- c(delta.mu0.all, delta.mu0[[levels.treat[i+1]]])
    }
  } else {
    delta.mu0.all <- rep(delta.mu0, tmax)
  }
  
  if (is.list(delta.A0)) {
    delta.A0.all <- NULL
    for (i in 1:tmax) {
      delta.A0.all <- c(delta.A0.all, as.double(delta.A0[[levels.treat[i+1]]]))
    }
  } else {
    delta.A0.all <- rep(as.double(delta.A0), tmax)
  }

  if (is.list(delta.tune)) {
    delta.tune.all <- NULL
    for (i in 1:tmax) {
      delta.tune.all <- c(delta.tune.all, as.double(delta.tune[[levels.treat[i+1]]]))
    }
  } else {
    delta.tune.all <- rep(as.double(delta.tune), tmax)
  }

  if (is.list(gamma.tune)) {
    gamma.tune.all <- NULL
    for (i in 1:tmax) {
      gamma.tune.all <- c(gamma.tune.all, as.double(gamma.tune[[levels.treat[i+1]]]))
    }
  } else {
    gamma.tune.all <- rep(as.double(gamma.tune), tmax)
  }
  
  if (is.list(Sigma.start)) {
    Sigma.start.all <- NULL
    for (i in 1:tmax) {
      Sigma.start.all <- c(Sigma.start.all, as.double(Sigma.start[[levels.treat[i+1]]]))
    }
  } else {
    Sigma.start.all <- rep(as.double(Sigma.start), tmax)
  }

  if (is.list(Sigma.scale)) {
    Sigma.scale.all <- NULL
    for (i in 1:tmax) {
      Sigma.scale.all <- c(Sigma.scale.all, as.double(Sigma.scale[[levels.treat[i+1]]]))
    }
  } else {
    Sigma.scale.all <- rep(as.double(Sigma.scale), tmax)
  }

  if (is.list(Sigma.df)) {
    Sigma.df.all <- NULL
    for (i in 1:tmax) {
      Sigma.df.all <- c(Sigma.df.all, as.double(Sigma.df[[levels.treat[i+1]]]))
    }
  } else {
    Sigma.df.all <- rep(as.double(Sigma.df), tmax)
  }  

  ## fixed effects, Sigma, Phi, random effects, acceptance ratios
  keep <- thin + 1
  alldraws <- floor((n.draws - burnin) / keep)
  n.par <- (tmax + 1) * (k + m*(m + 1)/2 + n.grp * m + 1 + n.grp)    
  if (!constrained)
    n.par <- n.par + tmax 

  res <- .C("ictregBinomMultiMixed", as.integer(Y), as.integer(J),
            as.integer(n), as.integer(n.draws), as.integer(treat),
            as.integer(tmax), as.double(X),
            as.double(delta.start.all), as.double(psi.start),
            as.integer(k), as.double(delta.mu0.all),
            as.double(psi.mu0), as.double(delta.A0.all),
            as.double(psi.A0), as.double(delta.tune.all),
            as.double(psi.tune), as.integer(!constrained),
            ##as.integer(ceiling), as.integer(floor),
            0,0,
            as.integer(grp-1),
            as.integer(n.grp), as.integer(max(table(grp))),
            as.double(t(Z)), as.integer(m), as.double(gamma.tune.all),
            as.double(zeta.tune), as.double(Sigma.start.all),
            as.double(Phi.start), as.integer(Sigma.df.all),
            as.double(Sigma.scale.all), as.integer(Phi.df),
            as.double(Phi.scale), as.integer(burnin),
            as.integer(keep), as.integer(verbose), allresults =
            double(n.par*alldraws), PACKAGE = "list")$allresults

  res <- matrix(res, byrow = TRUE, ncol = n.par)

  if (constrained) {
    return(list(delta = res[, 1:(k*tmax), drop = FALSE], psi = res[, (k*tmax + 1):(k*(tmax+1)), drop = FALSE], 
                Sigma = res[, (k*(tmax+1) + 1):(k*(tmax+1) + m*(m+1)/2*tmax), drop = FALSE],
                Phi = res[, (k*(tmax+1) + m*(m+1)/2*tmax + 1):((k + m*(m+1)/2)*(tmax+1)) , drop = FALSE], 
                gamma = res[, ((k + m*(m+1)/2)*(tmax+1) + 1):((k + m*(m+1)/2)*(tmax+1) + n.grp*m*tmax) , drop = FALSE], 
                zeta = res[, ((k + m*(m+1)/2)*(tmax+1) + n.grp*m*tmax + 1):((k + m*(m+1)/2 + n.grp*m)*(tmax+1)) , drop = FALSE], 
                delta.accept = res[alldraws, (n.par - (tmax+1)*n.grp - tmax):(n.par - (tmax+1)*n.grp - 1)], 
                psi.accept = res[alldraws, n.par - (tmax+1)*n.grp],
                gamma.accept = res[alldraws, (n.par - (tmax+1)*n.grp + 1):(n.par - n.grp)],
                zeta.accept = res[alldraws, (n.par - n.grp + 1):n.par]))
  } else {
    return(list(delta = res[, 1:((k+1)*tmax), drop = FALSE], psi = res[, ((k+1)*tmax + 1):((k+1)*tmax + k) , drop = FALSE], 
                Sigma = res[, ((k+1)*tmax + k + 1):((k+1)*tmax + k + m*(m+1)/2*tmax) , drop = FALSE],
                Phi = res[, ((k+1)*tmax + k + m*(m+1)/2*tmax + 1):((k+1)*tmax + k + m*(m+1)/2*(tmax+1)) , drop = FALSE], 
                gamma = res[, ((k+1)*tmax + k + m*(m+1)/2*(tmax+1) + 1):((k+1)*tmax + k + m*(m+1)/2*(tmax+1) + n.grp*m*tmax) , drop = FALSE], 
                zeta = res[, ((k+1)*tmax + k + m*(m+1)/2*(tmax+1) + n.grp*m*tmax + 1):((k+1)*tmax + k + m*(m+1)/2*(tmax+1) + n.grp*m*(tmax+1)) , drop = FALSE],
                delta.accept = res[alldraws, (n.par - (tmax+1)*n.grp - tmax):(n.par - (tmax+1)*n.grp - 1)], 
                psi.accept = res[alldraws, n.par - (tmax+1)*n.grp],
                gamma.accept = res[alldraws, (n.par - (tmax+1)*n.grp + 1):(n.par - n.grp)],
                zeta.accept = res[alldraws, (n.par - n.grp + 1):n.par]))
  }
}


###
### Debugging function for Robit regression
###
Robit <- function(Y, X, beta.start, sims, beta0, A0, df) {
  
  tmp <- .C("R2Robit", as.integer(Y), as.double(X),
            as.double(beta.start), as.integer(nrow(X)), as.integer(ncol(X)),
            as.double(beta0), as.double(A0),
            as.integer(df), as.integer(sims), 
            betaStore = double(sims*ncol(X)),
            PACKAGE = "list")

  return(matrix(tmp$betaStore, byrow = TRUE, ncol = ncol(X)))

}


###
### Combining endorsement and list experiments
###

ictregBayesEndorse <- function( ## START LIST OPTIONS
                    Y.list,
                    treat.list,
                    ##X.list,  ## use data from endorse
                    J.list,
                    ##constrained,  ## not available this version
                    ##ceiling, floor,
                    ##n.draws,
                    ##burnin,
                    ##thin,
                    ##verbose,
                    ##delta.start,
                    psi.start,
                    ##delta.mu0,
                    psi.mu0,
                    ##delta.A0,
                    psi.A0,
                    ##delta.tune,
                    psi.tune,
                    robit = FALSE,
                    df = 5, ## END LIST OPTIONS
                    
                    Y,
                    data,
                    treat = NA,
                    na.strings = 99,
                    covariates = FALSE,
                    identical.lambda = FALSE,
                    formula = NA,
                    x.start = 0,
                    s.start = 0,
                    beta.start = 1,
                    tau.start = NA,
                    lambda.start = 0,
                    omega2.start = .1,
                    ##theta.start = 0,
                    ##phi2.start = .1,
                    delta.start = 0,
                    mu.beta = 0,
                    mu.x = 0,
                    mu.lambda = 0,
                    mu.delta = 0,
                    precision.beta = 0.04,
                    precision.x = 1,
                    precision.lambda = 0.04,
                    precision.delta = 0.04,
                    s0.omega2= 3,
                    nu0.omega2 = 5,
                    s0.phi2 = 3,
                    nu0.phi2 = 5,
                    MCMC = 20000,
                    burn = 1000,
                    thin = 1,
                    mda = FALSE,
                    mh = TRUE,
                    prop = 0.001,
                    x.sd = TRUE,
                    tau.out = FALSE,
                    s.out = FALSE,
                    verbose = TRUE,
                    seed.store = FALSE,
                    update = FALSE,
                    update.start = NULL,
                    probit = FALSE) {

  ## stuff added for list experiment

  ##n <- length(Y)
  ##k <- ncol(X)

  levels.treat <- as.numeric(names(table(treat.list)))
  
  tmax <- length(levels.treat) - 1
  ##keep <- thin + 1

  cov.mat <- model.matrix(formula, data)

  M <- ncol(cov.mat)
  
  n.par <- (tmax + 1) * (M + 1) - tmax
  
  ##
  
  if (covariates == FALSE) formula <- ~ 1


  ##cov.mat <- model.matrix(formula, data)
  ##M <- ncol(cov.mat)

  data <- data[complete.cases(cov.mat),]
  cov.mat <- cov.mat[complete.cases(cov.mat), ]

  N <- nrow(data)
  J <- length(Y)
  
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
    tau.start <- update.start$tau.start
    x.start <- update.start$x.start
    s.start <- update.start$s.start
    lambda.start <- update.start$lambda.start
    omega2.start <- update.start$omega2.start
    psi.start <- update.start$psi.start
    delta.start <- update.start$delta.start
    .Random.seed <- update.start$seed

  } else {
    if (is.na(tau.start[1])) {
      tau.start <- matrix(-99, nrow = J, ncol = max.L - 1)

      for (j in 1:J){
        temp.tau <- seq(from = 0, to = .5 * (L[j] - 2), by = .5)
        for (i in 1:(L[j] - 1)) {
          tau.start[j, i] <- temp.tau[i]
        }
        ##tau.start[j, L[j]] <- max(temp.tau) + 1000
      }
    }
    
    if (length(x.start) == 1) {
      x.start <- rep(x.start, times = N)
    }

    if (length(s.start) == 1) {
      s.start <- matrix(s.start, nrow = N, ncol = J)
    }

    if (length(beta.start) == 1) {
      beta.start <- matrix(beta.start, nrow = J, ncol = 2)
    }

    if (length(lambda.start) == 1) {
      lambda.start <- matrix(lambda.start, nrow = M, ncol = K)
    }

    if (length(omega2.start) == 1) {
      omega2.start <- matrix(omega2.start, nrow = 1, ncol = K)
    }

    if (length(delta.start) == 1) {
      delta.start <- rep(delta.start, times = M)
    }

  }
  
  if (length(mu.beta) == 1) {
    mu.beta <- matrix(mu.beta, nrow = J, ncol = 2)
  }

  if (length(mu.lambda) == 1) {
    mu.lambda <- rep(mu.lambda, times = M) 
  }

  if (covariates) {
    if (length(mu.delta) == 1) {
      mu.delta <- rep(mu.delta, times = M)
    }
  } else {
    if (length(mu.x) == 1) {
      mu.delta <- rep(mu.x, times = N)
    } else {
      mu.delta <- mu.x
    }
  }

  if (length(precision.beta) == 1) {
    precision.beta <- diag(precision.beta, nrow = 2, ncol = 2)
  }

  if (length(precision.delta) == 1) {
    precision.delta <- diag(precision.delta, nrow = M, ncol = M)
  }
  
  if (length(precision.lambda) == 1) {
    precision.lambda <- diag(precision.lambda, nrow = M, ncol = M) 
  }

  if (length(prop) != J) {
    prop <- rep(prop, times = J)
  }

  printout <- floor( (MCMC - burn) / (thin + 1) )

  n.par <- M + 1

  if(probit == TRUE) {

    ## removed and b here
    
    temp <- .C("CombineEndorseListProbit",
               as.integer(response),
               as.integer(endorse),
               as.double(cov.mat),
               as.integer(N),
               as.integer(J),
               as.integer(M),
               as.integer(K),
               as.integer(L),
               as.integer(max.L),
               as.double(x.start),
               as.double(s.start),
               as.double(beta.start),
               as.double(tau.start),
               as.double(lambda.start),
               as.double(omega2.start),
               as.double(delta.start),
               as.double(mu.beta),
               as.double(precision.beta),
               as.double(precision.x),
               as.double(mu.lambda),
               as.double(precision.lambda),
               as.double(mu.delta),
               as.double(precision.delta),
               as.double(s0.omega2),
               as.integer(nu0.omega2),
               as.integer(MCMC),
               as.integer(burn),
               as.integer(thin),
               as.integer(mda),
               as.integer(mh),
               as.double(prop),
               as.integer(x.sd),
               as.integer(tau.out),
               as.integer(s.out),
               as.integer(covariates),
               betaStore = double(printout * 2 * J),
               tauStore = if (tau.out) double(printout * (max.L - 1) * J) else double(1),
               xStore = if (x.sd & seed.store) double(printout + N) else if (x.sd & !(seed.store)) double(printout) else double(printout * N),
               sStore = if (s.out & seed.store) double(printout * N * J) else if (!(s.out) & seed.store) double(N * J) else double(1),
               lambdaStore = double(printout * K * M),
               deltaStore = if (covariates) double(printout * M) else double(1),
               omega2Store = if (seed.store) double(K * printout) else double(1),
               accept.ratio = double(J),
               ## list parts
               as.integer(Y.list),
               as.integer(J.list),
               as.integer(treat.list),
               as.double(psi.start), 
               as.double(psi.mu0), 
               as.double(psi.A0), 
               as.double(psi.tune), 
               0,0, # ceiling and floor effects
               as.integer(verbose),
               psiStore = double((M+1) * printout),
               package = "list")
  }

  seedStore <- .Random.seed

  res <- list(beta = matrix(as.double(temp$betaStore), byrow = TRUE, ncol = 2 * J, nrow = printout),
              tau = if (tau.out) matrix(as.double(temp$tauStore), byrow = TRUE, ncol = (max.L - 1) * J, nrow = printout) else NULL,
              x = if (x.sd) matrix(as.double(temp$xStore), ncol = (N + 1), nrow = printout)
              else matrix(as.double(temp$xStore), byrow = TRUE, ncol = N, nrow = printout),
              s = if (s.out) matrix(as.double(temp$sStore), byrow = TRUE, ncol = N * J, nrow = printout) else NULL,
              lambda = matrix(as.double(temp$lambdaStore), byrow = TRUE, ncol = M * K, nrow = printout),
              psi = if (covariates) matrix(as.double(temp$psiStore), byrow = TRUE, ncol = M + 1, nrow = printout)[, 1:M] else NULL,
              delta = if (covariates) matrix(as.double(temp$deltaStore), byrow = TRUE, ncol = M, nrow = printout) else NULL,
              omega2 = matrix(as.double(temp$omega2Store), byrow = TRUE, ncol = K, nrow = printout),
              accept.ratio = if (mh) as.double(temp$accept.ratio) else NULL,
              seed = if (seed.store) seedStore else NULL,
              beta.start = if (seed.store) matrix(as.double(temp$betaStore)[(2 * J * (printout - 1) + 1):(2 * J * printout)],
                nrow = J, ncol = 2, byrow = TRUE) else NULL,
              tau.start = if (seed.store) matrix(as.double(temp$tauStore)[(length(as.double(temp$tauStore)) - (max.L - 1) * J + 1):(length(as.double(temp$tauStore)))],
                nrow = J, ncol = max.L - 1, byrow = TRUE) else NULL,
              x.start = if (seed.store) as.double(temp$xStore)[(length(as.double(temp$xStore)) - N + 1):(length(as.double(temp$xStore)))] else NULL,
              s.start = if (seed.store) matrix(as.double(temp$sStore)[(length(as.double(temp$sStore)) - (N * J) + 1):length(as.double(temp$sStore))],
                nrow = N, ncol = J, byrow = TRUE) else NULL,
              lambda.start = if (seed.store) matrix(as.double(temp$lambdaStore)[(K * M * (printout - 1) + 1):(K * M * printout)],
                nrow = M, ncol = K, byrow = FALSE)  else NULL,
              delta.start = if (seed.store & covariates) as.double(temp$deltaStore)[((printout - 1) * M + 1):(printout * M)] else NULL,
              psi.start = if (seed.store) as.double(temp$psiStore)[((printout - 1) * (M + 1) + 1):(printout * (M + 1) - 1)] else NULL,
              omega2.start = if (seed.store) matrix(as.double(temp$omega2Store)[(K * (printout - 1) + 1):(K * printout)],
                byrow = TRUE, ncol = K, nrow = 1) else NULL 
              )
  
  colnames(res$beta) <- paste(rep(c("alpha", "beta"), times = J),
                              rep(1:J, each = 2), sep = ".")
  res$beta <- mcmc(res$beta, start = 1, end = printout, thin = 1)

  if (tau.out) {
    temp.names <- paste("tau", 1:J, sep = "")
    colnames(res$tau) <- paste(rep(temp.names, each = (max.L - 1)),
                               rep(1:(max.L - 1), times = J), sep = ".")
    res$tau <- mcmc(res$tau, start =  1, end = printout, thin = 1)
  }

  if (x.sd) {
    colnames(res$x) <- "sd.x"
  } else {
    colnames(res$x) <- paste("x", 1:N, sep = ".")
    res$x <- mcmc(res$x, start = 1, end = printout, thin = 1)
  }

  if (s.out) {
    colnames(res$s) <- paste("s", rep(1:nrow(data), each = J),
                             rep(1:J, times = nrow(data)), sep = "")
  }

  
  if (identical.lambda) {
    temp.names <- paste("lambda", 1:K, sep = "")
    if (covariates) {
      colnames(res$lambda) <- paste(rep(temp.names, each = M),
                                    rep(1:M, times = K),
                                    sep = ".")
    } else {
      colnames(res$lambda) <- temp.names
    }
  } else {
    temp.names <- paste("lambda", rep(1:J, each = K), rep(1:K, times = J),
                        sep = "")
    if (covariates) {
      colnames(res$lambda) <- paste(rep(temp.names, each = M),
                                    rep(1:M, times = (J * K)),
                                    sep = ".")
    } else {
      colnames(res$lambda) <- temp.names
    }
  }
  res$lambda <- mcmc(res$lambda, start = 1, end = printout, thin = 1)


  if (identical.lambda) {
    colnames(res$omega2) <- paste("omega2.", 1:K, sep = "")
  } else {
    colnames(res$omega2) <- paste("omega2.", rep(1:J, each = K), rep(1:K,
                                                        times = J), sep = "")
  }
  res$omega2 <- mcmc(res$omega2, start = 1, end = printout, thin = 1)  

  if (identical.lambda == FALSE) { 
    temp.names <- paste("theta", 1:K, sep = "")
    if (covariates) {
      colnames(res$theta) <- paste(rep(temp.names, each = M), rep(1:M, times = K),
                                   sep = ".")
    } else {
      colnames(res$theta) <- temp.names
    }
    res$theta <- mcmc(res$theta, start = 1, end = printout, thin = 1)

    temp.names <- paste("phi2.", 1:K, sep = "")
    if (covariates) {
      colnames(res$phi2) <- paste(rep(temp.names, each = M), rep(1:M, times = K),
                                  sep = ".")
    } else {
      colnames(res$phi2) <- temp.names
    }
    res$phi2 <- mcmc(res$phi2, start = 1, end = printout, thin = 1)
  }

  if (covariates) {
    colnames(res$delta) <- paste("delta", 1:M, sep = ".")
    res$delta <- mcmc(res$delta, start = 1, end = printout, thin = 1)
  }
  
  if (covariates) {
    colnames(res$psi) <- paste("psi", 1:M, sep = ".")
    res$psi <- mcmc(res$psi, start = 1, end = printout, thin = 1)
  }

  names(res$accept.ratio) <- paste("Q", 1:J, sep = "")

  return(res)
}
