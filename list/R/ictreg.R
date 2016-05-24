
ictreg <- function(formula, data = parent.frame(), treat = "treat", J, method = "ml", weights, 
                   overdispersed = FALSE, constrained = TRUE, floor = FALSE, ceiling = FALSE, 
                   ceiling.fit = "glm", floor.fit = "glm", ceiling.formula = ~ 1, floor.formula = ~ 1, 
                   fit.start = "lm", fit.nonsensitive = "nls", multi.condition = "none", maxIter = 5000, verbose = FALSE, ...){

  ictreg.call <- match.call()

  # set up data frame, with support for standard and modified responses
  mf <- match.call(expand.dots = FALSE)
  
  # make all other call elements null in mf <- NULL in next line
  mf$method <- mf$maxIter <- mf$verbose <- mf$fit.start <- mf$J <- mf$design <- 
    mf$treat <- mf$weights <- mf$constrained <- mf$overdispersed <- mf$floor <- 
    mf$ceiling <- mf$ceiling.fit <- mf$fit.nonsensitive <- mf$floor.fit <- 
    mf$multi.condition <- mf$floor.formula <- mf$ceiling.formula <- NULL
  mf[[1]] <- as.name("model.frame")
  mf$na.action <- 'na.pass'
  mf <- eval.parent(mf)
  
  if(floor == TRUE | ceiling == TRUE) {
    boundary <- TRUE
  } else {
    boundary <- FALSE
  }
  
  # define design, response data frames
  x.all <- model.matrix.default(attr(mf, "terms"), mf)
  y.all <- model.response(mf)
  
  if(class(y.all)=="matrix") {
    design <- "modified"
  } else {
    design <- "standard"
  }
  
  if(missing("weights") == FALSE) {
      weighted <- TRUE
      w.all <- data[, paste(weights)]
  } else {
      weighted <- FALSE
      w.all <- rep(1, length(y.all))
  }

  if (method == "nls") fit.start <- "nls"
  
  # extract number of non-sensitive items from the response matrix
  if(design=="modified") J <- ncol(y.all) - 1
  
  # set -777 code for treatment status
#  if(design=="modified"){
#    for(j in 1:J) y.all[,j] <- ifelse(!is.na(y.all[,J+1]), -777, y.all[,j])
#    y.all[,J+1] <- ifelse(is.na(y.all[,J+1]), -777, y.all[,J+1])
#  }
  
  # list-wise missing deletion
  na.x <- apply(is.na(x.all), 1, sum)
  if(design=="standard") na.y <- is.na(y.all)
  if(design=="modified") {
    na.y <- apply(is.na(y.all[,1:J]), 1, sum)
    na.y[is.na(y.all[,J+1]) == FALSE] <- 0
  }
  na.w <- is.na(w.all)

  ## treatment indicator for subsetting the dataframe
  if(design=="standard") t <- data[na.x==0 & na.y==0 & na.w==0, paste(treat)]
  if(design=="modified") t <- as.numeric(is.na(y.all[na.x == 0 & na.y == 0 & na.w==0, J+1]) == FALSE)

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
    ##condition.labels <- as.character(paste("Sensitive Item",sort(unique(t))))
    ##condition.labels[which(condition.labels == "Sensitive Item 0")] <- "control"
    ##treatment.labels <- condition.labels[condition.labels != "control"]

    condition.labels <- sort(unique(t))
    treatment.labels <- condition.labels[condition.labels != 0]
    control.label <- 0
  }
  
  # list wise delete
  if(design=="standard") y.all <- y.all[na.x==0 & na.y==0 & na.w==0]
  if(design=="modified") y.all <- y.all[na.x==0 & na.y==0 & na.w==0,]
  x.all <- x.all[na.x==0 & na.y==0, , drop = FALSE]
  w.all <- w.all[na.x==0 & na.y==0 & na.w==0]

  ## so that the output data has the same dimension as x.all and y.all
  data <- data[na.x==0 & na.y==0 & na.w==0, , drop = FALSE]
  
  # set up data objects for y and x for each group from user input
  x.treatment <- x.all[t != 0, , drop = FALSE]
  y.treatment <- subset(y.all, t != 0)
  w.treatment <- subset(w.all, t != 0)
  x.control <- x.all[t == 0 , , drop = FALSE]
  y.control <- subset(y.all, t==0)
  w.control <- subset(w.all, t==0)

  if(design=="standard" & missing("J")) {
    J <- max(y.treatment) - 1
    cat(paste("Note: number of control items set to", J, "\n"))
  }

  condition.values <- sort(unique(t))
  treatment.values <- 1:length(condition.values[condition.values!=0])
  
  if(length(treatment.values) > 1) {
    multi <- TRUE
  } else {
    multi <- FALSE
  }

  if(weighted == TRUE) {
    if(design == "modified")
      stop("Weighted list experiment regression is not supported (yet) for the modified design.")
    if(boundary == TRUE)
      stop("Weighted list experiment regression is not supported (yet) for the ceiling/floor design.")
    if(multi == TRUE)
      stop("Weighted list experiment regression is not supported (yet) for the multi-item design.")
  }

  n <- nrow(x.treatment) + nrow(x.control)

  coef.names <- colnames(x.all)

  nPar <- ncol(x.all)

  intercept.only <- ncol(x.all)==1 & sum(x.all[,1]==1) == n

  if (intercept.only == TRUE) {
    x.vars <- "1"
  } else {
    x.vars <- coef.names[-1]
  }
  
  logistic <- function(x) exp(x)/(1+exp(x))
  logit <- function(x) return(log(x)-log(1-x))

  if (design == "standard" & method == "lm") {
    
    treat <- matrix(NA, nrow = n, ncol = length(treatment.values))
    
    for (m in 1:length(treatment.values)) 
      treat[ , m] <- as.numeric(t == m)
      
    x.all.noint <- x.all[, -1, drop = FALSE]

    x.all.lm <- x.all
    for (m in 1:length(treatment.values))
      x.all.lm <- cbind(x.all.lm, x.all * treat[, m])

    if (intercept.only == TRUE) {
      fit.lm <- lm(y.all ~ treat, weights = w.all)
    } else {
      ## fit.lm <- lm(y.all ~ x.all.noint * treat)
      fit.lm <- lm(y.all ~ x.all.lm - 1, weights = w.all)
    }
    
    vcov <- vcovHC(fit.lm)

    par.control <- coef(fit.lm)[1:length(coef.names)]
    se.control <- sqrt(diag(vcov))[1:length(coef.names)]

    names(par.control) <- names(se.control) <- coef.names
    
    par.treat <- se.treat <- list()
    
    for (m in 1:length(treatment.values)) {
      if (intercept.only == TRUE) {
        par.treat[[m]] <- coef(fit.lm)[nPar + m]
        se.treat[[m]] <- sqrt(diag(vcov))[nPar + m]
        names(par.treat[[m]]) <- names(se.treat[[m]]) <- "(Intercept)"
      } else {
        
        par.treat[[m]] <- coef(fit.lm)[(nPar + (m-1) * nPar + 1) : (nPar + m * nPar)]
        se.treat[[m]] <- sqrt(diag(vcov))[(nPar + (m-1) * nPar + 1) : (nPar + m * nPar)]
        
        names(par.treat[[m]]) <- names(se.treat[[m]]) <- coef.names
      }
      
    }

    if(multi == FALSE) {
      par.treat <- par.treat[[1]]
      se.treat <- se.treat[[1]]
    }
    
    sum.fit.lm <- summary(fit.lm)
      
    resid.se <- sum.fit.lm$sigma
    resid.df <- sum.fit.lm$df[2]
    
    if (multi == TRUE) {
      return.object <- list(par.treat=par.treat, se.treat=se.treat,
                            par.control=par.control, se.control=se.control,
                            vcov=vcov, treat.values = treatment.values, treat.labels = treatment.labels,
                            J=J, coef.names=coef.names,
                            resid.se = resid.se, resid.df = resid.df,  design = design, method = method,
                            multi = multi, boundary = boundary, call = match.call(),
                            data = data, x = x.all, y = y.all, treat = t)
      if(weighted == TRUE)
          return.object$weights <- w.all
    } else {
      return.object <- list(par.treat=par.treat, se.treat=se.treat,
                            par.control=par.control, se.control=se.control,
                            vcov=vcov, J=J,  coef.names=coef.names, resid.se = resid.se,
                            resid.df = resid.df,  design = design, method = method,
                            multi = multi, boundary = boundary, call = match.call(),
                            data = data, x = x.all, y = y.all, treat = t)
      if(weighted == TRUE)
          return.object$weights <- w.all
    }
    
  } else if (design=="standard" & method != "lm") {
    
    ##
    # Run two-step estimator
    
    vcov.twostep.std <- function(treat.fit, control.fit, J, y, treat, x) {
	
      n <- length(y)
      y1 <- y[treat == 1]
      y0 <- y[treat == 0]
      x1 <- x[treat == 1, , drop = FALSE]
      x0 <- x[treat == 0, , drop = FALSE]
      delta <- coef(treat.fit)
      gamma <- coef(control.fit)
      
      m1 <- c((y1 - J*logistic(x1 %*% gamma) - logistic(x1 %*% delta)) *
              logistic(x1 %*% delta)/(1+exp(x1 %*% delta))) * x1
      m0 <- c((y0 - J*logistic(x0 %*% gamma)) * J *
              logistic(x0 %*% gamma)/(1+exp(x0 %*% gamma))) * x0
      Em1 <- t(m1) %*% m1 / n
      Em0 <- t(m0) %*% m0 / n
      F <- adiag(Em1, Em0)
      Gtmp <- c(logistic(x1 %*% delta)/(1 + exp(x1 %*% delta))) * x1
      G1 <- -t(Gtmp) %*% Gtmp / n  
      Gtmp <- c(sqrt(J*logistic(x1 %*% delta)*logistic(x1 %*% gamma)/
                     ((1 + exp(x1 %*% delta))*(1 + exp(x1 %*% gamma))))) * x1
      G2 <- -t(Gtmp) %*% Gtmp / n
      Gtmp <- c(J*logistic(x0 %*% gamma)/(1 + exp(x0 %*% gamma))) * x0
      G3 <- -t(Gtmp) %*% Gtmp / n
      invG1 <- solve(G1)
      invG3 <- solve(G3)
      invG <- rbind(cbind(invG1, - invG1 %*% G2 %*% invG3),
                    cbind(matrix(0, ncol = ncol(G1), nrow = nrow(G3)), invG3))
      
      return(invG %*% F %*% t(invG) / n)
      
    }
    
    ## fit to the control group
    fit.glm.control <- glm(cbind(y.control, J-y.control) ~ x.control - 1,
                           family = binomial(logit), weights = w.control)
    
    coef.glm.control <- coef(fit.glm.control)
    names(coef.glm.control) <- paste("beta", 1:length(coef.glm.control), sep = "")

    if(fit.start == "nls") {
      if(is.null(ictreg.call$control)){
        fit.control <- nls( as.formula(paste("I(y.control/J) ~ logistic(x.control %*% c(",
                                             paste(paste("beta", 1:length(coef.glm.control),
                                                         sep=""), collapse= ","), "))")) ,
                           start = coef.glm.control, weights = w.control, control =
                           nls.control(maxiter=maxIter, warnOnly=TRUE), ... = ...)
      } else {
        fit.control <- nls( as.formula(paste("I(y.control/J) ~ logistic(x.control %*% c(",
                                             paste(paste("beta", 1:length(coef.glm.control),
                                                         sep=""), collapse= ","), "))")) ,
                           weights = w.control,
                           start = coef.glm.control, ... = ...)
      }
      
      fit.control.coef <- summary(fit.control)$parameters[,1]
      
    } else if (fit.start == "glm") {
      fit.control <- fit.glm.control
      fit.control.coef <- coef(fit.glm.control)
    } else if (fit.start == "lm") {
      fit.control <- lm(y.control ~ x.control - 1, weights = w.control)
      fit.control.coef <- coef(fit.glm.control)
    }
    
    fit.treat <- vcov.nls <- se.twostep <- par.treat.nls.std <- list()

    for(m in 1:length(treatment.values)){

      curr.treat <- t[t!=0] == treatment.values[m]

      x.treatment.curr <- x.treatment[curr.treat, , drop = F]
      w.treatment.curr <- w.treatment[curr.treat]

      ## calculate the adjusted outcome    
      y.treatment.pred <- y.treatment[curr.treat] - logistic(x.treatment.curr
                                                             %*% fit.control.coef)*J
      
      ## fit to the treated		
      y.treatment.pred.temp <- ifelse(y.treatment.pred>1, 1, y.treatment.pred)
      y.treatment.pred.temp <- ifelse(y.treatment.pred.temp<0, 0, y.treatment.pred.temp)
      
      alpha <- mean(y.treatment.pred.temp)
      
      y.treatment.start <- ifelse(y.treatment.pred.temp >
                                  quantile(y.treatment.pred.temp, alpha), 1, 0)
      
      try(fit.glm.treat <- glm(cbind(y.treatment.start, 1 - y.treatment.start) ~ x.treatment.curr - 1,
                               family = binomial(logit), weights = w.treatment.curr), silent = F)
      try(coef.glm.treat <- coef(fit.glm.treat), silent = T)
      try(names(coef.glm.treat) <- paste("beta", 1:length(coef.glm.treat), sep = ""), silent = T)

      if(fit.start == "nls") {
        if(exists("coef.glm.treat")) {
          if (is.null(ictreg.call$control)) {
            fit.treat[[m]] <- nls( as.formula(paste("y.treatment.pred ~ logistic(x.treatment.curr %*% c(",
                                                    paste(paste("beta", 1:length(coef.glm.treat),
                                                                sep=""), collapse= ","), "))")) ,
                                  start = coef.glm.treat, weights = w.treatment.curr, control =
                                  nls.control(maxiter = maxIter, warnOnly = TRUE), ... = ...)
          } else {
            fit.treat[[m]] <- nls( as.formula(paste("y.treatment.pred ~ logistic(x.treatment.curr %*% c(",
                                                    paste(paste("beta", 1:length(coef.glm.treat),
                                                                sep=""), collapse= ","), "))")) ,
                                  start = coef.glm.treat, weights = w.treatment.curr, ... = ...)
          }
        } else {
          if (is.null(ictreg.call$control)) {
            fit.treat[[m]] <- nls( as.formula(paste("y.treatment.pred ~ logistic(x.treatment.curr %*% c(",
                                                    paste(paste("beta", 1:length(coef.glm.treat),
                                                                sep=""), collapse= ","), "))")),
                                  weights = w.treatment.curr,
                                  control = nls.control(maxiter = maxIter, warnOnly = TRUE),
                                  ... = ...)
          } else {
            fit.treat[[m]] <- nls( as.formula(paste("y.treatment.pred ~ logistic(x.treatment.curr %*% c(",
                                                    paste(paste("beta", 1:length(coef.glm.treat),
                                                                sep=""), collapse= ","), "))")) ,
                                  weights = w.treatment.curr,
                                  ... = ...)
          }
        }
        
      } else if (fit.start == "glm") {
        fit.treat[[m]] <- fit.glm.treat
      } else if (fit.start == "lm") {
        fit.treat[[m]] <- lm(y.treatment.pred ~ x.treatment.curr - 1, weights = w.treatment.curr)
      }

      sample.curr <- t==0 | t==treatment.values[m]
      
      vcov.nls[[m]] <- vcov.twostep.std(fit.treat[[m]], fit.control, J,
                               y.all[sample.curr], t[sample.curr]>0, x.all[sample.curr, , drop = FALSE])

      se.twostep[[m]] <- sqrt(diag(vcov.nls[[m]]))
      par.treat.nls.std[[m]] <- coef(fit.treat[[m]])
      
    }

    if(multi == FALSE) {
      fit.treat <- fit.treat[[1]]
      vcov.nls <- vcov.nls[[1]]
      se.twostep <- se.twostep[[1]]
      par.treat.nls.std <- par.treat.nls.std[[m]]
    }
    
    par.control.nls.std <- coef(fit.control)

    if(method=="ml") {
      
      if(boundary == FALSE) {
      
        if (multi == FALSE) {

          coef.control.start <- par.control.nls.std
          coef.treat.start <- par.treat.nls.std
          
          ## ##
          ## Run EM estimator if requested
          
          ##
          ## observed-data log-likelihood for beta binomial
          ##
          obs.llik.std <- function(par, J, y, treat, x, wt, const = FALSE) {
            
            k <- ncol(x)
            if (const) {
              coef.h0 <- coef.h1 <- par[1:k]
              rho.h0 <- rho.h1 <- logistic(par[(k+1)])
              coef.g <- par[(k+2):(2*k+1)]
            } else {
              coef.h0 <- par[1:k]
              rho.h0 <- logistic(par[(k+1)])
              coef.h1 <- par[(k+2):(2*k+1)]
              rho.h1 <- logistic(par[(2*k+2)])
              coef.g <- par[(2*k+3):(3*k+2)]
            }
            
            h0X <- logistic(x %*% coef.h0)
            h1X <- logistic(x %*% coef.h1)
            gX <- logistic(x %*% coef.g)
            
            ind10 <- ((treat == 1) & (y == 0))
            ind1J1 <- ((treat == 1) & (y == (J+1)))
            ind1y <- ((treat == 1) & (y > 0) & (y < (J+1)))
            
            if (sum(ind10) > 0) {
              p10 <- sum(wt[ind10] * (log(1-gX[ind10]) + dBB(x = 0, bd = J, mu = h0X[ind10],
                                                sigma = rho.h0/(1-rho.h0), log = TRUE)))
            } else {
              p10 <- 0
            }
            
            if (sum(ind1J1) > 0) {
              p1J1 <- sum(wt[ind1J1] * (log(gX[ind1J1]) + dBB(J, bd = J, mu = h1X[ind1J1],
                                                sigma = rho.h1/(1-rho.h1), log = TRUE)))
            } else {
              p1J1 <- 0
            }
            
            if (sum(ind1y) > 0) {
              p1y <- sum(wt[ind1y] * (log(gX[ind1y]*dBB(y[ind1y]-1, bd = J, mu = h1X[ind1y], sigma =
                                           rho.h1/(1-rho.h1), log = FALSE) + (1-gX[ind1y])*
                             dBB(y[ind1y], bd = J, mu = h0X[ind1y],
                                 sigma = rho.h0/(1-rho.h0), log = FALSE))))
            } else {
              p1y <- 0
            }
            
            if (sum(treat == 0) > 0) {
              p0y <- sum(wt[!treat] * (log(gX[!treat]*dBB(y[!treat], bd = J, mu = h1X[!treat],
                                            sigma = rho.h1/(1-rho.h1), log = FALSE) +
                             (1-gX[!treat])*dBB(y[!treat], bd = J, mu = h0X[!treat],
                                                sigma = rho.h0/(1-rho.h0), log = FALSE))))
            } else {
              p0y <- 0
            }
            
            return(p10+p1J1+p1y+p0y)
            
          }
          
          ##
          ##  Observed data log-likelihood for binomial
          ##
          obs.llik.binom.std <- function(par, J, y, treat, x, wt, const = FALSE) {
            
            k <- ncol(x)
            if (const) {
              coef.h0 <- coef.h1 <- par[1:k]
              coef.g <- par[(k+1):(2*k)]
            } else {
              coef.h0 <- par[1:k]
              coef.h1 <- par[(k+1):(2*k)]
              coef.g <- par[(2*k+1):(3*k)]
            }
            
            h0X <- logistic(x %*% coef.h0)
            h1X <- logistic(x %*% coef.h1)
            gX <- logistic(x %*% coef.g)
            
            ind10 <- ((treat == 1) & (y == 0))
            ind1J1 <- ((treat == 1) & (y == (J+1)))
            ind1y <- ((treat == 1) & (y > 0) & (y < (J+1)))
            
            if (sum(ind10) > 0) {
              p10 <- sum(wt[ind10] * (log(1-gX[ind10]) + dbinom(x = 0, size = J, prob = h0X[ind10], log = TRUE)))
            } else {
              p10 <- 0
            }
            if (sum(ind1J1) > 0) {
              p1J1 <- sum(wt[ind1J1] * (log(gX[ind1J1]) + dbinom(J, size = J, prob = h1X[ind1J1], log = TRUE)))
            } else {
              p1J1 <- 0
            }
            if (sum(ind1y) > 0) {
              p1y <- sum(wt[ind1y] * (log(gX[ind1y]*dbinom(y[ind1y]-1, size = J, prob = h1X[ind1y],
                                              log = FALSE) + (1-gX[ind1y])*
                             dbinom(y[ind1y], size = J, prob = h0X[ind1y], log = FALSE))))
            } else {
              p1y <- 0
            }
            
            if (sum(treat == 0) > 0) {
              p0y <- sum(wt[!treat] * (log(gX[!treat]*dbinom(y[!treat], size = J, prob = h1X[!treat], log = FALSE) +
                             (1-gX[!treat])*dbinom(y[!treat], size = J, prob = h0X[!treat], log = FALSE))))
            } else {
              p0y <- 0
            }
            
            return(p10+p1J1+p1y+p0y)
          }
          
          ##
          ## EM algorithm
          ##
          
          ## Estep for beta-binomial
          Estep.std <- function(par, J, y, treat, x, const = FALSE) {
            
            k <- ncol(x)
            if (const) {
              coef.h0 <- coef.h1 <- par[1:k]
              rho.h0 <- rho.h1 <- logistic(par[(k+1)])
              coef.g <- par[(k+2):(2*k+1)]
            } else {
              coef.h0 <- par[1:k]
              rho.h0 <- logistic(par[(k+1)])
              coef.h1 <- par[(k+2):(2*k+1)]
              rho.h1 <- logistic(par[(2*k+2)])
              coef.g <- par[(2*k+3):(3*k+2)]
            }
            
            h0X <- logistic(x %*% coef.h0)
            h1X <- logistic(x %*% coef.h1)
            gX <- logistic(x %*% coef.g)
            
            ind <- !((treat == 1) & ((y == 0) | (y == (J+1))))
            w <- rep(NA, length(y))
            w[ind] <- gX[ind]*dBB((y-treat)[ind], bd = J, mu = h1X[ind], sigma = rho.h1/(1-rho.h1), log = FALSE) /
              (gX[ind]*dBB((y-treat)[ind], bd = J, mu = h1X[ind], sigma = rho.h1/(1-rho.h1), log = FALSE) +
               (1-gX[ind])*dBB(y[ind], bd = J, mu = h0X[ind], sigma = rho.h0/(1-rho.h0), log = FALSE))
            
            w[(treat == 1) & (y == 0)] <- 0
            w[(treat == 1) & (y == (J+1))] <- 1
            
            return(w)
          }
          
          
          ## Estep for binomial
          Estep.binom.std <- function(par, J, y, treat, x, const = FALSE) {
            
            k <- ncol(x)
            if (const) {
              coef.h0 <- coef.h1 <- par[1:k]
              coef.g <- par[(k+1):(2*k)]
            } else {
              coef.h0 <- par[1:k]
              coef.h1 <- par[(k+1):(2*k)]
              coef.g <- par[(2*k+1):(3*k)]
            }
            
            h0X <- logistic(x %*% coef.h0)
            h1X <- logistic(x %*% coef.h1)
            gX <- logistic(x %*% coef.g)
            
            ind <- !((treat == 1) & ((y == 0) | (y == (J+1))))
            w <- rep(NA, length(y))
            w[ind] <- gX[ind]*dbinom((y-treat)[ind], size = J, prob = h1X[ind], log = FALSE) /
              (gX[ind]*dbinom((y-treat)[ind], size = J, prob = h1X[ind], log = FALSE) +
               (1-gX[ind])*dbinom(y[ind], size = J, prob = h0X[ind], log = FALSE))
            w[(treat == 1) & (y == 0)] <- 0
            w[(treat == 1) & (y == (J+1))] <- 1
            
            return(w)
          }
          
          ## Mstep 1: weighted MLE for logistic regression
          wlogit.fit.std <- function(y, treat, x, w, par = NULL, wt) {
              ## wt is survey weights
              
            yrep <- rep(c(1,0), each = length(y))
            xrep <- rbind(x, x)
            wrep <- c(w, 1-w)
            wtrep <- c(wt, wt)
            return(glm(cbind(yrep, 1-yrep) ~ xrep - 1, weights = wrep * wtrep, family = binomial(logit), start = par))
            
          }
          
          ## Mstep 2: weighted MLE for beta-binomial regression
          wbb.fit.std <- function(J, y, treat, x, w, par0 = NULL, par1 = NULL) {
            
            Not0 <- ((treat == 1) & (y == (J+1)))
            y0 <- y[!Not0]
            x0 <- x[!Not0,]
            w0 <- 1-w[!Not0]
            fit0 <- vglm(cbind(y0, J-y0) ~ x0, betabinomial, weights = w0, coefstart = par0)
            
            Not1 <- ((treat == 1) & (y == 0))
            y1 <- y
            y1[treat == 1] <- y1[treat == 1] - 1
            y1 <- y1[!Not1]
            x1 <- x[!Not1]
            w1 <- w[!Not1]
            fit1 <- vglm(cbind(y1, J-y1) ~ x1, betabinomial, weights = w1, coefstart = par1)
            
            return(list(fit1 = fit1, fit0 = fit0))
          }
          
          ## Mstep2: weighted MLE for binomial regression
          wbinom.fit.std <- function(J, y, treat, x, w, par0 = NULL, par1 = NULL) {
            
            Not0 <- ((treat == 1) & (y == (J+1)))
            y0 <- y[!Not0]
            x0 <- x[!Not0,]
            w0 <- 1-w[!Not0]
            fit0 <- glm(cbind(y0, J-y0) ~ x0, family = binomial(logit), weights = w0, start = par0)
            
            Not1 <- ((treat == 1) & (y == 0))
            y1 <- y
            y1[treat == 1] <- y1[treat == 1] - 1
            y1 <- y1[!Not1]
            x1 <- x[!Not1,]
            w1 <- w[!Not1]
            fit1 <- glm(cbind(y1, J-y1) ~ x1, family = binomial(logit), weights = w1, start = par1)
            
            return(list(fit1 = fit1, fit0 = fit0))
          }
          
          ##
          ## Running the EM algorithm
          ##
                    
          if(constrained==F) {
            
            ## Run unconstrained model
            
            if (overdispersed==T) {
              par <- c(rep(c(coef.control.start, 0.5), 2), coef.treat.start)
            } else {
              par <- c(rep(coef.control.start, 2), coef.treat.start)
            }
            
            ## max number of iterations
            pllik <- -Inf
            
            if (overdispersed==T) {
              llik <- obs.llik.std(par, J = J, y = y.all, treat = t, x = x.all, wt = w.all)
            } else {
              llik <- obs.llik.binom.std(par, J = J, y = y.all, treat = t, x = x.all, wt = w.all)
            }
            
            Not0 <- (t & (y.all == (J+1)))
            Not1 <- (t & (y.all == 0))
            
            counter <- 0
            while (((llik - pllik) > 10^(-8)) & (counter < maxIter)) {
              
              if(overdispersed==T) {
                
                w <- Estep.std(par, J, y.all, t, x.all)
                
                lfit <- wlogit.fit.std(y.all, t, x.all, w, par = par[(nPar*2+3):length(par)], wt = w.all)
              
                y1 <- y.all
                
                y1[t==1] <- y1[t==1]-1

                if (intercept.only == TRUE) {
                  fit0 <- vglm(cbind(y.all[!Not0], J-y.all[!Not0]) ~ 1,
                               betabinomial, weights = (1-w[!Not0]) * w.all[!Not0], coefstart = par[1:2])
                  fit1 <- vglm(cbind(y1[!Not1], J-y1[!Not1]) ~ 1,
                               betabinomial, weights = w[!Not1] * w.all[!Not1],
                               coefstart = par[3:4])
                  par <- c(coef(fit0), coef(fit1), coef(lfit))
                  
                } else {
                  fit0 <- vglm(cbind(y.all[!Not0], J-y.all[!Not0]) ~ x.all[!Not0, -1, drop = FALSE],
                               betabinomial, weights = (1-w[!Not0]) * w.all[!Not0],
                               coefstart = par[c(1,(nPar+1),2:(nPar))])
                  fit1 <- vglm(cbind(y1[!Not1], J-y1[!Not1]) ~ x.all[!Not1, -1, drop = FALSE] ,
                               betabinomial, weights = w[!Not1] * w.all[!Not1],
                               coefstart = par[c((nPar+2),(2*nPar+2),(nPar+3):(2*nPar+1))])
                  par <- c(coef(fit0)[c(1,3:(nPar+1),2)], coef(fit1)[c(1,3:(nPar+1),2)], coef(lfit))
                }
                                
              } else {
                
                w <- Estep.binom.std(par, J, y.all, t, x.all)

                lfit <- wlogit.fit.std(y.all, t, x.all, w, par = par[(nPar*2+1):length(par)], wt = w.all)
                
                fit0 <- glm(cbind(y.all[!Not0], J-y.all[!Not0]) ~ x.all[!Not0,] - 1,
                            family = binomial(logit), weights = (1-w[!Not0]) * w.all[!Not0], start = par[1:(nPar)])
                
                y1 <- y.all
                
                y1[t==1] <- y1[t==1]-1
                
                fit1 <- glm(cbind(y1[!Not1], J-y1[!Not1]) ~ x.all[!Not1,] - 1,
                            family = binomial(logit), weights = w[!Not1] * w.all[!Not1],
                            start = par[(nPar+1):(2*nPar)])
                
                par <- c(coef(fit0), coef(fit1), coef(lfit))
                
              }
              
              pllik <- llik
              
              if(verbose==T)
                cat(paste(counter, round(llik, 4), "\n"))
              
              if (overdispersed==T) {
                llik <- obs.llik.std(par, J = J, y = y.all, treat = t, x = x.all, wt = w.all)
              } else {
                llik <- obs.llik.binom.std(par, J = J, y = y.all, treat = t, x = x.all, wt = w.all)
              }
              
              counter <- counter + 1
              
              if (llik < pllik) 
                warning("log-likelihood is not monotonically increasing.")
              ## NOTE CHANGED FROM STOP()
              
              if(counter == (maxIter-1))
                warning("number of iterations exceeded maximum in ML.")
              
            }
            
            ## getting  standard errors
            
            if (overdispersed==T) {
              MLEfit <- optim(par, obs.llik.std, method = "BFGS", J = J, y = y.all,
                              treat = t, x = x.all, wt = w.all, hessian = TRUE, control = list(maxit = 0))
              vcov.mle <- solve(-MLEfit$hessian)
              se.mle <- sqrt(diag(vcov.mle))
            } else {
              MLEfit <- optim(par, obs.llik.binom.std, method = "BFGS", J = J, y = y.all,
                              treat = t, x = x.all, wt = w.all, hessian = TRUE, control = list(maxit = 0))
              vcov.mle <- solve(-MLEfit$hessian)
              se.mle <- sqrt(diag(vcov.mle))
            }
            
          } else { # end of unconstrained model
            
            ## constrained model

            if (overdispersed==T) {
              par <- c(coef.control.start, 0.5, coef.treat.start)
            } else {
              par <- c(coef.control.start, coef.treat.start)
            }
            
            pllik.const <- -Inf
            
            if (overdispersed==T) {
              llik.const <- obs.llik.std(par, J = J, y = y.all, treat = t,
                                         x = x.all, wt = w.all, const = TRUE)
            } else {
              llik.const <- obs.llik.binom.std(par, J = J, y = y.all, treat = t,
                                               x = x.all, wt = w.all, const = TRUE)
            }
            
            counter <- 0
            while (((llik.const - pllik.const) > 10^(-8)) & (counter < maxIter)) {
              
              if (overdispersed==T) {
                
                w <- Estep.std(par, J, y.all, t, x.all, const = TRUE)
                lfit <- wlogit.fit.std(y.all, t, x.all, w, par = par[(nPar+2):(nPar*2+1)], wt = w.all)
                
              } else {
                
                w <- Estep.binom.std(par, J, y.all, t, x.all, const = TRUE)
                lfit <- wlogit.fit.std(y.all, t, x.all, w, par = par[(nPar+1):(nPar*2)], wt = w.all)
                
              }
              
              y.var <- as.character(formula)[[2]]
 
              data.all <- as.data.frame(cbind(y.all, x.all))
              names(data.all)[1] <- y.var
              
              dtmp <- rbind(cbind(data.all, w, t, w.all), cbind(data.all, w, t, w.all)[t==1, ])
              
              dtmp[((dtmp$t == 1) & ((1:nrow(dtmp)) <= n)), paste(y.var)] <-
                dtmp[((dtmp$t == 1) & ((1:nrow(dtmp)) <= n)), paste(y.var)] - 1
              
              dtmp$w[((dtmp$t == 1) & ((1:nrow(dtmp)) > n))] <-
                1 - dtmp$w[((dtmp$t == 1) & ((1:nrow(dtmp)) > n))]
              
              dtmp$w[dtmp$t == 0] <- 1

              dtmp$w <- dtmp$w * dtmp$w.all
              
              dtmp <- dtmp[dtmp$w > 0, ]
              
              if (overdispersed==T) {

                if (intercept.only == TRUE) {

                  fit <- vglm(as.formula(paste("cbind(", y.var, ", J-", y.var, ") ~ 1")),
                              betabinomial, weights = dtmp$w,
                              coefstart = par[1:2], data = dtmp)
                  
                  par <- c(coef(fit), coef(lfit))

                } else {
                  
                  fit <- vglm(as.formula(paste("cbind(", y.var, ", J-", y.var, ") ~ ",
                                               paste(x.vars, collapse=" + "))),
                              betabinomial, weights = dtmp$w,
                              coefstart = par[c(1,(nPar+1),2:(nPar))], data = dtmp)
                  
                  par <- c(coef(fit)[c(1,3:(nPar+1),2)], coef(lfit))

                }
                  
              } else {

                ## GB: This is where the error comes in for survey weights.

                fit <- glm(as.formula(paste("cbind(", y.var, ", J-", y.var, ") ~ ",
                                            paste(x.vars, collapse=" + "))),
                           family = binomial(logit), weights = dtmp$w,
                           start = par[1:(nPar)], data = dtmp)
                
                par <- c(coef(fit), coef(lfit))
                
              }
              
              pllik.const <- llik.const
              
              if(verbose==T)
                cat(paste(counter, round(llik.const, 4), "\n"))
              
              if (overdispersed==T) {
                llik.const <- obs.llik.std(par, J = J, y = y.all, treat = t,
                                           x = x.all, wt = w.all, const = TRUE)
              } else {
                llik.const <- obs.llik.binom.std(par, J = J, y = y.all, treat = t,
                                                 x = x.all, wt = w.all, const = TRUE)
              }
              
              counter <- counter + 1
              if (llik.const < pllik.const)
                warning("log-likelihood is not monotonically increasing.")
              ## NOTE CHANGED FROM STOP()
              
              if(counter == (maxIter-1))
                warning("number of iterations exceeded maximum in ML")
              
            }
            
            if (overdispersed==T) {
              
              MLEfit <- optim(par, obs.llik.std, method = "BFGS", J = J,
                              y = y.all, treat = t, x = x.all, wt = w.all, const = TRUE,
                              hessian = TRUE, control = list(maxit = 0))
              
              vcov.mle <- solve(-MLEfit$hessian)
              se.mle <- sqrt(diag(vcov.mle))
              
            } else {
              
              MLEfit <- optim(par, obs.llik.binom.std, method = "BFGS", J = J,
                              y = y.all, treat = t,  x = x.all, wt = w.all, const = TRUE,
                              hessian = TRUE, control = list(maxit = 0))
              
              vcov.mle <- solve(-MLEfit$hessian)
              se.mle <- sqrt(diag(vcov.mle))
              
            }

          } # end of constrained model

        } else if (multi == TRUE) {

          ##
          ## Start multi-treatment standard design model
          
          ## test start values
          ## par <- rep(0, 3*k + 2*(J+1))
            
          obs.llik.multi <- function(par, J, y, treat, x, multi = "none") {

            treat.values <- sort(unique(treat[treat!=0]))
            k <- ncol(x)

            ## pull out all g coefs, including alpha_yt and beta_t
            par.g <- par[(k+1):length(par)]

            ## pull out h coefs and construct h vector
            coef.h <- par[1:k]
            hX <- logistic(x %*% coef.h)

            ## separate coefs for alpha_yt and beta_t and construct
            ## g_t vector for each treatment condition by filling in a
            ## matrix with columns representing possible values of Y_i(0),
            ## rows individuals, and it is in a list where g_t = gX[[t]]
            coef.g.b <- coef.g.a <- gX <- list()
            for (m in 1:length(treat.values)) {

              if (multi == "indicators") {
                coef.g.b[[m]] <- par.g[((m-1)*(k+J)+1):((m-1)*(k+J)+k)]
                coef.g.a[[m]] <- c(par.g[((m-1)*(k+J)+k+1):((m-1)*(k+J)+k+J)],0)
              } else if (multi == "level") {
                coef.g.b[[m]] <- par.g[((m-1)*(k+1)+1):((m-1)*(k+1)+k)]
                coef.g.a[[m]] <- par.g[((m-1)*(k+1)+k+1)]
              } else if (multi == "none") {
                coef.g.b[[m]] <- par.g[((m-1)*k+1):((m-1)*k+k)]
              }

              gX[[m]] <- matrix(NA, nrow = length(y), ncol = J+1)

              if (multi == "indicators") {
                for(y.iter in 0:J) {
                  gX[[m]][,y.iter + 1] <- logistic(coef.g.a[[m]][y.iter + 1] + x %*% coef.g.b[[m]])
                }
              } else if (multi == "level") {
                for(y.iter in 0:J) {
                  gX[[m]][,y.iter + 1] <- logistic(y.iter * coef.g.a[[m]] + x %*% coef.g.b[[m]])
                }
              } else if (multi == "none") {
                for(y.iter in 0:J) {
                  gX[[m]][,y.iter + 1] <- logistic(x %*% coef.g.b[[m]])
                }
              }
            }

            ## contribution from control group
            
            ind0y <- (treat==0)
            p0y <- sum(dbinom(x = y[ind0y], size = J, prob = hX[ind0y], log = TRUE))

            ## contribution from treatment groups
            
            pm0 <- pmJ1 <- rep(NA, length(treat.values))
            pmy <- matrix(NA, nrow = length(treat.values), ncol = J)

            for (m in 1:length(treat.values)) {
              
              indm0 <- ((treat == m) & (y == 0))
              pm0[m] <- sum(dbinom(x = 0, size = J, prob = hX[indm0], log = TRUE) +
                            log(1 - gX[[m]][indm0, 0 + 1]))

              indmJ1 <- ((treat == m) & (y == (J+1)))
              pmJ1[m] <- sum(dbinom(x = J, size = J, prob = hX[indmJ1], log = TRUE) +
                             log(gX[[m]][indmJ1, J + 1]))

              for (y.iter in 1:J) {
                indmy <- ((treat == m) & (y == y.iter))
                pmy[m,y.iter] <- sum(log(gX[[m]][indmy, (y.iter - 1) + 1] *
                                         dbinom(x = y.iter - 1, size = J,
                                                prob = hX[indmy], log = FALSE) +
                 (1-gX[[m]][indmy, y.iter + 1]) * dbinom(x = y.iter, size = J,
                                                         prob = hX[indmy], log = FALSE)))
              }
                         
            }

            return(p0y + sum(pm0) + sum(pmJ1) + sum(pmy))
            
          }

          Estep.multi <- function(par, J, y, treat, x, m, multi = "none") {

            treat.values <- sort(unique(treat[treat!=0]))
            k <- ncol(x)

            ## pull out all g coefs, including alpha_yt and beta_t
            par.g <- par[(k+1):length(par)]

            ## pull out h coefs and construct h vector
            coef.h <- par[1:k]
            hX <- logistic(x %*% coef.h)

            if (multi == "indicators") {
              coef.g.b <- par.g[((m-1)*(k+J)+1):((m-1)*(k+J)+k)]
              coef.g.a <- c(par.g[((m-1)*(k+J)+k+1):((m-1)*(k+J)+k+J)],0)
            } else if (multi == "level") {
              coef.g.b <- par.g[((m-1)*(k+1)+1):((m-1)*(k+1)+k)]
              coef.g.a <- par.g[((m-1)*(k+1)+k+1)]
            } else if (multi == "none") {
              coef.g.b <- par.g[((m-1)*k+1):((m-1)*k+k)]
            }

            ## construct g_t matrix, where columns represent values of y_i(0),
            ## since g_t varies in y_i(0) through alpha_ty. Rows are indivs.

            ## in constrained model, g_t is constant across y_i(0) values
            
            gX <- matrix(NA, nrow = length(y), ncol = J+1)

            if (multi == "indicators") {
              for(y.iter in 0:J) {
                gX[,y.iter + 1] <- logistic(coef.g.a[y.iter + 1] + x %*% coef.g.b)
              }
            } else if (multi == "level") {
              for(y.iter in 0:J) {
                gX[,y.iter + 1] <- logistic(y.iter * coef.g.a + x %*% coef.g.b)
              }
            } else if (multi == "none") {
              for(y.iter in 0:J) {
                gX[,y.iter + 1] <- logistic(x %*% coef.g.b)
              }
            }

            w <- rep(NA, length(y))

            ind0 <- ((y==0) & (treat==m))
            w[ind0] <- 0

            indJ1 <- ((y==(J+1)) & (treat==m))
            w[indJ1] <- 1

            ## for values of Y_i(0) not 0 or (J+1), construct the weight by
            ## looping through all possible values of Y_i(0) and adding weights
            ## for individuals with that Y_i(0)
            
            for (y.iter in 1:J) {
              
              indother <- ((ind0==FALSE) & (indJ1==FALSE) & (y==y.iter) & (treat==m))
              w[indother] <- (gX[indother, (y.iter - 1) + 1] *
                              dbinom(x = y.iter - 1, size = J, prob = hX[indother],
                                     log = FALSE)) /
                                (gX[indother, (y.iter - 1) + 1] *
                                 dbinom(x = y.iter - 1, size = J, prob = hX[indother],
                                        log = FALSE) +
                                 (1 - gX[indother, y.iter + 1]) *
                                 dbinom(x = y.iter, size = J, prob = hX[indother],
                                        log = FALSE) )
              
            }
            
            return(w)
            
          }

          ##
          ## start EM algorithm

          ## get starting values from NLS
          
          par <- par.control.nls.std
          if (multi.condition == "none") {
            par <- c(par, do.call(c, par.treat.nls.std))
          } else if (multi.condition == "indicators") {
            for (m in 1:length(treatment.values)){
              par <- c(par, par.treat.nls.std[[m]], rep(0, J))
            }
          } else if (multi.condition == "level") {
            for (m in 1:length(treatment.values)){
              par <- c(par, par.treat.nls.std[[m]], 0)
            }
          }
           
          pllik <- -Inf
          
          llik <- obs.llik.multi(par, J = J, y = y.all, treat = t, x = x.all,
                                 multi = multi.condition)

          ##
          ## create temporary datasets for M step

          ## start fit data for H including control group data
          
          ytmpH <- y.all[t==0]
          xtmpH <- x.all[t==0, , drop = FALSE]
          
          ## now loop over treatment conditions to construct data for g_t fit
          ## and add treatment condition data (m times) to data for h fit

          ytmpG <- xtmpG <- dvtmpG <- list()

          for (m in 1:length(treatment.values)){

            ##
            ## set up temporary data for g fit
            
            ytmpG[[m]] <- c(y.all[t==m]-1, y.all[t==m])
            xtmpG[[m]] <- rbind(x.all[t==m, , drop = FALSE], x.all[t==m, , drop = FALSE])
            dvtmpG[[m]] <- c(rep(1, sum(t==m)), rep(0, sum(t==m)))

            if (multi.condition == "indicators") {
              
              ## create matrix of indicator variables for non-sensitive y_i(0) value
              ## omits indicator for value == J              
              y.ind <- matrix(NA, nrow = sum(t==m) * 2, ncol = J)
              for (y.iter in 0:(J-1))
                y.ind[, y.iter + 1] <- ytmpG[[m]] == y.iter
              
              ## x vector for g_t fit includes covariates and y_i(0) indicator matrix
              ## the constrained model does not include these indicators in the g_t M fits              
              xtmpG[[m]] <- cbind(xtmpG[[m]], y.ind)
              
            } else if (multi.condition == "level") {

              xtmpG[[m]] <- cbind(xtmpG[[m]], ytmpG[[m]])
              
            }
            
            ##
            ## set up temporary data for h fit
            
            ytmpH <- c(ytmpH, y.all[t==m]-1, y.all[t==m])
            xtmpH <- rbind(xtmpH, x.all[t==m, , drop = FALSE], x.all[t==m, , drop = FALSE])
            
          }

          counter <- 0
          while (((llik - pllik) > 10^(-8)) & (counter < maxIter)) {
            
            wtmpH <- rep(1, sum(t==0))
            
            lfit <- list()
            
            for (m in 1:length(treatment.values)){

              ## get E step values for this treatment condition
              w <- Estep.multi(par, J, y.all, t, x.all, m, multi = multi.condition)

              ## create weight vectors for temporary data for fits

              wtmpG <- c(w[t==m], 1-w[t==m])

              ## for h, we are adding weights to the vector of weights = 1 created earlier
              
              wtmpH <- c(wtmpH, w[t==m], 1-w[t==m])

              if (multi.condition == "none")
                par.start <- par[(nPar + (m-1) * nPar + 1):(nPar + (m-1) * nPar + nPar)]
              else if (multi.condition == "indicators")
                par.start <- par[(nPar + (m - 1) * (nPar + J) + 1):
                                 (nPar + (m - 1) * (nPar + J) + nPar + J)]
              else if (multi.condition == "level")
                par.start <- par[(nPar + (m - 1) * (nPar + 1) + 1):
                                 (nPar + (m - 1) * (nPar + 1) + nPar + 1)]

              ## run g_t fits
              lfit[[m]] <- glm(cbind(dvtmpG[[m]][wtmpG>0], 1 - dvtmpG[[m]][wtmpG>0]) ~
                               xtmpG[[m]][wtmpG>0, , drop = FALSE] - 1, weights = wtmpG[wtmpG>0],
                               family = binomial(logit), start = par.start,
                               control = glm.control(maxit = maxIter))

            }

            ## run h fit
            fit <- glm(cbind(ytmpH[wtmpH>0], J - ytmpH[wtmpH>0]) ~ xtmpH[wtmpH>0, , drop = FALSE] - 1,
                       family = binomial(logit), weights = wtmpH[wtmpH>0], start = par[1:(nPar)])

            ## put coefficients back together
            par <- coef(fit)
            par.treat <- list()
            for (m in 1:length(treatment.values)) {
               par <- c(par, coef(lfit[[m]]))
               par.treat[[m]] <- coef(lfit[[m]])
            }

            pllik <- llik
            
            if(verbose==T)
              cat(paste(counter, round(llik, 4), "\n"))
            
            llik <- obs.llik.multi(par, J = J, y = y.all, treat = t, x = x.all,
                                   multi = multi.condition)
            
            counter <- counter + 1
            if (llik < pllik)
              warning("log-likelihood is not monotonically increasing.")
            ## NOTE CHANGED FROM STOP()
            
            if(counter == (maxIter-1))
              warning("number of iterations exceeded maximum in ML")
            
          }

          MLEfit <- optim(par, obs.llik.multi, method = "BFGS", J = J, y = y.all, treat = t,
                          x = x.all, hessian = TRUE, control = list(maxit = 0),
                          multi = multi.condition)
          
          vcov.mle <- solve(-MLEfit$hessian)
          
          se.mle <- sqrt(diag(vcov.mle))
  
        } # end of multi treat regression
        
      } else if (boundary == TRUE) {
                
        coef.control.start <- par.control.nls.std
        coef.treat.start <- par.treat.nls.std

        ##
        ##  Observed data log-likelihood for binomial
        ##
        
        obs.llik.binom.boundary <- function(par, J, y, treat, x, x.floor, x.ceiling, ceiling.fit,
                                            floor.fit, ceiling, floor) {
          
          k <- ncol(x)
          k.ceiling <- ncol(x.ceiling)
          k.floor <- ncol(x.floor)
          
          coef.h <- par[1:k]
          coef.g <- par[(k+1):(2*k)]
          coef.ql <- par[(2*k+1):(2*k + k.floor)] 
          coef.qu <- par[(2*k+k.floor+1):(2*k+k.floor+k.ceiling)] 
          
          hX <- logistic(x %*% coef.h)
          gX <- logistic(x %*% coef.g)
          quX <- logistic(x.ceiling %*% coef.qu)
          qlX <- logistic(x.floor %*% coef.ql)
          
          ind0y <- ((treat == 0) & (y < (J+1)))
          ind10 <- ((treat == 1) & (y == 0))
          ind11 <- ((treat == 1) & (y == 1))
          ind1y <- ((treat == 1) & (y > 1) & (y < J))
          ind1J <- ((treat == 1) & (y == J))
          ind1J1 <- ((treat == 1) & (y == (J+1)))
          
          p0y <- p10 <- p11 <- p1y <- p1J <- p1J1 <- 0
          
          if(sum(ind0y) > 0){
            p0y <- sum(dbinom(x = y[ind0y], size = J, prob = hX[ind0y], log = TRUE))
          }
          
          if(sum(ind10) > 0){
            p10 <- sum(log((1-gX[ind10]) + gX[ind10] * qlX[ind10])
                       + dbinom(x = 0, size = J, prob = hX[ind10], log = TRUE))
          }
          
          if (sum(ind11) > 0) {
            p11 <- sum(log( gX[ind11] * (1-qlX[ind11]) *
                           dbinom(x = 0, size = J, prob = hX[ind11], log = FALSE) +
                           (1-gX[ind11]) * dbinom(x = 1, size = J, prob = hX[ind11], log = FALSE) ))
          }
          
          if (sum(ind1y) > 0) {
            p1y <- sum(log( gX[ind1y] * dbinom(x = y[ind1y]-1, size = J, prob = hX[ind1y], log = FALSE) +
                           (1-gX[ind1y]) * dbinom(x = y[ind1y], size = J, prob = hX[ind1y], log = FALSE) ))
          }
          
          if (sum(ind1J)) {
            p1J <- sum(log( gX[ind1J] * (dbinom(x = J-1, size = J, prob = hX[ind1J], log = FALSE) +
                                         quX[ind1J] * dbinom(x = J, size = J, prob = hX[ind1J], log = FALSE)) +
                           (1-gX[ind1J]) * dbinom(x = J, size = J, prob = hX[ind1J], log = FALSE) ))
          }
          
          if (sum(ind1J1) > 0) {
            p1J1 <- sum(log(1-quX[ind1J1]) + log(gX[ind1J1]) + dbinom(x = J, size = J, prob = hX[ind1J1], log = TRUE))
          }
          
          if (ceiling.fit=="bayesglm" & ceiling == TRUE) {
            p.prior.ceiling <- sum(dcauchy(x = coef.qu,
                                           scale = rep(2.5, length(coef.qu)), log = TRUE))
          } else {
            p.prior.ceiling <- 0
          }
          
          if (floor.fit=="bayesglm" & floor == TRUE) {
            p.prior.floor <- sum(dcauchy(x = coef.ql,
                                         scale = rep(2.5, length(coef.ql)), log = TRUE))
          } else {
            p.prior.floor <- 0
          }
					
          return(p0y  + p10  + p11  + p1y  + p1J  + p1J1 + p.prior.ceiling + p.prior.floor)
        }
				
        ##
        ##  Observed data log-likelihood for binomial for optim -- shortened par
        ##
        obs.llik.binom.optim.boundary <- function(par, J, y, treat, x, x.floor, x.ceiling,
                                                  ceiling.fit, floor.fit, floor, ceiling) {
          
          k <- ncol(x)
          k.ceiling <- ncol(x.ceiling)
          k.floor <- ncol(x.floor)
          coef.h <- par[1:k]
          coef.g <- par[(k+1):(2*k)]
          if(floor==TRUE) {
            coef.ql <- par[(2*k+1):(2*k + k.floor)]
            if(ceiling==TRUE){
              coef.qu <- par[(2*k+k.floor+1):(2*k+k.floor+k.ceiling)] 
            } else {
              coef.qu <- c(-Inf, rep(0, (k.ceiling-1)))
            }
          } else {
            coef.ql <- c(-Inf, rep(0, (k.floor-1)))
            if(ceiling==TRUE){
              coef.qu <- par[(2*k+1):(2*k + k.ceiling)] 
            } else {
              coef.qu <- c(-Inf, rep(0, (k.ceiling-1)))
            }
          }
          
          hX <- logistic(x %*% coef.h)
          gX <- logistic(x %*% coef.g)
          quX <- logistic(x.ceiling %*% coef.qu)
          qlX <- logistic(x.floor %*% coef.ql)
          
          ind0y <- ((treat == 0) & (y < (J+1)))
          ind10 <- ((treat == 1) & (y == 0))
          ind11 <- ((treat == 1) & (y == 1))
          ind1y <- ((treat == 1) & (y > 1) & (y < J))
          ind1J <- ((treat == 1) & (y == J))
          ind1J1 <- ((treat == 1) & (y == (J+1)))
          
          p0y <- p10 <- p11 <- p1y <- p1J <- p1J1 <- 0
          
          if(sum(ind0y) > 0){
            p0y <- sum(dbinom(x = y[ind0y], size = J, prob = hX[ind0y], log = TRUE))
          }
          
          if(sum(ind10) > 0){
            p10 <- sum(log((1-gX[ind10]) + gX[ind10] * qlX[ind10]) + dbinom(x = 0, size = J, prob = hX[ind10], log = TRUE))
          }
          
          if(sum(ind11) > 0){
            p11 <- sum(log( gX[ind11] * (1-qlX[ind11]) * dbinom(x = 0, size = J, prob = hX[ind11], log = FALSE) + (1-gX[ind11]) * dbinom(x = 1, size = J, prob = hX[ind11], log = FALSE) ))
          }
          
          if(sum(ind1y) > 0){
            p1y <- sum(log( gX[ind1y] * dbinom(x = y[ind1y]-1, size = J, prob = hX[ind1y], log = FALSE) + (1-gX[ind1y]) * dbinom(x = y[ind1y], size = J, prob = hX[ind1y], log = FALSE) ))
          }
          
          if(sum(ind1J)){
            p1J <- sum(log( gX[ind1J] * (dbinom(x = J-1, size = J, prob = hX[ind1J], log = FALSE) + quX[ind1J] * dbinom(x = J, size = J, prob = hX[ind1J], log = FALSE)) + (1-gX[ind1J]) * dbinom(x = J, size = J, prob = hX[ind1J], log = FALSE) ))
          }
          
          if(sum(ind1J1) > 0){
            p1J1 <- sum(log(1-quX[ind1J1]) + log(gX[ind1J1]) + dbinom(x = J, size = J, prob = hX[ind1J1], log = TRUE))
          }
          
          if(ceiling.fit=="bayesglm" & ceiling==TRUE) {
            p.prior.ceiling <- sum(dcauchy(x = coef.qu,
                                         scale = rep(2.5, length(coef.qu)), log = TRUE))
          } else {
            p.prior.ceiling <- 0
          }
          
          if(floor.fit=="bayesglm" & floor==TRUE) {
            p.prior.floor <- sum(dcauchy(x = coef.ql,
                                         scale = rep(2.5, length(coef.ql)), log = TRUE))
          } else {
            p.prior.floor <- 0
          }
          
          return(p0y  + p10  + p11  + p1y  + p1J  + p1J1 + p.prior.ceiling + p.prior.floor)
          
        }
        
        ##
        ## EM algorithm
        ##
        
        ## Estep for beta-binomial
        Estep.boundary <- function(par, J, y, treat, x, product) {

          ## not currently used
          ## needs to be updated for various changes, including ceiling/floor formulaes
          
          k <- ncol(x)
          coef.h <- par[1:k]
          rho.h <- logistic(par[(k+1)])
          coef.g <- par[(k+2):(2*k+1)]
          coef.ql <- par[(2*k+2):(3*k+1)] 
          coef.qu <- par[(3*k+2):(4*k+1)] 
          
          hX <- logistic(x %*% coef.h)
          gX <- logistic(x %*% coef.g)
          quX <- logistic(x %*% coef.qu)
          qlX <- logistic(x %*% coef.ql)
          
          ind0 <- (y==0)
          ind1 <-  (y==1)
          indJ <- (y==J)
          indJ1 <- (y==(J+1))
          indother <- (y > 1 & y < J)
          
          w <- rep(NA, length(y))
          
          if(product == FALSE){
            
            w[ind0] <- dBB(0, bd = J, mu = hX[ind0], sigma = rho.h/(1-rho.h), log = FALSE)*gX[ind0]*qlX[ind0] /
              (dBB(0, bd = J, mu = hX[ind0], sigma = rho.h/(1-rho.h), log = FALSE) * (gX[ind0]*qlX[ind0] + (1-gX[ind0]) ) )
            
            w[ind1] <- dBB(0, bd = J, mu = hX[ind1], sigma = rho.h/(1-rho.h), log = FALSE)*gX[ind1]*(1-qlX[ind1]) / ( dBB(0, bd = J, mu = hX[ind1], sigma = rho.h/(1-rho.h), log = FALSE) * gX[ind1] * (1-qlX[ind1]) + dBB(1, bd = J, mu = hX[ind1], sigma = rho.h/(1-rho.h), log = FALSE)*(1-gX[ind1]) )
            
            w[indother] <- gX[indother]*dBB((y-1)[indother], bd = J, mu = hX[indother], sigma = rho.h/(1-rho.h), log = FALSE) /
              (gX[indother]*dBB((y-1)[indother], bd = J, mu = hX[indother], sigma = rho.h/(1-rho.h), log = FALSE) +
               (1-gX[indother])*dBB(y[indother], bd = J, mu = hX[indother], sigma = rho.h/(1-rho.h), log = FALSE))
            
            w[indJ] <- gX[indJ]*( dBB(J-1, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE) + quX[indJ]*dBB(J, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE) ) / ( gX[indJ] * ( dBB(J-1, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE) + quX[indJ] * dBB(J, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE)) + (1-gX[indJ])*dBB(J, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE) )
            
            w[indJ1] <- 1
            
          }
          
          if(product == TRUE){
            
            w[ind0] <- 0
            
            w[ind1] <- dBB(0, bd = J, mu = hX[ind1], sigma = rho.h/(1-rho.h), log = FALSE)*gX[ind1]*(1-qlX[ind1]) / ( dBB(0, bd = J, mu = hX[ind1], sigma = rho.h/(1-rho.h), log = FALSE) * gX[ind1] * (1-qlX[ind1]) + dBB(1, bd = J, mu = hX[ind1], sigma = rho.h/(1-rho.h), log = FALSE)*(1-gX[ind1]) )
            
            w[indother] <- gX[indother]*dBB((y-1)[indother], bd = J, mu = hX[indother], sigma = rho.h/(1-rho.h), log = FALSE) / (gX[indother]*dBB((y-1)[indother], bd = J, mu = hX[indother], sigma = rho.h/(1-rho.h), log = FALSE) + (1-gX[indother])*dBB(y[indother], bd = J, mu = hX[indother], sigma = rho.h/(1-rho.h), log = FALSE))
            
            w[indJ] <- gX[indJ]*dBB(J - 1, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE) / (gX[indJ]*(dBB(J - 1, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE) + quX[indJ]*dBB(J, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE)) + (1-gX[indJ])*dBB(J, bd = J, mu = hX[indJ], sigma = rho.h/(1-rho.h), log = FALSE))
            
            w[indJ1] <- 1
            
          }
          
          return(w)
          
        }
        
        ## Estep for binomial
        Estep.binom.boundary <- function(par, J, y, treat, x, x.floor, x.ceiling, product) {
          
          k <- ncol(x)
          k.ceiling <- ncol(x.ceiling)
          k.floor <- ncol(x.floor)
          coef.h <- par[1:k]
          coef.g <- par[(k+1):(2*k)]
          coef.ql <- par[(2*k+1):(2*k+k.floor)] 
          coef.qu <- par[(2*k+k.floor+1):(2*k+k.floor+k.ceiling)] 
          
          hX <- logistic(x %*% coef.h)
          gX <- logistic(x %*% coef.g)
          quX <- logistic(x.ceiling %*% coef.qu)
          qlX <- logistic(x.floor %*% coef.ql)
          
          ind0 <- (y==0)
          ind1 <-  (y==1)
          indJ <- (y==J)
          indJ1 <- (y==(J+1))
          indother <- (y > 1 & y < J)
          
          w <- rep(NA, length(y))
          
          if(product == FALSE){
            
            w[ind0] <- (dbinom(x = 0, size = J, prob = hX[ind0], log = FALSE) * gX[ind0] * qlX[ind0]) / (dbinom(x = 0, size = J, prob = hX[ind0], log = FALSE) * (gX[ind0] * qlX[ind0] + (1-gX[ind0])))
            
            w[ind1] <- (dbinom(x = 0, size = J, prob = hX[ind1], log = FALSE) * gX[ind1] * (1-qlX[ind1])) / (dbinom(x = 0, size = J, prob = hX[ind1], log = FALSE) * gX[ind1] * (1-qlX[ind1]) + dbinom(x = 1, size = J, prob = hX[ind1], log = FALSE) * (1-gX[ind1]))
            
            w[indother] <- (gX[indother] * dbinom(x = y[indother]-1, size = J, prob = hX[indother], log = FALSE)) / (gX[indother] * dbinom(x = y[indother]-1, size = J, prob = hX[indother], log = FALSE) + (1-gX[indother]) * dbinom(x = y[indother], size = J, prob = hX[indother], log = FALSE) )
            
            w[indJ] <- (gX[indJ] * (dbinom(x = J-1, size = J, prob = hX[indJ], log = FALSE) + quX[indJ] * dbinom(x = J, size = J, prob = hX[indJ], log = FALSE))) / (gX[indJ] * (dbinom(x = J-1, size = J, prob = hX[indJ], log = FALSE) + quX[indJ] * dbinom(x = J, size = J, prob = hX[indJ], log = FALSE)) + (1-gX[indJ]) * dbinom(x = J, size = J, prob = hX[indJ], log = FALSE))
            
            w[indJ1] <- 1
            
          }
          
          if(product == TRUE){
            
            w[ind0] <- 0
            
            w[ind1] <- (dbinom(x = 0, size = J, prob = hX[ind1], log = FALSE) * gX[ind1] * (1-qlX[ind1])) / (dbinom(x = 0, size = J, prob = hX[ind1], log = FALSE) * gX[ind1] * (1-qlX[ind1]) + dbinom(x = 1, size = J, prob = hX[ind1], log = FALSE) * (1-gX[ind1]))
            
            w[indother] <- (gX[indother] * dbinom(x = y[indother]-1, size = J, prob = hX[indother], log = FALSE)) / (gX[indother] * dbinom(x = y[indother]-1, size = J, prob = hX[indother], log = FALSE) + (1-gX[indother]) * dbinom(x = y[indother], size = J, prob = hX[indother], log = FALSE) )
            
            w[indJ] <- (gX[indJ] * dbinom(x = J-1, size = J, prob = hX[indJ], log = FALSE)) / (gX[indJ] * (dbinom(x = J-1, size = J, prob = hX[indJ], log = FALSE) + quX[indJ] * dbinom(x = J, size = J, prob = hX[indJ], log = FALSE)) + (1-gX[indJ]) * dbinom(x = J, size = J, prob = hX[indJ], log = FALSE))
        
            w[indJ1] <- 1
            
          }
          
          return(w)
          
        }
        
        ## Mstep 1: weighted MLE for logistic regression
        wlogit.fit.boundary <- function(y, x, w, par = NULL, maxIter = 50000) {
          yrep <- rep(c(1,0), each = length(y))
          xrep <- rbind(x, x)
          wrep <- c(w, 1-w)
          return(glm(cbind(yrep, 1-yrep) ~ xrep - 1, weights = wrep, family = binomial(logit), start = par, control = glm.control(maxit = maxIter)))
          
        }
        
        ## Mstep2: weighted MLE for binomial regression
        wbinom.fit.boundary <- function(J, y, treat, x, w, par0 = NULL, par1 = NULL) {
          
          Not0 <- ((treat == 1) & (y == (J+1)))
          y0 <- y[!Not0]
          x0 <- x[!Not0,]
          w0 <- 1-w[!Not0]
          fit0 <- glm(cbind(y0, J-y0) ~ x0, family = binomial(logit), weights = w0, start = par0)
          
          Not1 <- ((treat == 1) & (y == 0))
          y1 <- y
          y1[treat == 1] <- y1[treat == 1] - 1
          y1 <- y1[!Not1]
          x1 <- x[!Not1,]
          w1 <- w[!Not1]
          fit1 <- glm(cbind(y1, J-y1) ~ x1, family = binomial(logit), weights = w1, start = par1)
          
          return(list(fit1 = fit1, fit0 = fit0))
        }
        
        ##
        ## Running the EM algorithm
        ##

        if (overdispersed == T)
          stop("The ceiling-floor liar model does not support overdispersion. Set overdispersed = F.")
                
        y.var <- as.character(formula)[[2]]

        x.ceiling.all <- model.matrix(as.formula(ceiling.formula), as.data.frame(x.all))
        x.ceiling.treatment <- model.matrix(as.formula(ceiling.formula), as.data.frame(x.treatment))
        x.ceiling.control <- model.matrix(as.formula(ceiling.formula), as.data.frame(x.control))
        
        nPar.ceiling <- ncol(x.ceiling.all)
        
        intercept.only.ceiling <- ncol(x.ceiling.all)==1 & sum(x.ceiling.all[,1]==1) == n
        
        coef.names.ceiling <- colnames(x.ceiling.all)
        
        if (intercept.only.ceiling == TRUE) {
          x.vars.ceiling <- "1"
        } else {
          x.vars.ceiling <- coef.names.ceiling[-1]
        }
        
        x.floor.all <- model.matrix(as.formula(floor.formula), as.data.frame(x.all))
        x.floor.treatment <- model.matrix(as.formula(floor.formula), as.data.frame(x.treatment))
        x.floor.control <- model.matrix(as.formula(floor.formula), as.data.frame(x.control))
        nPar.floor <- ncol(x.floor.all)
        
        intercept.only.floor <- ncol(x.floor.all)==1 & sum(x.floor.all[,1]==1) == n
        
        coef.names.floor <- colnames(x.floor.all)
        
        if (intercept.only.floor == TRUE) {
          x.vars.floor <- "1"
        } else {
          x.vars.floor <- coef.names.floor[-1]
        }
        
        data.all <- as.data.frame(cbind(y.all, x.all))
        names(data.all)[1] <- y.var
        
        data.treatment <- as.data.frame(cbind(y.treatment, x.treatment))
        names(data.treatment)[1] <- y.var
        
        ## set start values
        
        data.ceiling.all <- as.data.frame(cbind(y.all, x.ceiling.all))
        names(data.ceiling.all)[1] <- y.var
       
        if (ceiling == TRUE){
    
          data.ceiling.treatment <- as.data.frame(cbind(y.treatment, x.ceiling.treatment))
          names(data.ceiling.treatment)[1] <- y.var
          
          dtmpS <- data.treatment[data.ceiling.treatment[,c(y.var)]==J |
                                  data.ceiling.treatment[,c(y.var)]==(J+1), ]
          
          dtmpS$dv <- dtmpS[,c(y.var)] == J
          
          coef.ceiling.start <- coef(glm(as.formula(paste("cbind(dv, 1-dv) ~ ",
                                                          paste(x.vars.ceiling, collapse=" + "))),
                                         family = binomial(logit), data = dtmpS,
                                         control = glm.control(maxit = maxIter)))
        } else {
          coef.ceiling.start <- c(-Inf, rep(0, (nPar.ceiling-1)))
        }

        data.floor.all <- as.data.frame(cbind(y.all, x.floor.all))
        names(data.floor.all)[1] <- y.var
        
        if (floor == TRUE){
          
          data.floor.treatment <- as.data.frame(cbind(y.treatment, x.floor.treatment))
          names(data.floor.treatment)[1] <- y.var
          
          dtmpS <- data.floor.treatment[data.floor.treatment[,c(y.var)]==0 |
                                        data.floor.treatment[,c(y.var)]==1, ]
          
          dtmpS$dv <- dtmpS[,c(y.var)] == 1
          
          coef.floor.start <- coef(glm(as.formula(paste("cbind(dv, 1-dv) ~ ",
                                                        paste(x.vars.floor, collapse=" + "))),
                                       family = binomial(logit), data = dtmpS,
                                       control = glm.control(maxit = maxIter)))
          
        } else {
          coef.floor.start <- c(-Inf, rep(0, (nPar.floor-1)))
        }
        
        par <- c(coef.control.start, coef.treat.start, coef.floor.start, coef.ceiling.start)
       
        ## calculate starting log likelihood

        pllik <- -Inf
        
        llik <- obs.llik.binom.boundary(par, J = J, y = y.all, treat = t, x = x.all,
                                        x.ceiling = x.ceiling.all, x.floor = x.floor.all,
                                        ceiling.fit = ceiling.fit, floor.fit = floor.fit,
                                        ceiling = ceiling, floor = floor)
        
        ## begin EM iterations
        
        counter <- 0
        while (((llik - pllik) > 10^(-8)) & (counter < maxIter)) {
          
          w <- Estep.binom.boundary(par, J, y.all, t, x.all,
                                    x.floor.all, x.ceiling.all, product = FALSE)
          w.prod <- Estep.binom.boundary(par, J, y.all, t, x.all,
                                         x.floor.all, x.ceiling.all, product = TRUE)
          
          lfit <- wlogit.fit.boundary(y.all[t==1], x.all[t==1, , drop = FALSE], w[t==1],
                             par = par[(nPar+1):(nPar*2)], maxIter = maxIter)
          
          dtmp <- rbind(data.all[t==0, ], data.all[t==1, ], data.all[t==1, ], data.all[t==1, ])
          
          dtmp$w <- c(rep(1, sum(t==0)), (w-w.prod)[t==1], (1-w)[t==1], w.prod[t==1])
          
          dtmp[(sum(c(t==0, t==1, t==1))+1):nrow(dtmp),paste(y.var)] <-
            dtmp[(sum(c(t==0, t==1, t==1))+1):nrow(dtmp),paste(y.var)] - 1
          
          dtmp <- dtmp[dtmp$w > 0, ]
          
          ## temporary data for ceiling logit
          dtmpC <- rbind(data.ceiling.all[t==1 & y.all==(J+1),], data.ceiling.all[t==1 & y.all==J, ])
          
          dtmpC[,paste(y.var)] <- c(rep(0,sum(t==1 & y.all==(J+1))), rep(1,sum(t==1 & y.all==J)))
          
          dtmpC$w <- c( w.prod[t==1 & y.all==(J+1)], (w - w.prod)[t==1 & y.all==J])
          
          dtmpC <- dtmpC[dtmpC$w > 0, ]
          
          ## temporary data for floor logit
          dtmpF <- rbind(data.floor.all[t==1 & y.all==1,], data.floor.all[t==1 & y.all==0, ])
          
          dtmpF[,paste(y.var)] <- c(rep(0,sum(t==1 & y.all==1)), rep(1,sum(t==1 & y.all==0)))
          
          dtmpF$w <- c(w.prod[t==1 & y.all==1], (w-w.prod)[t==1 & y.all==0])
          
          dtmpF <- dtmpF[dtmpF$w > 0, ]
                      
          fit <- glm(as.formula(paste("cbind(", y.var, ", J-", y.var, ") ~ ",
                                      paste(x.vars, collapse=" + "))),
                     family = binomial(logit), weights = dtmp$w, start = par[1:(nPar)],
                     data = dtmp)
          
          if (ceiling==TRUE) {
            ##coef.qufit.start <- par[(3*nPar+1):(4*nPar)]
            coef.qufit.start <- par[(2*nPar + nPar.floor + 1):(2*nPar + nPar.floor + nPar.ceiling)]
          } else {
            coef.qufit.start <- rep(0, nPar.ceiling)
          }

          if (ceiling == TRUE) {
            if (ceiling.fit=="glm") {
              qufit <- glm(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ ",
                                            paste(x.vars.ceiling, collapse=" + "))),
                           weights = dtmpC$w, family = binomial(logit),
                           start = coef.qufit.start, data = dtmpC,
                           control = glm.control(maxit = maxIter))
            } else if(ceiling.fit=="bayesglm") {
              
              if (intercept.only.ceiling == F) {
                qufit <- bayesglm.internal(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ ",
                                                   paste(x.vars.ceiling, collapse=" + "))),
                                  weights = dtmpC$w, family = binomial(logit),
                                  start = coef.qufit.start, data = dtmpC,
                                  control = glm.control(maxit = maxIter), scaled = F)
                
              } else {
                qufit <- bayesglm.internal(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ 1")),
                                  weights = dtmpC$w, family = binomial(logit),
                                  start = coef.qufit.start, data = dtmpC,
                                  control = glm.control(maxit = maxIter), scaled = F)
                
              }
            }
          }
          
          if(floor==TRUE) {
            ##coef.qlfit.start <- par[(2*nPar+1):(3*nPar)]
            coef.qlfit.start <- par[(2*nPar+1):(2*nPar + nPar.floor)]
          } else {
            coef.qlfit.start <- rep(0, nPar.floor)
          }

          if (floor == TRUE) {
            if(floor.fit=="glm") {
              qlfit <- glm(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ ",
                                            paste(x.vars.floor, collapse=" + "))), weights = dtmpF$w,
                           family = binomial(logit), start = coef.qlfit.start, data = dtmpF,
                           control = glm.control(maxit = maxIter))
            } else if(floor.fit=="bayesglm") {
              if (intercept.only.floor == F) {
                qlfit <- bayesglm.internal(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ ",
                                                   paste(x.vars.floor, collapse=" + "))),
                                  weights = dtmpF$w, family = binomial(logit),
                                  start = coef.qlfit.start, data = dtmpF,
                                  control = glm.control(maxit = maxIter), scaled = F)
              } else {
                qlfit <- bayesglm.internal(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ 1")),
                                  weights = dtmpF$w, family = binomial(logit),
                                  start = coef.qlfit.start, data = dtmpF,
                                  control = glm.control(maxit = maxIter), scaled = F)
              }
            }
          }
          
          ## reset parameters based on floor / ceiling options
          
          if(floor==TRUE & ceiling==TRUE) {
            par <- c(coef(fit), coef(lfit),coef(qlfit), coef(qufit))
          } else if(floor==FALSE & ceiling==TRUE) {
            par <- c(coef(fit), coef(lfit),  c(-Inf, rep(0, (nPar.floor-1))), coef(qufit))
          } else if(floor==TRUE & ceiling==FALSE) {
            par <- c(coef(fit), coef(lfit),coef(qlfit),  c(-Inf, rep(0, (nPar.ceiling-1))))
          }
          
          pllik <- llik
          
          if(verbose==T) cat(paste(counter, round(llik, 4), "\n"))
          
          llik <- obs.llik.binom.boundary(par, J = J, y = y.all, treat = t, x = x.all,
                                          x.ceiling = x.ceiling.all, x.floor = x.floor.all,
                                          ceiling.fit = ceiling.fit, floor.fit = floor.fit,
                                          ceiling = ceiling, floor = floor)
          			  
          counter <- counter + 1
          if (llik < pllik)
            warning("log-likelihood is not monotonically increasing.")
          ## NOTE CHANGED FROM STOP()
          
          if(counter == (maxIter-1))
            warning("number of iterations exceeded maximum in ML")
				  				  
        }
          
        if(floor==FALSE & ceiling==TRUE) {
          par <- par[c((1:(2*nPar)), (((2*nPar + nPar.floor)+1):(2*nPar + nPar.floor + nPar.ceiling)))]
        } else if(floor==TRUE & ceiling==FALSE) {
          par <- par[1:(2*nPar + nPar.floor)]
        }
        
        MLEfit <- optim(par, obs.llik.binom.optim.boundary, method = "BFGS", J = J, y = y.all,
                        treat = t,  x = x.all, x.ceiling = x.ceiling.all, x.floor = x.floor.all,
                        ceiling.fit = ceiling.fit,
                        floor.fit  = floor.fit, floor = floor, ceiling = ceiling,
                        hessian = TRUE, control = list(maxit = 0))
        
        if(ceiling.fit=="bayesglm" & ceiling==TRUE) {
          p.prior.ceiling <- sum(dcauchy(x = coef(qufit),
                                         scale = rep(2.5, length(coef(qufit))), log = TRUE))
        } else {
          p.prior.ceiling <- 0
        }
        
        if(floor.fit=="bayesglm" & floor==TRUE) {
          p.prior.floor <- sum(dcauchy(x = coef(qlfit),
                                       scale = rep(2.5, length(coef(qlfit))), log = TRUE))
        } else {
          p.prior.floor <- 0
        }
        
        llik <- llik - p.prior.floor - p.prior.ceiling
        
        vcov.mle <- solve(-MLEfit$hessian, tol = 1e-20)
        se.mle <- sqrt(diag(vcov.mle))
        
      } # end of boundary ml
      
    } # end of em algorithm
    
  } else if (design=="modified") {

    ## end of standard design
    
    ##
    ## Run two-step estimator for modified design
    
    ## fit to the control group
    
    fit.control <- list()
    par.control.list <- list()
    pi.pred <- matrix(NA, nrow(x.treatment), J)
    
    for (j in 1:J) {
      
      fit.glm.control <- glm(y.control[,j] ~ x.control - 1, family = binomial(logit), weights = w.control)
      coef.glm.control <- coef(fit.glm.control)
      names(coef.glm.control) <- paste("beta", 1:length(coef.glm.control), sep = "")
      
      if (fit.nonsensitive == "glm") {
        fit.control[[j]] <- fit.glm.control
      } else if (fit.nonsensitive == "nls") {

        tmpY <- y.control[, j]
        
        fit.control[[j]] <- nls( as.formula(paste("tmpY ~ logistic( x.control %*% c(",
                                                  paste(paste("beta", 1:length(coef.glm.control),
                                                              sep=""), collapse= ","), "))")) ,
                                start = coef.glm.control,
                                weights = w.control,
                                control = nls.control(maxiter=maxIter, warnOnly=TRUE))
      }
      
      par.control.list[[j]] <- coef(fit.control[[j]])
      
      pi.pred[,j] <- logistic(x.treatment %*% par.control.list[[j]])
      
    }
    
    y.treatment.pred <- y.treatment[,J+1] - apply(pi.pred, 1, sum)
    
    y.treatment.pred.temp <- ifelse(y.treatment.pred > 1, 1, y.treatment.pred)
    y.treatment.pred.temp <- ifelse(y.treatment.pred.temp < 0, 0, y.treatment.pred.temp)
    
    alpha <- mean(y.treatment.pred.temp)
    
    y.treatment.start <- ifelse(y.treatment.pred.temp >=
                                quantile(y.treatment.pred.temp, alpha), 1, 0)
    
    try(fit.glm.control <- glm(cbind(y.treatment.start, 1 - y.treatment.start) ~ x.treatment - 1, family = binomial(logit),
                               weights = w.treatment))
    try(coef.glm.control <- coef(fit.glm.control), silent = T)
    try(names(coef.glm.control) <- paste("beta", 1:length(coef.glm.control), sep = ""), silent = T)
    
    if(exists("coef.glm.control")) {
      try(fit.treat <- nls( as.formula(paste("y.treatment.pred ~ logistic( x.treatment %*% c(",
                                             paste(paste("beta", 1:length(coef.glm.control),sep=""),
                                                   collapse= ","), "))")) , start = coef.glm.control,
                           weights = w.treatment,
                           control = nls.control(maxiter=maxIter, warnOnly=TRUE)), silent = T)
    } else {
      try(fit.treat <- nls( as.formula(paste("y.treatment.pred ~ logistic( x.treatment %*% c(",
                                             paste(paste("beta", 1:length(coef.glm.control),sep=""),
                                                   collapse= ","), "))")) ,
                           weights = w.treatment,
                           control = nls.control(maxiter=maxIter, warnOnly=TRUE)), silent = T)
    }

    if(!exists("fit.treat"))
      fit.treat <- lm(y.treatment.pred ~ x.treatment - 1, weights = w.treatment)
      
    vcov.twostep.modified <- function(coef.treat, coef.control, J,
                                      x1, y1, x0, y0, fit.nonsensitive = "nls") {
      
      nPar <- length(coef.treat)
      
      n <- length(c(y1,y0))
      
      Htmp <- c(logistic(x1 %*% coef.treat) / (1 + exp(x1 %*% coef.treat))) * x1
      
      H <- - t(Htmp) %*% Htmp / n
      
      L <- matrix(NA, ncol=(J*nPar), nrow=nPar)
      
      for(j in 1:J){
        
        Ltmp <- c(sqrt(logistic(x1 %*% coef.control[[j]])*
                       logistic(x1 %*% coef.treat)/((1 + exp(x1 %*% coef.control[[j]]))*
                                                    (1 + exp(x1 %*% coef.treat))))) * x1
        
        L[,((j-1)*nPar+1):((j)*nPar)] <- -t(Ltmp) %*% Ltmp / n
        
      }
      
      M <- matrix(0, ncol = (J*nPar), nrow = (J*nPar))
      
      if (fit.nonsensitive=="glm") {
        
        for(j in 1:J){
        			
          Mtmp <- sqrt(c(exp(x0 %*% coef.control[[j]])/(1+exp(x0 %*% coef.control[[j]])) -
                         exp(2*x0 %*% coef.control[[j]])/(1+exp(x0 %*% coef.control[[j]]))^2))*
                           x0
          
          M[((j-1)*nPar+1):((j)*nPar),((j-1)*nPar+1):((j)*nPar)] <- -t(Mtmp) %*% Mtmp / n
          
        }
        
      } else if (fit.nonsensitive == "nls") {
        
        for(j in 1:J){
          
          Mtmp <- c(logistic(x0 %*% coef.control[[j]])/(1 + exp(x0 %*% coef.control[[j]]))) * x0
          
          M[((j-1)*nPar+1):((j)*nPar),((j-1)*nPar+1):((j)*nPar)] <- -t(Mtmp) %*% Mtmp / n
          
        }
        
      }
      
      invH <- solve(H)
      invM <- solve(M)
      
      invG <- rbind( cbind(invH, -invH %*% L %*% invM),
                    cbind(matrix(0, nrow(invM), ncol(invH)), invM)  ) 
      
      gtmp <- matrix(NA, nrow=nrow(x1), ncol=J)
      for(j in 1:J)
        gtmp[,j] <- logistic(x1 %*% coef.control[[j]])
      
      gtmpsum <- as.vector(apply(gtmp, 1, sum))
      
      g <- c((y1[,J+1] - gtmpsum - logistic(x1 %*% coef.treat))*
             logistic(x1 %*% coef.treat)/(1+exp(x1 %*% coef.treat))) * x1
      
      gg <- t(g) %*% g / n
      
      h <- list()
      
      if (fit.nonsensitive == "glm"){
        for (j in 1:J) {
          h[[j]] <- c((y0[,j] - logistic(x0 %*% coef.control[[j]]))) * x0
        }
      } else if (fit.nonsensitive == "nls") {
        for (j in 1:J) {
          h[[j]] <- c((y0[,j] - logistic(x0 %*% coef.control[[j]]))*
                      logistic(x0 %*% coef.control[[j]])/(1+exp(x0 %*% coef.control[[j]]))) * x0
        }	
      }
      
      Fh <- matrix(NA, nrow = J*nPar, ncol = J*nPar)
      for(j in 1:J) {
        for(k in 1:J){
          Fh[((j-1)*nPar+1):(j*nPar), ((k-1)*nPar+1):(k*nPar)] <- t(h[[j]]) %*% h[[k]] / n
        }
      }
      
      F <- adiag(gg, Fh)
    
      return(invG %*% F %*% t(invG) / n)
      
    }
    
    par.treat.nls.mod <- coef(fit.treat)
    par.control.nls.mod <- do.call(c, par.control.list)
    
    if(method == "nls") {
      vcov.twostep <- vcov.twostep.modified(par.treat.nls.mod, par.control.list,
                                            J, x.treatment, y.treatment, x.control, y.control,
                                            fit.nonsensitive)
      se.twostep <- sqrt(diag(vcov.twostep))
    }
    
    ##
    ## Run maximum likelihood method
    
    if(method == "ml"){

      coef.control.start <- par.control.nls.mod
      coef.treat.start <- par.treat.nls.mod
      
      R.poisbinomR <- function(k, p) {								
        colSum <- 0
        for(i in 1:k){
          if((k-i)>=1) colSum <- colSum + (-1)^(i+1) * sum( (p/(1-p))^i ) * R.poisbinomR(k - i, p)
          if((k-i)==0) colSum <- colSum + (-1)^(i+1) * sum( (p/(1-p))^i )
          if((k-i)<0) colSum <- colSum + 0						
        }
        if(k>0) return((1/k) * colSum)	
        else if(k==0) return(1)
      }
      
      ## Poisson-Binomial density function
      ## takes y scalar and p vector of Bernoulli success probabilities
      dpoisbinomR <- function(y, p) {
        
        return( R.poisbinomR(y, p) * prod(1 - p) )
        
      }
      
      dpoisbinomC <- function(y, p) {
        
        k <- length(p)
        return(.C("dpoisbinom", as.integer(y), as.integer(k), as.double(p),
                  res = double(1), PACKAGE = "list")$res)
        
      }
      
      R.poisbinomC <- function(k, p) {
        
        return(.C("RpoisbinomReturn", as.integer(k), as.double(p),
                  as.integer(length(p)), res = double(1), PACKAGE = "list")$res)
        
      }
      
      R.poisbinomE <- function(k, p) {
        
        maxk <- max(k)
        return(.C("RpoisbinomEff", as.integer(maxk), 
                  as.double(p), as.integer(length(p)), Rs = double(maxk+1),
                  PACKAGE = "list")$Rs[k+1])
        
      }
      
      dpoisbinomE <- function(y, p) {
        
        return(R.poisbinomE(y, p)*prod(1-p))
        
      }
      
      R.poisbinomM <- function(k, p) {
        
        m <- ncol(p)
        if (m != length(k))
          stop("the dimension of k match with the column dimension of p")
        return(.C("RpoisbinomEffMatrix", as.integer(k), as.integer(max(k)),
                  as.double(p), as.integer(nrow(p)), as.integer(m),
                  Rs = double(m), PACKAGE = "list")$Rs)
        
      }
      
      dpoisbinomM <- function(y, p) {
        
        return(R.poisbinomM(y, p)*apply(1-p, 2, prod))
        
      }
      
      obs.llik.modified <- function(par, y, X, treat, wt) {
        
        J <- ncol(y) - 1
        
        nPar <- length(par) / (J+1)
        
        n.treat <- nrow(y[treat==1,])
        
        x.treat <- X[treat==1,]
        y.treat <- y[treat==1,]
        
        pi <- matrix(NA, ncol = nrow(X), nrow = J+1)
        for(j in 1:(J+1))
          pi[j,] <- logistic(X %*% par[((j-1)*nPar+1):(j*nPar)])
        
        llik.treat <- sum(wt[t==1] * log(dpoisbinomM(y.treat[,J+1], pi[,t==1])))
        
        x.control <- X[treat==0,]
        y.control <- y[treat==0,]
        
        llik.control <- 0
        for(j in 1:J)
          llik.control <- llik.control + sum(wt[t==0] * (y.control[,j]*log(pi[j,t==0]) +
                                             (1-y.control[,j])*log(1-pi[j,t==0])))
        
        llik <- llik.treat + llik.control
        
        return(llik)
        
      }
      
      Estep.modified <- function(par, y, X, treat) {
        
        J <- ncol(y) - 1
        
        nPar <- length(par) / (J+1)
        
        x.treat <- X[treat==1, , drop = FALSE]
        y.treat <- y[treat==1, , drop = FALSE]
        
        n.treat <- nrow(x.treat)
        
        pi <- matrix(NA, ncol = n.treat, nrow = J+1)
        for(j in 1:(J+1))
          pi[j,] <- logistic(x.treat %*% par[((j-1)*nPar+1):(j*nPar)])
        
        all0.cond <- y.treat[,J+1]==0
        all1.cond <- y.treat[,J+1]==(J+1)
        either.cond <- all0.cond==TRUE | all1.cond==TRUE
        
        rpb <- matrix(NA, nrow = n.treat, ncol = J+1)
        rpb.inv <- matrix(NA, nrow = n.treat, ncol = J+1)
        
        rpb.inv.vector <- R.poisbinomM(y.treat[!either.cond, J+1], pi[,!either.cond])
        
        for(j in 1:(J+1)){
          rpb[!either.cond ,j] <- R.poisbinomM(y.treat[!either.cond,J+1]-1, pi[-j,!either.cond])
          rpb.inv[!either.cond,j] <- rpb.inv.vector
        }
        
        w <- t(pi) * rpb / ((1-t(pi)) * rpb.inv)
        w[all0.cond, ] <- matrix(0, ncol = ncol(w), nrow = sum(all0.cond))
        w[all1.cond, ] <- matrix(1, ncol = ncol(w), nrow = sum(all1.cond))
        
        return(w)
        
      }
      
      par <- c(coef.control.start, coef.treat.start)
      
      nPar <- length(par) / (J+1)
      
      pllik <- -Inf
      
      llik <- obs.llik.modified(par, y.all, x.all, t, wt = w.all)
      
      counter <- 0
      while (((llik - pllik) > 10^(-8)) & (counter < maxIter)) {
        
        w <- Estep.modified(par, y.all, x.all, t)
        
        y.var <- as.character(formula)[[2]]
                
        coef.list <- list()
        
        for(j in 1:(J+1)){
          
          dtmp <- as.data.frame(rbind(cbind(1, x.treatment, w[,j], 1),
                                      cbind(y.control[,j], x.control, 1, 0),
                                      cbind(0, x.treatment, 1-w[,j], 1)))

          dtmp.weights <- c(w.treatment, w.control, w.treatment)

          if (intercept.only == TRUE)
            names(dtmp) <- c("y","(Intercept)","w","treat")
          else
            names(dtmp) <- c("y","(Intercept)",x.vars,"w","treat")
          
          if(j==(J+1)) {
              dtmp <- dtmp[dtmp$treat==1,]
              dtmp.weights <- dtmp.weights[dtmp$treat==1]
          }
          
          if (j<(J+1)) {
            trials <- 1
          } else {
            trials <- j
          }
          
          fit <- glm(as.formula(paste("cbind(y, 1-y) ~ ", paste(x.vars, collapse=" + "))),
                     family = binomial(logit), weights = dtmp$w * dtmp.weights,
                     start = par[((j-1)*nPar+1):(j*nPar)], data = dtmp)
          
          coef.list[[j]] <- coef(fit)
          
        }
        
        par <- do.call(c, coef.list)
        
        pllik <- llik
        
        if(verbose==T)
          cat(paste(counter, round(llik, 4), "\n"))
        
        llik <- obs.llik.modified(par, y.all, x.all, t)
        
        counter <- counter + 1
        
        if (llik < pllik)
          warning("log-likelihood is not monotonically increasing.")
        ## NOTE CHANGED FROM STOP()
        
        if(counter == (maxIter-1))
          warning("number of iterations exceeded maximum in ML.")
        
      } # end ml while loop
      
      MLEfit <- optim(par, fn = obs.llik.modified, method = "BFGS", y = y.all, X = x.all,
                      treat = t, hessian = TRUE, control = list(maxit = 0))
      
      vcov.mle <- solve(-MLEfit$hessian)
      se.mle <- sqrt(diag(vcov.mle))
      
    } # end ml for modified design
    
  } # end modified design
  
  ## 
  ## Set up return object
  
  if (method == "nls") {

    if (multi == FALSE) {
      
      if(design=="standard")
        par.treat <- par.treat.nls.std
      if(design=="modified")
        par.treat <- par.treat.nls.mod
      se.treat <- se.twostep[1:(length(par.treat))]
      
      if(design=="standard")
        par.control <- par.control.nls.std
      if(design=="modified")
        par.control <- par.control.nls.mod
      se.control <- se.twostep[(length(par.treat)+1):(length(se.twostep))]
      
      names(par.treat) <- names(se.treat) <- coef.names
      
      if (design=="standard")
        names(par.control) <- names(se.control) <- coef.names
      if (design=="modified")
        names(par.control) <- names(se.control) <- rep(coef.names, J)
      
      sum.fit.treat <- summary(fit.treat)
      
      resid.se <- sum.fit.treat$sigma
      resid.df <- sum.fit.treat$df[2]
      
      if(design=="standard") {
        return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.nls, resid.se=resid.se, resid.df=resid.df, coef.names=coef.names,  J=J, design = design, method = method, fit.start = fit.start, overdispersed=overdispersed, boundary = boundary, multi = multi, data = data, x = x.all, y = y.all, treat = t, call = match.call())
        if(weighted == TRUE)
            return.object$weights <- w.all

      } else if (design=="modified") {
        return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.twostep, resid.se=resid.se, resid.df=resid.df, coef.names=coef.names, J=J, design = design, method = method, fit.nonsensitive = fit.nonsensitive, data = data, x = x.all, y = y.all, treat = t, boundary = FALSE, multi = FALSE, call = match.call())
        if(weighted == TRUE)
            return.object$weights <- w.all
        
      }
      
    } else if (multi == TRUE) {

      par.treat <- par.treat.nls.std
      
      se.treat <- list()
      for (m in 1:length(treatment.values))
        se.treat[[m]] <- se.twostep[[m]][1:(length(par.treat[[m]]))]
      
      par.control <- par.control.nls.std
      se.control <- se.twostep[[1]][(length(par.treat[[1]])+1):(length(se.twostep[[1]]))]
      
      for (m in 1:length(treatment.values)) {
        names(par.treat[[m]]) <- coef.names
        names(se.treat[[m]]) <- coef.names
      }
      
      names(par.control) <- names(se.control) <- coef.names

      
      resid.se <- resid.df <- rep(NA, length(treatment.values))
      for (m in 1:length(treatment.values)) {
        sum.fit.treat <- summary(fit.treat[[m]])
        resid.se[m] <- sum.fit.treat$sigma
        resid.df[m] <- sum.fit.treat$df[2]
      }
      
      return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.nls, treat.values = treatment.values, treat.labels = treatment.labels, control.label = control.label, resid.se=resid.se, resid.df=resid.df, J=J,  coef.names=coef.names, design = design, method = method, overdispersed=overdispersed, boundary = boundary, multi = multi, data = data, x = x.all, y = y.all, treat = t, call = match.call())
      if(weighted == TRUE)
          return.object$weights <- w.all

    }
  }
  
  if (method == "ml"){
    
    if(design == "standard") {
      
      if (multi == FALSE) {
        
        if (constrained == T) {

          par.control <- MLEfit$par[1:(nPar)]

          if (overdispersed == T){
            par.treat <- MLEfit$par[(nPar+2):(nPar*2+1)]
            se.treat <- se.mle[(nPar+2):(nPar*2+1)]
            se.control <- se.mle[1:(nPar)]
            par.overdispersion <- MLEfit$par[nPar+1]
            se.overdispersion <- se.mle[nPar+1]
            names(par.overdispersion) <- names(se.overdispersion) <- "overdispersion"
          } else {
            par.treat <- MLEfit$par[(nPar+1):(nPar*2)]
            se.treat <- se.mle[(nPar+1):(nPar*2)]
            se.control <- se.mle[1:(nPar)]
          }
          
          if(floor==TRUE)
            par.floor <- par[(nPar*2+1):(nPar*2 + nPar.floor)]
          if(floor==TRUE & ceiling==TRUE)
            par.ceiling <- par[(nPar*2 + nPar.floor + 1):(nPar*2 + nPar.floor + nPar.ceiling)]
          if(floor==FALSE & ceiling==TRUE)
            par.ceiling <- par[(nPar*2+1):(nPar*2 + nPar.ceiling)]
          
          if(floor==TRUE)
            se.floor <- se.mle[(nPar*2+1):(nPar*2+ nPar.floor)]
          if(floor==TRUE & ceiling==TRUE)
            se.ceiling <- se.mle[(nPar*2 + nPar.floor + 1):(nPar*2 + nPar.floor + nPar.ceiling)]
          if(floor==FALSE & ceiling==TRUE)
            se.ceiling <- se.mle[(nPar*2+1):(nPar*2 + nPar.ceiling)]
          
          names(par.treat) <- names(se.treat) <- names(par.control) <- names(se.control)  <- coef.names
          
          if(floor==TRUE)
            names(par.floor) <- names(se.floor) <- coef.names.floor
          if(ceiling==TRUE)
            names(par.ceiling) <- names(se.ceiling) <- coef.names.floor
          
          if(boundary == TRUE | multi == TRUE)
            llik.const <- llik
          
          if(boundary==F)
            if (overdispersed == T)
              return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, par.overdispersion=par.overdispersion, se.overdispersion=se.overdispersion, vcov=vcov.mle, pred.post = w, treat.labels = treatment.labels, control.label = control.label, llik=llik.const, J=J,  coef.names=coef.names, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
            else
              return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.mle, pred.post = w, treat.labels = treatment.labels, control.label = control.label, llik=llik.const, J=J,  coef.names=coef.names, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
          
          if(floor==FALSE & ceiling==TRUE)
            return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, par.ceiling = par.ceiling, se.ceiling = se.ceiling, vcov=vcov.mle, pred.post = w, treat.labels = treatment.labels, control.label = control.label, llik=llik.const, J=J,  coef.names=coef.names, coef.names.ceiling = coef.names.ceiling, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
          
          if(floor==TRUE & ceiling==FALSE)
            return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, par.floor = par.floor, se.floor = se.floor, vcov=vcov.mle, pred.post = w, llik=llik.const, treat.labels = treatment.labels, control.label = control.label, J=J, coef.names=coef.names, coef.names.floor = coef.names.floor, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
          
          if(floor==TRUE & ceiling==TRUE)
            return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, par.floor = par.floor, se.floor = se.floor, par.ceiling = par.ceiling, se.ceiling = se.ceiling, pred.post = w, vcov=vcov.mle, treat.labels = treatment.labels, control.label = control.label, llik=llik.const, J=J,  coef.names=coef.names, coef.names.floor = coef.names.floor, coef.names.ceiling = coef.names.ceiling, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, ceiling = ceiling, floor = floor, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
          
        } else if (constrained == FALSE) { 
          
          par.treat <- MLEfit$par[(nPar*2+1):(nPar*3)]
          par.control.psi0 <- MLEfit$par[1:(nPar)]
          par.control.psi1 <- MLEfit$par[(nPar+1):(nPar*2)]
          
          if (overdispersed == T){
            se.treat <- se.mle[(nPar*2+3):(nPar*3 + 2)]
            se.control.psi0 <- se.mle[1:(nPar)]
            se.control.psi1 <- se.mle[(nPar*2+2):(nPar*2+1)]
            par.overdispersion <- MLEfit$par[nPar*2+2]
            se.overdispersion <- se.mle[nPar*2+2]
            names(par.overdispersion) <- names(se.overdispersion) <- "overdispersion"
          } else {
            se.treat <- se.mle[(nPar*2+1):(nPar*3)]
            se.control.psi0 <- se.mle[1:(nPar)]
            se.control.psi1 <- se.mle[(nPar+1):(nPar*2)]
          }
          
          names(par.treat) <- names(se.treat) <- names(par.control.psi0) <- names(se.control.psi0) <- names(par.control.psi1) <- names(se.control.psi1) <- coef.names

          if(overdispersed==T)
            return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control.psi0=par.control.psi0, se.control.psi0=se.control.psi0, par.control.psi1=par.control.psi1, se.control.psi1=se.control.psi1, par.overdispersion=par.overdispersion, se.overdispersion=se.overdispersion, vcov=vcov.mle, pred.post = w, treat.labels = treatment.labels, control.label = control.label, llik=llik, J=J,  coef.names=coef.names, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, call = match.call(), data=data, x = x.all, y = y.all, treat = t)
          else
            return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control.psi0=par.control.psi0, se.control.psi0=se.control.psi0, par.control.psi1=par.control.psi1, se.control.psi1=se.control.psi1, vcov=vcov.mle, pred.post = w, treat.labels = treatment.labels, control.label = control.label, llik=llik, J=J,  coef.names=coef.names, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, call = match.call(), data=data, x = x.all, y = y.all, treat = t)
	
        }

        if(weighted == TRUE)
            return.object$weights <- w.all

      } else if (multi == TRUE) {

        par.control <- MLEfit$par[1:(nPar)]
        se.control <- se.mle[1:(nPar)]

        if (multi.condition == "none") {
          se.treat <- list()
          for (m in 1:length(treatment.values))
            se.treat[[m]] <- se.mle[(nPar + (m-1) * nPar + 1) : (nPar + m * nPar)]
          
          for (m in 1:length(treatment.values)) {
            names(par.treat[[m]]) <- coef.names
            names(se.treat[[m]]) <- coef.names
          }
        } else if (multi.condition == "level") {
          se.treat <- list()
          for (m in 1:length(treatment.values))
            se.treat[[m]] <- se.mle[(nPar + (m-1) * (nPar+1) + 1) :
                                    (nPar + m * (nPar+1))]
          
          for (m in 1:length(treatment.values)) {
            names(par.treat[[m]]) <- c(coef.names, "y_i(0)")
            names(se.treat[[m]]) <- c(coef.names, "y_i(0)")
          }
        }
        
        names(par.control) <- names(se.control) <- coef.names
        
        return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.mle, pred.post = w, treat.values = treatment.values, treat.labels = treatment.labels, control.label = control.label, multi.condition = multi.condition, llik=llik, J=J,  coef.names=coef.names, design = design, method = method, overdispersed=overdispersed, constrained=constrained, boundary = boundary, multi = multi, call = match.call(), data = data, x = x.all, y = y.all, treat = t)
        if(weighted == TRUE)
            return.object$weights <- w.all

      }
      
    } else if (design == "modified") {
      
      par.treat <- MLEfit$par[(nPar*J+1):(nPar*(J+1))]
      par.control <- MLEfit$par[1:(nPar*J)]
      
      se.treat <- se.mle[(nPar*J+1):(nPar*(J+1))]
      se.control <- se.mle[1:(nPar*J)]
      
      names(par.treat) <- names(se.treat) <- coef.names
      
      names(par.control) <- names(se.control) <- rep(coef.names, J)
      
      return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.mle, llik=llik, treat.labels = treatment.labels, control.label = control.label, coef.names=coef.names,  J=J, design = design, method = method, boundary = FALSE, multi = FALSE, call = match.call(), data=data, x = x.all, y = y.all, treat = t)
      if(weighted == TRUE)
          return.object$weights <- w.all
    }
  }
  
  
  class(return.object) <- "ictreg"
  
  return.object
  
}


print.ictreg <- function(x, ...){
  
  cat("\nItem Count Technique Regression \n\nCall: ")
  
  dput(x$call)

  cat("\nCoefficient estimates\n")

  tb <- as.matrix(x$par.treat)
  colnames(tb) <- "est."
  
  cat("\n")

  print(coef(x))

  cat("\n")

  treat.print <- c()
  for (i in 1:length(x$treat.labels)) {
    treat.print <- c(treat.print, "'", x$treat.labels[i], "'", sep = "")
    if (i != length(x$treat.labels))
      treat.print <- c(treat.print, " and ")
  }
  
  cat("Number of control items J set to ", x$J, ". Treatment groups were indicated by ", sep = "")
  cat(treat.print, sep ="")
  cat(" and the control group by '", x$control.label, "'.\n\n", sep = "")
     
  invisible(x)
  
}

print.predict.ictreg <- function(x, ...){
  
  cat("\nList Experiment Prediction\n")
  
  cat("\nProportion of affirmative responses to the sensitive item\n")

  if (class(x$fit) == "data.frame")
    n <- nrow(x$fit)
  else
    n <- length(x$fit)

  for (i in 1:n) {
    cat(paste("\nPrediction:", rownames(as.matrix(x$fit))[i], "\nEst. = ",
              round(as.matrix(x$fit)[i,1], 4)))
    if (is.null(x$se.fit) == FALSE) {
        cat(paste(" (s.e. = ", round(as.matrix(x$se.fit)[i], 4), ")", sep = ""))
    }
    if (class(x$fit) == "data.frame") {
      cat(paste("\n95% confidence interval: (", round(as.matrix(x$fit)[i,2], 4), ", ",
                round(as.matrix(x$fit)[i,3], 4), ")", sep = ""))
    }
    cat("\n")
  }
  
  cat("\n")

  invisible(x)
  
}

predict.ictreg <- function(object, newdata, newdata.diff, direct.glm, se.fit = FALSE,
                           interval = c("none","confidence"), level = .95, avg = FALSE, sensitive.item, ...){

  ##if(object$design != "standard")
  ##  stop("The predict function only currently works for standard design, single sensitive item experiments without ceiling or floor effects.")
  
  if(missing(interval)) interval <- "none"
  
  nPar <- length(object$coef.names)
  
  ## dummy value for constraint
  if(object$method!="ml") object$constrained <- F
  
  logistic <- function(object) exp(object)/(1+exp(object))

  if(missing(sensitive.item)) {
    sensitive.item <- 1
    if(object$multi==TRUE)
      warning("Using the first sensitive item for predictions. Change with the sensitive.item option.")
  }

  if (class(object$par.treat) != "list") {
    ## extract only treatment coef's, var/covar
    beta <- object$par.treat ## was coef(object)[1:nPar]
    var.beta <- vcov(object)[1:nPar, 1:nPar]
  } else {
    ## extract only treatment coef's, var/covar
    beta <- object$par.treat[[sensitive.item]] ## was coef(object)[1:nPar]
    var.beta <- vcov(object)[(sensitive.item*nPar + 1):((sensitive.item+1)*nPar), (sensitive.item*nPar + 1):((sensitive.item+1)*nPar)]
  }
  
  
  if(missing(newdata)) {
    xvar <- object$x
  } else { 
    if(nrow(newdata)==0)
      stop("No data in the provided data frame.")
    ##formula <- do.call(c, as.list(as.character(print(object$call$formula[[3]]))))
    xvar <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))), newdata)
  }

  if (missing(direct.glm) == TRUE) {
    
    if (!missing(newdata.diff)) {
      data.list <- list(xvar,
                        model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))),
                                     newdata.diff))
    } else {
      data.list <- list(xvar)
    }
    
    k <- 1
    
    multi.return.object <- list()
    
    for (xvar in data.list) {
      
      n <- nrow(xvar)
      
      if (object$method == "lm")
        pix <- xvar %*% beta
      else
        pix <- logistic(xvar %*% beta)
      
      v <- pix/(1+exp(xvar %*% beta))
      
      if(avg==FALSE){
        var.pred <- rep(NA, n)
        for (i in 1:n) {
          if (object$method == "lm")
            var.pred[i] <- xvar[i,, drop = FALSE] %*% var.beta %*% t(xvar[i,, drop = FALSE])
          else
            var.pred[i] <- v[i] * v[i] * (xvar[i,, drop = FALSE] %*% var.beta %*%
                                          t(xvar[i,, drop = FALSE]))
        }
        
        se <- sqrt(var.pred)
        
      } else {
        
        pix <- mean(pix)
        
        var.pred <- 0
        for (i in 1:n) {
          for (j in 1:n) {
            if (object$method == "lm")
              var.pred <- var.pred + xvar[i,, drop = FALSE] %*% var.beta %*%
                t(xvar[j,, drop = FALSE])
            else
              var.pred <- var.pred + v[i] * v[j] * (xvar[i,, drop = FALSE] %*%
                                                    var.beta %*% t(xvar[j,, drop = FALSE]))
          }
        }
        
        se <- c(sqrt(var.pred)/n)
        
      }
      
      return.object <- list(fit=as.vector(pix))
      
      if(interval=="confidence"){
        
        critical.value <- qt(1-(1-level)/2, df = n)
        ci.upper <- pix + critical.value * se
        ci.lower <- pix - critical.value * se
        
        return.object <- list(fit=data.frame(fit = pix, lwr = ci.lower, upr = ci.upper))
        
      }
      
      if(se.fit==T)
        return.object$se.fit <- as.vector(se)
      
      if (length(data.list) == 2)
        multi.return.object[[k]] <- return.object
      
      k <- k + 1
      
    }
    
    if (length(data.list) == 2) {
      
      if (object$method == "lm")
        stop("Currently the difference option does not work with the lm method. Try ml.")
      
      xa <- data.list[[1]]
      xb <- data.list[[2]]
      
      n <- nrow(xa)
      var <- est <- 0
      pixa <- logistic(xa %*% beta)
      pixb <- logistic(xb %*% beta)
      va <- pixa/(1+exp(xa %*% beta))
      vb <- pixb/(1+exp(xb %*% beta))
      
      for (i in 1:n) {
        for (j in 1:n) {
          var <- var + va[i] * va[j] * (xa[i,, drop = FALSE] %*% var.beta %*% t(xa[j,, drop = FALSE]))
          - 2 * va[i] * vb[j] * (xa[i,, drop = FALSE] %*% var.beta %*% t(xb[j,, drop = FALSE]))
          + vb[i] * vb[j] * (xb[i,, drop = FALSE] %*% var.beta %*% t(xb[j,, drop = FALSE]))
        }
      }
      
      est <- mean(pixa-pixb)
      
      se <- as.vector(sqrt(var)) / n
      
      return.object <- list(fit = est)
      
      if(interval=="confidence"){
        
        critical.value <- qt(1-(1-level)/2, df = n)
        ci.upper <- est + critical.value * se
        ci.lower <- est - critical.value * se
        
        return.object <- list(fit=data.frame(fit = est, lwr = ci.lower, upr = ci.upper))
        
      }
      
      if(se.fit==T)
        return.object$se.fit <- as.vector(se)
      
      multi.return.object[[3]] <- return.object
      
      ## now put them together
      
      fit <- se <- list()
      for (j in 1:3) {
        fit[[j]] <- multi.return.object[[j]]$fit
        se[[j]] <- multi.return.object[[j]]$se.fit
      }
      
      return.object <- c()
      return.object$fit <- as.data.frame(do.call(rbind, fit))

      rownames(return.object$fit) <- c("newdata", "newdata.diff", "Difference (newdata - newdata.diff)")
      if (se.fit == T)
        return.object$se.fit <- as.vector(do.call(c, se))

      attr(return.object, "concat") <- TRUE
      
    }


  } else {

    ## direct versus list Monte Carlo prediction

    beta.direct <- coef(direct.glm)
    vcov.direct <- vcov(direct.glm)

    n.draws <- 10000
        
    xvar.direct <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))),
                                direct.glm$data)
    
    if (ncol(xvar.direct) != nPar)
      stop("Different number of covariates in direct and list regressions.")
     
    draws.list <- mvrnorm(n = n.draws, mu = beta, Sigma = var.beta)
    draws.direct <- mvrnorm(n = n.draws, mu = beta.direct, Sigma = vcov.direct)

    pred.list.mean <- pred.direct.mean <- pred.diff.mean <- rep(NA, n.draws)
    
    for (d in 1:n.draws) {
      
      par.g <- draws.list[d, ]

      if (object$method == "lm")   
        pred.list <- xvar %*% par.g
      else
        pred.list <- logistic(xvar %*% par.g)

      pred.direct <- logistic(xvar.direct %*% draws.direct[d,])
      
      pred.list.mean[d] <- mean(pred.list)
      pred.direct.mean[d] <- mean(pred.direct)
      pred.diff.mean[d] <- pred.list.mean[d] - pred.direct.mean[d]
      
    }

    est.list <- mean(pred.list.mean)
    se.list <- sd(pred.list.mean)
    est.direct <- mean(pred.direct.mean)
    se.direct <- sd(pred.direct.mean)
    est.diff <- mean(pred.diff.mean)
    se.diff <- sd(pred.diff.mean)
    
    critical.value <- qt(1-(1-level)/2, df = nrow(xvar))
    ci.upper.list <- est.list + critical.value * se.list
    ci.lower.list <- est.list - critical.value * se.list
    ci.upper.direct <- est.direct + critical.value * se.direct
    ci.lower.direct <- est.direct - critical.value * se.direct
    ci.upper.diff <- est.diff + critical.value * se.diff
    ci.lower.diff <- est.diff - critical.value * se.diff
       
    fit.matrix <- as.data.frame(rbind(c(est.list, ci.lower.list, ci.upper.list),
                                         c(est.direct, ci.lower.direct, ci.upper.direct),
                                         c(est.diff, ci.lower.diff, ci.upper.diff)))
    names(fit.matrix) <- c("fit","lwr","upr")
    rownames(fit.matrix) <- c("List", "Direct", "Difference (list - direct)")

    return.object <- c()
    return.object$fit <- fit.matrix
    if (se.fit == T)
      return.object$se.fit <- c(se.list, se.direct, se.diff)

    attr(return.object, "concat") <- TRUE
    
  }
  
  class(return.object) <- "predict.ictreg"
  
  return.object
  
}

coef.ictreg <- function(object, ...){
  
  nPar <- length(object$coef.names)
  
  if((object$method=="lm" | object$method=="nls" | object$design=="modified") & object$boundary == FALSE & object$multi == FALSE){
    coef <- c(object$par.treat,object$par.control)
    
    if(object$design=="standard") {
      names(coef) <- c(paste("sensitive.",object$coef.names,sep=""), paste("control.",object$coef.names,sep=""))
    } else if(object$design=="modified") {
      coef.names <- paste("sensitive.",object$coef.names,sep="")
      for(j in 1:object$J)
        coef.names <- c(coef.names, paste("control.", j, ".", object$coef.names, sep=""))
      names(coef) <- coef.names
    }
    
  } else if(object$method=="ml" & object$boundary == FALSE & object$multi == FALSE) {
    
    if(object$design=="standard"){
      if(object$constrained==F) {
        coef <- c(object$par.treat,object$par.control.psi0,object$par.control.psi1)
        names(coef) <- c(paste("sensitive.",object$coef.names,sep=""), paste("control.psi0.",object$coef.names,sep=""),paste("control.psi1.",object$coef.names,sep=""))
      } else if (object$design=="modified") {
        coef <- c(object$par.treat,object$par.control)
        names(coef) <- c(paste("sensitive.",object$coef.names,sep=""), paste("control.",object$coef.names,sep=""))
      } else if (object$constrained == T) {
        coef <- c(object$par.treat,object$par.control)
        names(coef) <- c(paste("sensitive.",object$coef.names,sep=""),
                         paste("control.",object$coef.names,sep=""))
      } 
    } 
  }

  if (object$boundary == TRUE) {
    if (object$floor == TRUE & object$ceiling == TRUE) {
      coef <- c(object$par.control, object$par.treat,
                object$par.floor, object$par.ceiling)
      coef.names <- paste("control.",object$coef.names,sep="")
      coef.names <- c(coef.names, paste("sensitive.",object$coef.names,sep=""))
      coef.names <- c(coef.names, paste("floor.",object$coef.names.floor,sep=""))
      coef.names <- c(coef.names, paste("ceiling.",object$coef.names.ceiling,sep=""))
      names(coef) <- coef.names
    } else if (object$floor == TRUE){
      coef <- c(object$par.control, object$par.treat, object$par.floor)
      coef.names <- paste("control.",object$coef.names,sep="")
      coef.names <- c(coef.names, paste("sensitive.",object$coef.names,sep=""))
      coef.names <- c(coef.names, paste("floor.",object$coef.names.floor,sep=""))
      names(coef) <- coef.names
    } else  if (object$ceiling == TRUE) {
      coef <- c(object$par.control, object$par.treat, object$par.ceiling)
      coef.names <- paste("control.",object$coef.names,sep="")
      coef.names <- c(coef.names, paste("sensitive.",object$coef.names,sep=""))
      coef.names <- c(coef.names, paste("ceiling.",object$coef.names.ceiling,sep=""))
      names(coef) <- coef.names
    }
  }

  if(object$multi == TRUE) {

    if (object$method == "nls" | object$method == "lm")
      object$multi.condition <- "none"
    
    coef <- c(do.call(c, object$par.treat), object$par.control)
    coef.names <- c() 
    if(object$multi.condition == "none")
      for(j in 1:length(object$treat.labels)) coef.names <- c(coef.names, paste("sensitive.", object$treat.labels[j], ".", object$coef.names,sep=""))
    if(object$multi.condition == "level")
      for(j in 1:length(object$treat.labels)) coef.names <- c(coef.names, paste("sensitive.", object$treat.labels[j], ".", c(object$coef.names, "y_i(0)"),sep=""))
    coef.names <- c(coef.names, paste("control.",object$coef.names,sep=""))
    names(coef) <- coef.names

  }
  
  return(coef)
  
}

vcov.ictreg <- function(object, ...){
  
  vcov <- object$vcov
  
  nPar <- length(object$coef.names)
  
  ## dummy value for constraint
  if(object$method=="nls" | object$design=="modified" | object$method == "lm")
    object$constrained <- F
  
  if (object$method == "lm" | (object$method=="ml" & object$constrained==T &
     object$boundary == F & object$multi == F)) {
    vcov <- rbind( cbind( vcov[(nPar+1):(nPar*2), (nPar+1):(nPar*2)],
                         vcov[(nPar+1):(nPar*2), 1:nPar]  ),
                  cbind(vcov[1:nPar, (nPar+1):(nPar*2)] ,vcov[1:nPar, 1:nPar])  )
  } else if (object$method=="ml" & object$constrained==F & object$boundary == F & object$multi == F) {
    vcov <- rbind( cbind(vcov[(nPar*2+1):(nPar*3),(nPar*2+1):(nPar*3)],
                         vcov[(nPar*2+1):(nPar*3),1:nPar],
                         vcov[(nPar*2+1):(nPar*3), (nPar+1):(nPar*2)] ),
                  cbind(vcov[1:nPar, (nPar*2+1):(nPar*3)] , vcov[1:nPar,1:nPar] ,
                        vcov[1:nPar, (nPar+1):(nPar*2)]),
                  cbind(vcov[(nPar+1):(nPar*2), (nPar*2+1):(nPar*3)] ,
                        vcov[(nPar+1):(nPar*2), 1:nPar] ,
                        vcov[(nPar+1):(nPar*2),(nPar+1):(nPar*2)]) )
  }

  if((object$method=="nls" | object$method == "lm") &
     object$boundary == FALSE & object$multi == FALSE){
    if(object$design=="standard") rownames(vcov) <-
      colnames(vcov)<- c(paste("sensitive.",object$coef.names,sep=""),
                         paste("control.",object$coef.names,sep=""))
    else if(object$design=="modified") {
      coef.names<- paste("sensitive.",object$coef.names,sep="")
      for(j in 1:object$J) coef.names <-
        c(coef.names, paste("control.", j, ".", object$coef.names, sep=""))
      rownames(vcov) <- colnames(vcov) <- coef.names
    }
  } else if(object$method=="ml" & object$boundary == FALSE & object$multi == FALSE) {
    if(object$constrained==F) {
      rownames(vcov) <- colnames(vcov) <- c(paste("sensitive.",object$coef.names,sep=""),
                                            paste("control.psi0.",object$coef.names,
                                                  sep=""),
                                            paste("control.psi1.",object$coef.names,
                                                  sep=""))
    } else {
      rownames(vcov) <- colnames(vcov) <-
        c(paste("sensitive.",object$coef.names,sep=""),
          paste("control.",object$coef.names,sep=""))
    }
  }

  if (object$boundary == TRUE) {
    if (object$floor == TRUE & object$ceiling == TRUE) {
      coef.names <- paste("control.",object$coef.names,sep="")
      coef.names <- c(coef.names, paste("sensitive.",object$coef.names,sep=""))
      coef.names <- c(coef.names, paste("floor.",object$coef.names.floor,sep=""))
      coef.names <- c(coef.names, paste("ceiling.",object$coef.names.ceiling,sep=""))
      rownames(vcov) <- colnames(vcov) <- coef.names
    } else if (object$floor == TRUE){
      coef.names <- paste("control.",object$coef.names,sep="")
      coef.names <- c(coef.names, paste("sensitive.",object$coef.names,sep=""))
      coef.names <- c(coef.names, paste("floor.",object$coef.names.floor,sep=""))
      rownames(vcov) <- colnames(vcov) <- coef.names
    } else  if (object$ceiling == TRUE) {
      coef.names <- paste("control.",object$coef.names,sep="")
      coef.names <- c(coef.names, paste("sensitive.",object$coef.names,sep=""))
      coef.names <- c(coef.names, paste("ceiling.",object$coef.names.ceiling,sep=""))
      rownames(vcov) <- colnames(vcov) <- coef.names
    }
  }

  if(object$multi == TRUE) {

    ##if (object$method == "nls")
    ##  stop("This function does not yet work for the NLS estimator for the multiple sensitive item design.")

    coef.names <- paste("control.",object$coef.names,sep="")
    if(object$method != "ml") {
      for(j in 1:length(object$treat.values)) coef.names <-
        c(coef.names, paste("sensitive.", object$treat.labels[j], ".", object$coef.names,sep=""))
    } else {
      if(object$multi.condition == "none")
        for(j in 1:length(object$treat.values)) coef.names <-
          c(coef.names, paste("sensitive.", object$treat.labels[j], ".", object$coef.names,sep=""))
      else if(object$multi.condition == "level")
        for(j in 1:length(object$treat.values)) coef.names <-
          c(coef.names, paste("sensitive.", object$treat.labels[j], ".", c(object$coef.names, "y_i(0)"),sep=""))
      rownames(vcov) <- colnames(vcov) <- coef.names
    }
  }
  
  return(vcov)
  
}

c.predict.ictreg <- function(...){

  x <- list(...)
  
  for (j in 1:length(x)) 
    x[[j]] <- x[[j]]$fit
  
  return.object <- c()
  return.object$fit <- as.matrix(do.call(rbind, x))

  class(return.object) <- "predict.ictreg"

  return.object
  
}

plot.predict.ictreg <- function(x, labels = NA, axes.ict = TRUE,
                                xlim = NULL, ylim = NULL, xlab = NULL, ylab = "Estimated Proportion",
                                axes = F, pch = 19, xvec = NULL, ...){

  x <- as.matrix(x$fit)

  if (is.null(xlim))
    xlim <- c(0.5,nrow(x)+0.5)
  if (is.null(ylim))
    ylim <- c(min(x[,2]-.05,0), max(x[,3]+.05))  
  if (is.null(xvec))
    xvec <- 1:nrow(x)  
  if (is.null(xlab))
    xlab <- ""

  plot(xvec, x[,1], pch = pch, xlim = xlim, ylim = ylim, 
       xlab = xlab, ylab = ylab, axes = axes, ...)
  
  critical <- abs(qnorm(0.025))
  for (j in 1:nrow(x))
    lines(rep(xvec[j], 2), x[j, 2:3]) 
  
  abline(h = 0, lty = "dashed")
  
  if (axes.ict == TRUE) {
    axis(side = 2, at = seq(from = round(min(x[,2],0),1),
                     to = round(max(x[,3]),1), by = 0.1))
    axis(side = 1, at = xvec, tick = FALSE, labels = labels)
  }
  
}

summary.ictreg <- function(object, boundary.proportions = FALSE, n.draws = 10000, ...) {
  object$boundary.proportions <- boundary.proportions
  object$n.draws <- n.draws
  structure(object, class = c("summary.ictreg", class(object)))
}

print.summary.ictreg <- function(x, ...){
  
  cat("\nItem Count Technique Regression \n\nCall: ")
  
  dput(x$call)

  cat("\n")
  
  if(x$method=="nls" | x$design=="modified" | x$method == "lm")
    x$constrained <- T

  if (x$design == "standard") {
  
    if((x$method=="nls" | x$method == "lm") & x$multi == FALSE){
      ##cat(rep(" ", max(nchar(x$coef.names))+8), "Sensitive Item", rep(" ", 5),
      ##    "Control Items \n", sep="")
      tb.treat <- tb.control <- matrix(NA, ncol = 2, nrow = length(x$par.control))
      colnames(tb.treat) <- colnames(tb.control) <- c("Est.", "S.E.")
      if(x$design=="standard") {
        rownames(tb.treat) <- rownames(tb.control) <- x$coef.names
      } else if(x$design=="modified") {
        rownames(tb.treat) <- rownames(tb.control) <- rep(x$coef.names, x$J)
      }
      tb.treat[,1] <- x$par.treat
      tb.treat[,2] <- x$se.treat
      tb.control[,1] <- x$par.control
      tb.control[,2] <- x$se.control
      
      cat("Sensitive item \n")
      print(as.matrix(round(tb.treat,5)))
      cat("\nControl items \n")
      print(as.matrix(round(tb.control,5)))    
      
      summ.stat <- paste("Residual standard error:",signif(x$resid.se, digits=6),
                         "with",x$resid.df,"degrees of freedom")
      
      cat("\n",summ.stat,"\n\n", sep="")
      
    } else if(x$method=="ml" | x$multi == TRUE) {
      
      if (x$multi == FALSE) {
        if(x$constrained==F) {
          ##cat(rep(" ", max(nchar(x$coef.names))+7), "Sensitive Item", rep(" ", 3),
          ##    "Control (psi0)", rep(" ",2), "Control (psi1)\n", sep="")
          tb.treat <- tb.control.psi0 <- tb.control.psi1 <- matrix(NA, ncol = 2, nrow = length(x$par.treat))
          colnames(tb.treat) <- colnames(tb.control.psi0) <- colnames(tb.control.psi1) <- c("Est.", "S.E.")
          rownames(tb.treat) <- rownames(tb.control.psi0) <- rownames(tb.control.psi1) <- x$coef.names
          tb.treat[,1] <- x$par.treat
          tb.treat[,2] <- x$se.treat
          tb.control.psi0[,1] <- x$par.control.psi0
          tb.control.psi0[,2] <- x$se.control.psi0
          tb.control.psi1[,1] <- x$par.control.psi1
          tb.control.psi1[,2] <- x$se.control.psi1

          cat("Sensitive item \n")
          print(as.matrix(round(tb.treat, 5)))
          cat("\nControl items (psi0) \n")
          print(as.matrix(round(tb.control.psi0,5)))
          cat("\nControl items (psi1) \n")
          print(as.matrix(round(tb.control.psi1,5)))
          
        } else {
          #cat(rep(" ", max(nchar(x$coef.names))+8), "Sensitive Item", rep(" ", 5),
          #    "Control Items \n", sep="")
          tb.treat <- tb.control <- matrix(NA, ncol = 2, nrow = length(x$par.treat))
          colnames(tb.treat) <- colnames(tb.control) <- c("Est.", "S.E.")
          rownames(tb.treat) <- rownames(tb.control) <- x$coef.names
          tb.treat[,1] <- x$par.treat
          tb.treat[,2] <- x$se.treat
          tb.control[,1] <- x$par.control
          tb.control[,2] <- x$se.control
          
          cat("Sensitive item \n")
          print(as.matrix(round(tb.treat,5)))
          cat("\nControl items \n")
          print(as.matrix(round(tb.control,5)))
          
        }
        
        if (x$overdispersed == TRUE)
          cat(paste("\nOverdispersion parameter: ",
                    round(x$par.overdispersion, 4), " (s.e. = ", round(x$se.overdispersion, 4), ").\n",
                    sep = ""))
        
      } else {
        ## multi
        
        if (x$method == "nls" | x$method == "lm")
          x$multi.condition <- "none"
        
        ##cat(rep(" ", max(nchar(x$coef.names))+5), x$treat.labels[1], sep="")
        ##for(k in 2:length(x$treat.values))
        ##  cat(rep(" ", 7), x$treat.labels[k], "\n", sep = "")
        
        tb.treat <- tb.control <- matrix(NA, ncol = 2, nrow = length(x$par.treat[[1]]))
        colnames(tb.treat) <- c("Est.", "S.E.")
        if(x$multi.condition == "none")
          rownames(tb.treat) <- x$coef.names
        else if (x$multi.condition == "level")
          rownames(tb.treat) <- c(x$coef.names, "y_i(0)")
        else if (x$multi.condition == "indicators")
          rownames(tb.treat) <- rep("varnames", length(x$par.treat[[1]]))

        for (k in 1:length(x$treat.values)) {

          tb.treat[,1] <- x$par.treat[[k]]
          tb.treat[,2] <- x$se.treat[[k]]

          cat(paste("Sensitive item (", x$treat.labels[k], ")", "\n", sep = ""))
          print(as.matrix(round(tb.treat,5)))
   
          cat("\n")

        }
                
        ##cat("\n",rep(" ", max(nchar(x$coef.names))+16), "Control\n", sep="")
        cat("Control items\n")
        tb.control <- matrix(NA, ncol = 2, nrow = length(x$par.control))
        colnames(tb.control) <- rep(c("Est.", "S.E."), 1)
        rownames(tb.control) <- x$coef.names
        tb.control[,1] <- x$par.control
        tb.control[,2] <- x$se.control
        
        print(as.matrix(round(tb.control, 5)))
        
        if (x$method == "lm")
          summ.stat <- paste("Residual standard error:",signif(x$resid.se, digits=6),
                             "with",x$resid.df,"degrees of freedom")
        
      }
      
      if(x$method == "ml")
        summ.stat <- paste("Log-likelihood:", round(x$llik, 3))
      else if (x$method == "nls") {
        summ.stat <- "Residual standard error of "
        for (m in 1:length(x$treat.values)) {
          if(m==1)
            summ.stat <- paste(summ.stat, signif(x$resid.se[m], digits=6), " with ",x$resid.df[m]," degrees of freedom for sensitive item ",m,"", sep = "")
          else
            summ.stat <- paste(summ.stat, "; residual standard error of ", signif(x$resid.se[m], digits=6), " with ",x$resid.df[m]," degrees of freedom for sensitive item ",m,"", sep ="")
        }
        summ.stat <- paste(summ.stat, ".", sep = "")
      }
      
      if(x$boundary == FALSE)
        cat("\n",summ.stat,"\n\n", sep="")
      
    }
    
    if (x$boundary == TRUE & x$method == "ml" & x$design == "standard") {
      
      tb.ceiling <- tb.floor <- matrix(NA, ncol = 2, nrow = length(x$par.treat))
      colnames(tb.ceiling) <- colnames(tb.floor) <- rep(c("Est.", "S.E."), 1)
      rownames(tb.ceiling) <- rownames(tb.floor) <- x$coef.names
      
      if (x$ceiling == TRUE) {
        cat("\nCeiling\n", sep="")
        tb.ceiling[,1] <- x$par.ceiling
        tb.ceiling[,2] <- x$se.ceiling
        print(as.matrix(round(tb.ceiling, 5)))
      }
      if (x$floor == TRUE) {
        cat("\nFloor\n", sep="")
        tb.floor[,1] <- x$par.floor
        tb.floor[,2] <- x$se.floor
        print(as.matrix(round(tb.floor, 5)))
     }
      
      cat("\n",summ.stat,"\n\n", sep="")
      
      if (x$boundary.proportions == TRUE) {
        
        logistic <- function(x) exp(x)/(1+exp(x))
        logit <- function(x) return(log(x)-log(1-x))
        
        nPar <- ncol(x$x)
        n <- nrow(x$x)
                
        mu.par <- c(x$par.control, x$par.treat)
        if (x$floor == TRUE)
          mu.par <- c(mu.par, x$par.floor)
        if (x$ceiling == TRUE)
          mu.par <- c(mu.par, x$par.ceiling)
        
        try(draws <- mvrnorm(n = x$n.draws, mu = mu.par, Sigma = x$vcov, tol = 1e-1))
        
        if (exists("draws")) {
          
          coef.h.draws <- draws[,1:nPar, drop = FALSE]
          coef.g.draws <- draws[,(nPar+1):(2*nPar), drop = FALSE]
          
          if(x$floor==TRUE) {
            coef.ql.draws <- draws[,(2*nPar+1):(3*nPar), drop = FALSE]
            if(x$ceiling==TRUE){
              coef.qu.draws <- draws[,(3*nPar+1):(4*nPar), drop = FALSE] 
            } 
          } else {
            if(x$ceiling==TRUE){
              coef.qu.draws <- draws[,(2*nPar+1):(3*nPar), drop = FALSE] 
            } 
          }
          
          liar.ceiling.sims <- liar.floor.sims <-
            liar.ceiling.pop.sims <- liar.floor.pop.sims <- rep(NA, x$n.draws)
          
          for (i in 1:x$n.draws) {
            if (x$floor==TRUE)
              qlX <- logistic(x$x %*% as.vector(coef.ql.draws[i, , drop = FALSE]))
            if (x$ceiling==TRUE)
              quX <- logistic(x$x %*% as.vector(coef.qu.draws[i, , drop = FALSE]))
            
            hX <- logistic(x$x %*% as.vector(coef.h.draws[i,, drop = FALSE]))
            gX <- logistic(x$x %*% as.vector(coef.g.draws[i,, drop = FALSE]))
            
            hX0 <- dbinom(x = 0, size = x$J, prob = hX, log = FALSE)
            hXJ <- dbinom(x = x$J, size = x$J, prob = hX, log = FALSE)
            
            if (x$ceiling==TRUE) {
              liar.ceiling.sims[i] <- sum(quX * hXJ * gX) / sum(hXJ * gX)
              liar.ceiling.pop.sims[i] <- (1/n) * sum(quX * hXJ * gX)
            }
            
            if (x$floor==TRUE) {
              liar.floor.sims[i] <- sum(qlX * hX0 * gX) / sum(hX0 * gX)
              liar.floor.pop.sims[i] <- (1/n) * sum(qlX * hX0 * gX)
            }
            
          }
          
          if(x$ceiling==TRUE | x$floor == TRUE)
            cat("Quasi-Bayesian approximation estimates\n")
          
          if(x$ceiling==TRUE) {
            liar.ceiling <- mean(liar.ceiling.sims)
            liar.ceiling.se <- sd(liar.ceiling.sims)/sqrt(nrow(x$x))
            liar.ceiling.pop <- mean(liar.ceiling.pop.sims)
            liar.ceiling.pop.se <- sd(liar.ceiling.pop.sims)/sqrt(nrow(x$x))
            cat("Ceiling liar cond. prob. est. = ", round(liar.ceiling,5), " (s.e. = ",
                round(liar.ceiling.se,5), ")\nCeiling liar pop. prop. est. = ", round(liar.ceiling.pop,5), " (s.e. = ",
                round(liar.ceiling.pop.se,5), ")\n", sep = "")
          } 
          
          if(x$floor==TRUE) {
            liar.floor <- mean(liar.floor.sims)
            liar.floor.se <- sd(liar.floor.sims)/sqrt(nrow(x$x))
            liar.floor.pop <- mean(liar.floor.pop.sims)
            liar.floor.pop.se <- sd(liar.floor.pop.sims)/sqrt(nrow(x$x))
            cat("Floor liar cond. prob. est. = ",round(liar.floor,5), " (s.e. = ",
                round(liar.floor.se,5), ")\nFloor liar pop. prop. est. = ",round(liar.floor.pop,5), " (s.e. = ",
                round(liar.floor.pop.se,5), ")\n", sep = "")
          } 
          
        ##} else {
        }
        
          ## cannot use mvrnorm()
          
          coef.h <- x$par.control
          coef.g <- x$par.treat
          
          if(x$floor == TRUE) {
            coef.ql <- x$par.floor
            qlX <- logistic(x$x %*% coef.ql)
          }
          if(x$ceiling == TRUE) {
            coef.qu <- x$par.ceiling
            quX <- logistic(x$x %*% coef.qu)
          }
          
          hX <- logistic(x$x %*% coef.h)
          gX <- logistic(x$x %*% coef.g)
          
          hX0 <- dbinom(x = 0, size = x$J, prob = hX, log = FALSE)
          hXJ <- dbinom(x = x$J, size = x$J, prob = hX, log = FALSE)
          
           if(x$ceiling == TRUE | x$floor == TRUE)
          cat("\nMaximum likelihood estimates\n")
    
          if (x$ceiling == TRUE) {
            liar.ceiling <- sum(quX * hXJ * gX) / sum(hXJ * gX)
            liar.ceiling.pop <- (1/n) * sum(quX * hXJ * gX)
            cat("Ceiling liar cond. prob. est. =", round(liar.ceiling,5), "\n")
            cat("Ceiling liar pop. prop. est. =", round(liar.ceiling.pop,5),"\n")
          }
          if (x$floor == TRUE) {
            liar.floor <- sum(qlX * hX0 * gX) / sum(hX0 * gX)
            liar.floor.pop <- (1/n) * sum(qlX * hX0 * gX)
            cat("Floor liar cond. prob. est. =", round(liar.floor,5), "\n")
            cat("Floor liar pop. prop. est. =", round(liar.floor.pop,5),"\n")
          }
          if(!exists("draws")) {
          cat("\nNote: covariance matrix was not computationally singular,\n")
          cat("so std. errors for the proportion estimates could not be calculated.\n\n")
          } else {
          cat("\n")
          }
        }
        
      }
    
  } else {

    ## modified design
    
    cat("Sensitive Item \n", sep = "")
    tb.control <- matrix(NA, ncol = 2, nrow = length(x$par.control))
    tb.treat <- matrix(NA, ncol = 2, nrow = length(x$par.treat))
    colnames(tb.control) <- colnames(tb.treat) <- c("Est.", "S.E.")
    rownames(tb.control) <- rep(x$coef.names, x$J)
    rownames(tb.treat) <- x$coef.names
    
    tb.control[,1] <- x$par.control
    tb.control[,2] <- x$se.control

    tb.treat[,1] <- x$par.treat
    tb.treat[,2] <- x$se.treat
    
    print(as.matrix(round(tb.treat,5)))
    
    cat("Control Items \n", sep="")
    print(as.matrix(round(tb.control,5)))

    if (x$method == "nls")
      summ.stat <- paste("Residual standard error:",signif(x$resid.se, digits=6),
                         "with",x$resid.df,"degrees of freedom")
    else if (x$method == "ml")
      summ.stat <- paste("Log-likelihood:", round(x$llik,5))
    
    cat("\n",summ.stat,"\n\n", sep="")
      
  }


  treat.print <- c()
  for (i in 1:length(x$treat.labels)) {
    treat.print <- c(treat.print, "'", x$treat.labels[i], "'", sep = "")
    if (i != length(x$treat.labels))
      treat.print <- c(treat.print, " and ")
  }
  
  cat("Number of control items J set to ", x$J, ". Treatment groups were indicated by ", sep = "")
  cat(treat.print, sep ="")
  cat(" and the control group by '", x$control.label, "'.\n\n", sep = "")
   
  invisible(x)
  
}
