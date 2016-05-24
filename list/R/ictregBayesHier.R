
ictregBayesHier <- function(formula, data = parent.frame(), group.level.2, group.level.3, group.level.4, formula.level.2, formula.level.3, formula.level.4, treat="treat", J, fit.start = "lm", n.draws = 10000, burnin = 5000, thin = 0, delta.start.level.1, delta.mu0.level.1, delta.A0.level.1, delta.start.level.2, delta.mu0.level.2, delta.A0.level.2, delta.start.level.3, delta.mu0.level.3, delta.A0.level.3, delta.start.level.4, delta.mu0.level.4, delta.A0.level.4, sigma.start.level.1, sigma.df.level.1, sigma.scale.level.1, sigma.start.level.2, sigma.df.level.2, sigma.scale.level.2, sigma.start.level.3, sigma.df.level.3, sigma.scale.level.3, sigma.start.level.4, sigma.df.level.4, sigma.scale.level.4, delta.tune, alpha.tune, ##ceiling = FALSE, floor = FALSE,
                            verbose = TRUE, ...){

  ictreg.call <- match.call()

  ## removed to fix error in package compiling
  ##require(magic)
  
  ## set up data frame, with support for standard and modified responses
  mf <- match.call(expand.dots = FALSE)
  
  ## make all other call elements null in mf <- NULL in next line
  mf$group.level.2 <- mf$group.level.3 <- mf$group.level.4 <- mf$formula.level.2 <- mf$formula.level.3 <- mf$formula.level.4 <- mf$treat <- mf$J <- mf$fit.start <- mf$n.draws <- mf$burnin <- mf$thin <- mf$delta.start <- mf$delta.mu0 <- mf$delta.A0 <- mf$delta.tune <- mf$delta.start.level.1 <- mf$delta.mu0.level.1 <- mf$delta.A0.level.1 <- mf$delta.start.level.2 <- mf$delta.mu0.level.2 <- mf$delta.A0.level.2 <- mf$delta.start.level.3 <- mf$delta.mu0.level.3 <- mf$delta.A0.level.3 <- mf$delta.start.level.4 <- mf$delta.mu0.level.4 <- mf$delta.A0.level.4 <- mf$sigma.start.level.1 <- mf$sigma.df.level.1 <- mf$sigma.scale.level.1 <- mf$sigma.start.level.2 <- mf$sigma.df.level.2 <- mf$sigma.scale.level.2 <- mf$sigma.start.level.3 <- mf$sigma.df.level.3 <- mf$sigma.scale.level.3 <- mf$sigma.start.level.4 <- mf$sigma.df.level.4 <- mf$sigma.scale.level.4 <- mf$alpha.tune <- mf$verbose <- NULL
  ## mf$ceiling <- mf$floor <- NULL

  mf[[1]] <- as.name("model.frame")
  mf$na.action <- 'na.pass'
  mf <- eval.parent(mf)
  
  ## define design, response data frames
  x.all <- model.matrix.default(attr(mf, "terms"), mf)
  y.all <- model.response(mf)
  
  if(sum(x.all[,1]==1) == nrow(x.all))
    x.all <- x.all[, -1, drop = FALSE]
  
  levels <- 4

  if (missing("formula.level.4") == TRUE) {
    formula.level.4 <- ~ 1
    levels <- 3
  }  
  
  z.all <- model.matrix(formula.level.4, data)

  if (missing("formula.level.3") == TRUE) {
    formula.level.3 <- ~ 1
    levels <- 2
  }
  
  w.all <- model.matrix(formula.level.3, data)

  if (missing("formula.level.2") == TRUE) {
    formula.level.2 <- ~ 1
    stop("You must have at least a two-level hierarchical model. Please specify the second level formula.")
  }
  ## x v w z
  v.all <- model.matrix(formula.level.2, data)
  
  if(sum(v.all[,1]==1) == nrow(v.all) & ncol(v.all) > 1 & levels != 2)
    v.all <- v.all[, -1, drop = FALSE]
  if(sum(w.all[,1]==1) == nrow(w.all) & ncol(w.all) > 1 & levels != 3) 
    w.all <- w.all[, -1, drop = FALSE]
  if(sum(z.all[,1]==1) == nrow(z.all) & ncol(z.all) > 1)
    z.all <- z.all[, -1, drop = FALSE]
  
  # list-wise missing deletion
  na.x <- apply(is.na(as.matrix(x.all)), 1, sum)
  na.y <- is.na(y.all)
  
  na.v <- apply(is.na(as.matrix(v.all)), 1, sum)
  na.w <- apply(is.na(as.matrix(w.all)), 1, sum)
  na.z <- apply(is.na(as.matrix(z.all)), 1, sum)
  na.cond <- na.x==0 & na.y==0 & na.v==0 & na.w==0 & na.z==0
 
  ## group indicator for mixed effects regression
  grp2 <- data[na.cond == TRUE, paste(group.level.2)]
  if (levels >= 3) grp3 <- data[na.cond == TRUE, paste(group.level.3)]
  if (levels == 4) grp4 <- data[na.cond == TRUE, paste(group.level.4)]
  
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
    ##condition.labels <- as.character(paste("Sensitive Item",sort(unique(t))))
    ##condition.labels[which(condition.labels == "Sensitive Item 0")] <- "Control Items"
    ##treatment.labels <- condition.labels[condition.labels != "Control Items"]
    ##control.label <- "Control Items"
    condition.labels <- sort(unique(t))
    treatment.labels <- condition.labels[condition.labels != 0]
    control.label <- 0
  }
  
  ## list wise delete
  y.all <- y.all[na.cond == TRUE]
  x.all <- x.all[na.cond == TRUE, , drop = FALSE]
  v.all <- v.all[na.cond == TRUE, , drop = FALSE]
  w.all <- w.all[na.cond == TRUE, , drop = FALSE]
  z.all <- z.all[na.cond == TRUE, , drop = FALSE]

  ## now chop the covariate dataframes for the hierarchical levels so they are
  ## one row per group, not all indivs duplicated
  v.all <- na.omit(aggregate(v.all, by = list(grp2), FUN = mean))[,-1,drop=F]
  if (levels >= 3) w.all <- na.omit(aggregate(w.all, by = list(grp3), FUN = mean))[,-1,drop=F]
  if (levels == 4) z.all <- na.omit(aggregate(z.all, by = list(grp4), FUN = mean))[,-1,drop=F]
    
  ## set up data objects for y and x for each group from user input
  x.treatment <- x.all[t != 0, , drop = FALSE]
  y.treatment <- subset(y.all, t != 0)
  x.control <- x.all[t == 0 , , drop = FALSE]
  y.control <- subset(y.all, t==0)

  ##v.treatment <- v.all[t != 0, , drop = FALSE]
  ##v.control <- v.all[t == 0 , , drop = FALSE]
  ##w.treatment <- w.all[t != 0, , drop = FALSE]
  ##w.control <- w.all[t == 0 , , drop = FALSE]
  ##z.treatment <- z.all[t != 0, , drop = FALSE]
  ##z.control <- z.all[t == 0 , , drop = FALSE]

  ## J is defined by user for the standard design
  if(missing("J")) {
    J <- max(y.treatment) - 1
  }

  condition.values <- sort(unique(t))
  treatment.values <- 1:length(condition.values[condition.values!=0])
 
  n <- nrow(x.treatment) + nrow(x.control)

  n.grp2 <- length(unique(grp2))
  if (levels >= 3) n.grp3 <- length(unique(grp3))
  if (levels == 4) n.grp4 <- length(unique(grp4))

  coef.names <- colnames(x.all)

  coef.names.level.2 <- colnames(v.all)
  coef.names.level.3 <- colnames(w.all)
  coef.names.level.4 <- colnames(z.all)

  print(summary(v.all))
  
  nPar.level.1 <- ncol(x.all)
  nPar.level.2 <- ncol(v.all)
  nPar.level.3 <- ncol(w.all)
  nPar.level.4 <- ncol(z.all)

  ## now switch the group indicators to be of the shorter length - of length the # of groups in the level below
  if (levels >= 3) {
    grp3.long <- grp3
    grp3 <- unique(cbind(grp2, grp3))[order(unique(cbind(grp2, grp3))[,1]),][,2]
  }
  if (levels == 4) {
    grp4.long <- grp4
    grp4 <- unique(cbind(grp3.long, grp4))[order(unique(cbind(grp3.long, grp4))[,1]),][,2]
  }
  
  intercept.only <- ncol(x.all)==1 & sum(x.all[,1]==1) == n

  if (intercept.only == TRUE) {
    x.vars <- "1"
  } else {
    x.vars <- coef.names[-1]
  }
  
  ## NB: MAY NEED TO CHECK THIS
  intercept.only.level.2 <- ncol(v.all)==1 & sum(v.all[,1]==1) == n
  
  if (intercept.only.level.2 == TRUE) {
    x.vars.level.2 <- "1"
  } else {
    if (levels == 2)
      x.vars.level.2 <- coef.names.level.2[-1]
    else
      x.vars.level.2 <- coef.names.level.2
  }

  intercept.only.level.3 <- ncol(w.all)==1 & sum(w.all[,1]==1) == n
  
  if (intercept.only.level.3 == TRUE) {
    x.vars.level.3 <- "1"
  } else {
    if (levels == 3)
      x.vars.level.3 <- coef.names.level.3[-1]
    else
      x.vars.level.3 <- coef.names.level.3
  }

  intercept.only.level.4 <- ncol(z.all)==1 & sum(z.all[,1]==1) == n
  
  if (intercept.only.level.4 == TRUE) {
    x.vars.level.4 <- "1"
  } else {
    if (levels == 4)
      x.vars.level.4 <- coef.names.level.4[-1]
    else
      x.vars.level.4 <- coef.names.level.4
  }
  
  logistic <- function(x) exp(x)/(1+exp(x))
  logit <- function(x) return(log(x)-log(1-x))

  if (missing("delta.tune")) {
    stop("The Metropolis tuning input object delta.tune is required.")
  }
    
  ## multi code
 
  ##if (length(ceiling) == 1)
  ##  ceiling <- rep(ceiling, length(treatment.values))
  ##if (length(floor) == 1)
  ##  floor <- rep(floor, length(treatment.values))
 
  ## mixed multi model
  
  ##if (missing("delta.start.level.1"))
  ##  delta.start.level.1 = rep(0, nPar.level.1)

  if (missing("delta.start.level.1")) {
    ictreg.fit <- ictreg(formula, data = data, treat = treat, J = J, method = "ml", multi.condition = "none")
    
    if(sum(x.all[,1]==1) == n)
      subset <- -1
    else
      subset <- 1:nPar.level.1
    
    delta.start <- list("0" = ictreg.fit$par.control[subset])
    for (m in 1:length(treatment.values))
      delta.start[[as.character(m)]] <- ictreg.fit$par.treat[[m]][subset]
  }
  
  if (missing("delta.start.level.2"))
    delta.start.level.2 = rep(0, nPar.level.2)
  if (missing("delta.start.level.3"))
    delta.start.level.3 = rep(0, nPar.level.3)
  if (missing("delta.start.level.4"))
    delta.start.level.4 = rep(0, nPar.level.4)
  
  if (missing("sigma.start.level.1"))
    sigma.start.level.1 = rep(1, length(treatment.values)+1) 
  if (missing("sigma.start.level.2"))
    sigma.start.level.2 = rep(1, length(treatment.values)+1)
  if (missing("sigma.start.level.3"))
    sigma.start.level.3 = rep(1, length(treatment.values)+1)
  if (missing("sigma.start.level.4"))
    sigma.start.level.4 = rep(1, length(treatment.values)+1)

  print(paste("npar level 1", nPar.level.1))

  ##if (missing("delta.mu0")) {
  ##  delta.mu0 <- rep(0, nPar.level.1)
  ##}
  ##if (missing("delta.A0")) {
  ##  delta.A0 <- list(diag(1, nPar.level.1), diag(1, nPar.level.1), diag(1, nPar.level.1))
  ##}
  
  if (missing("delta.mu0.level.1"))
    delta.mu0.level.1 = rep(0, nPar.level.1)
  if (missing("delta.A0.level.1"))
    delta.A0.level.1 = list(diag(1, nPar.level.1), diag(1, nPar.level.1), diag(1, nPar.level.1))
  if (missing("sigma.df.level.1"))
    sigma.df.level.1 = 5
  if (missing("sigma.scale.level.1"))
    sigma.scale.level.1 = 1
  
  if (missing("delta.mu0.level.2"))
    delta.mu0.level.2 = rep(0, nPar.level.2)
  if (missing("delta.A0.level.2"))
    delta.A0.level.2 = list(diag(1, nPar.level.2), diag(1, nPar.level.2), diag(1, nPar.level.2))
  if (missing("sigma.df.level.2"))
    sigma.df.level.2 = 5
  if (missing("sigma.scale.level.2"))
    sigma.scale.level.2 = 1

  if (missing("delta.mu0.level.3"))
    delta.mu0.level.3 = rep(0, nPar.level.3)
  if (missing("delta.A0.level.3"))
    delta.A0.level.3 = list(diag(1, nPar.level.3), diag(1, nPar.level.3), diag(1, nPar.level.3))
  if (missing("sigma.df.level.3"))
    sigma.df.level.3 = 5
  if (missing("sigma.scale.level.3"))
    sigma.scale.level.3 = 1

if (missing("delta.mu0.level.4"))
    delta.mu0.level.4 = rep(0, nPar.level.4)
  if (missing("delta.A0.level.4"))
    delta.A0.level.4 = list(diag(1, nPar.level.4), diag(1, nPar.level.4), diag(1, nPar.level.4))
  if (missing("sigma.df.level.4"))
    sigma.df.level.4 = 5
  if (missing("sigma.scale.level.4"))
    sigma.scale.level.4 = 1
  
  ##if (missing("delta.tune"))
  ##  delta.tune = rep(0.001, nPar.level.1)
  
  if (missing("alpha.tune"))
    alpha.tune = rep(0.001, n.grp2)

  if (levels==2) {

##delta.grp.mu0 = 0, delta.grp.A0 = diag(1), 
##delta.grp.mu0 = delta.mu0.level.1, delta.grp.A0 = delta.A0.level.1,
    
    fit <- ictregBayesMulti2Level.fit(Y = y.all, treat = t, X = as.matrix(x.all), V = as.matrix(v.all), J = J, grp = as.numeric(grp2),
                                      ##ceiling = ceiling, floor = floor,
                                      n.draws = n.draws, burnin = burnin, thin = thin, verbose = verbose, delta.start = delta.start.level.1,
                                      delta.grp.start = delta.start.level.2, sigma.start = sigma.start.level.2, delta.mu0 = delta.mu0.level.1, 
                                      delta.A0 = delta.A0.level.1, delta.grp.mu0 = delta.mu0.level.2, delta.grp.A0 = delta.A0.level.2,
                                      sigma.df = sigma.df.level.1, sigma.scale = sigma.scale.level.1,
                                      delta.tune = delta.tune, alpha.tune = alpha.tune, ...)

  } else if (levels==3){

    fit <- ictregBayesMulti3Level.fit(Y = y.all, treat = t, X = as.matrix(x.all), V = as.matrix(v.all), W = as.matrix(w.all), J = J,
                                      grp1 = as.numeric(grp2), grp2 = as.numeric(grp3),
                                      ##ceiling = ceiling, floor = floor,
                                      n.draws = n.draws, burnin = burnin, thin = thin, verbose = verbose, delta.start = delta.start.level.1,
                                      delta.grp1.start = delta.start.level.2, delta.grp2.start = delta.start.level.3,
                                      sigma.grp1.start = sigma.start.level.2, sigma.grp2.start = sigma.start.level.3,
                                      delta.mu0 = delta.mu0.level.1, delta.A0 = delta.A0.level.1,
                                      delta.grp1.mu0 = delta.mu0.level.2, delta.grp2.mu0 = delta.mu0.level.3,
                                      delta.grp1.A0 = delta.A0.level.2, delta.grp2.A0 = delta.A0.level.3,
                                      sigma.grp1.df = sigma.df.level.2, sigma.grp1.scale = sigma.scale.level.2,
                                      sigma.grp2.df = sigma.df.level.3, sigma.grp2.scale = sigma.scale.level.3,
                                      delta.tune = delta.tune, alpha.tune = alpha.tune, ...)
  } else if (levels==4){

    fit <- ictregBayesMulti4Level.fit(Y = y.all, treat = t, X = as.matrix(x.all), V = as.matrix(v.all), W = as.matrix(w.all), Z = as.matrix(z.all), J = J,
                                      grp1 = as.numeric(grp2), grp2 = as.numeric(grp3), grp3 = as.numeric(grp4),
                                     ## ceiling = ceiling, floor = floor,
                                      n.draws = n.draws, burnin = burnin, thin = thin, verbose = verbose, delta.start = delta.start.level.1,
                                      delta.grp1.start = delta.start.level.2, delta.grp2.start = delta.start.level.3,
                                      delta.grp3.start = delta.start.level.4,
                                      sigma.grp1.start = sigma.start.level.2, sigma.grp2.start = sigma.start.level.3,
                                      sigma.grp3.start = sigma.start.level.4,
                                      delta.mu0 = delta.mu0.level.1, delta.grp1.mu0 = delta.mu0.level.2, delta.grp2.mu0 = delta.mu0.level.3,
                                      delta.grp3.mu0 = delta.mu0.level.4,
                                      delta.A0 = delta.A0.level.1, delta.grp1.A0 = delta.A0.level.2, delta.grp2.A0 = delta.A0.level.3,
                                      delta.grp3.A0 = delta.A0.level.4,
                                      sigma.grp1.df = sigma.df.level.2, sigma.grp1.scale = sigma.scale.level.2,
                                      sigma.grp2.df = sigma.df.level.3, sigma.grp2.scale = sigma.scale.level.3,
                                      sigma.grp3.df = sigma.df.level.4, sigma.grp3.scale = sigma.scale.level.4,
                                      delta.tune = delta.tune, alpha.tune = alpha.tune, ...)
   }
    
  rownames(fit$delta) <- rownames(fit$alpha) <- rownames(fit$delta.grp) <- rownames(fit$sigma)
      seq(from = n.draws - burnin + 1 + thin + 1,
          to = n.draws - burnin + thin + 1 + floor((n.draws - burnin)/(thin + 1)) * (thin + 1), by = thin + 1)                   
  ##delta <- list()
  ##for (m in 1:length(treatment.labels)) {
  ##      delta[[m]] <- mcmc(data = fit$delta[, ((m-1)*nPar + 1):(m*nPar)],
  ##                         start = 1, thin = 1, end = nrow(fit$delta))
  ##    }
  ##fit$delta <- delta
  ##fit$psi <- mcmc(data = fit$psi, start = 1, thin = 1, end = nrow(fit$psi))
  ##fit$Sigma <- mcmc(data = fit$Sigma, start = 1, thin = 1, end = nrow(fit$Sigma))
  ##fit$Phi <- mcmc(data = fit$Phi, start = 1, thin = 1, end = length(fit$Phi))
  ##fit$gamma <- mcmc(data = fit$gamma, start = 1, thin = 1, end = nrow(fit$gamma))
  ##fit$zeta <- mcmc(data = fit$zeta, start = 1, thin = 1, end = nrow(fit$zeta))
  
  fit$x <- x.all
  fit$delta.start.level.1 <- delta.start.level.1
  fit$delta.mu0 <- delta.mu0.level.1
  fit$delta.A0 <- delta.A0.level.1
  fit$delta.tune <- delta.tune
  fit$J <- J
  fit$treat.labels <- treatment.labels
  fit$control.label <- control.label
  fit$levels <- levels
  fit$coef.names <- coef.names
  fit$coef.names.level.2 <- coef.names.level.2
  if(levels >= 3) fit$coef.names.level.3 <- coef.names.level.3
  if(levels >= 4) fit$coef.names.level.4 <- coef.names.level.4
  fit$data.level.1 <- as.matrix(x.all)
  fit$data.level.2 <- as.matrix(v.all)
  if(levels >= 3) fit$data.level.3 <- as.matrix(w.all)
  if(levels >= 4) fit$data.level.4 <- as.matrix(z.all)
  fit$group.level.2 <- grp2
  if(levels >= 3) fit$group.level.3 <- grp3
  if(levels >= 4) fit$group.level.4 <- grp4
  fit$call <- match.call()
  fit$random.seed <- .Random.seed

  if(levels==2)
    names(fit)[2:4] <- c("alpha.grp1","delta.grp1","sigma.grp1")
  
  class(fit) <- "ictregBayesHier"
  
  return(fit)
  
}

as.list.ictregBayesHier <- function(...) {
  
  x <- list(...)

  delta.list <- list()
  for (i in 1:length(x))
    delta.list[[i]] <- x[[i]]$delta
  
  alpha.grp1.list <- list()
  for (i in 1:length(x))
    alpha.grp1.list[[i]] <- x[[i]]$alpha.grp1

  delta.grp1.list <- list()
  for (i in 1:length(x))
    delta.grp1.list[[i]] <- x[[i]]$delta.grp1

  sigma.grp1.list <- list()
  for (i in 1:length(x))
    sigma.grp1.list[[i]] <- x[[i]]$sigma.grp1

  if (x[[1]]$levels >= 3) {
    alpha.grp2.list <- list()
    for (i in 2:length(x))
      alpha.grp2.list[[i]] <- x[[i]]$alpha.grp2
    
    delta.grp2.list <- list()
    for (i in 2:length(x))
      delta.grp2.list[[i]] <- x[[i]]$delta.grp2
    
    sigma.grp2.list <- list()
    for (i in 2:length(x))
      sigma.grp2.list[[i]] <- x[[i]]$sigma.grp2
  }
  
  if (x[[1]]$levels >= 4) {
    alpha.grp3.list <- list()
    for (i in 3:length(x))
      alpha.grp3.list[[i]] <- x[[i]]$alpha.grp3
    
    delta.grp3.list <- list()
    for (i in 3:length(x))
      delta.grp3.list[[i]] <- x[[i]]$delta.grp3
    
    sigma.grp3.list <- list()
    for (i in 3:length(x))
      sigma.grp3.list[[i]] <- x[[i]]$sigma.grp3
  }
  
  return.object <- x[[1]]
  return.object$delta <- delta.list
  return.object$alpha.grp1 <- alpha.grp1.list
  return.object$delta.grp1 <- delta.grp1.list
  return.object$sigma.grp1 <- sigma.grp1.list
  if (return.object$levels >= 3) {
    return.object$alpha.grp2 <- alpha.grp2.list
    return.object$delta.grp2 <- delta.grp2.list
    return.object$sigma.grp2 <- sigma.grp2.list
  }
  if (return.object$levels == 4) {
    return.object$alpha.grp3 <- alpha.grp3.list
    return.object$delta.grp3 <- delta.grp3.list
    return.object$sigma.grp3 <- sigma.grp3.list
  }
  class(return.object) <- "ictregBayesHier.list"

  return.object
  
}

coef.ictregBayesHier.list <- function(object, ranef = FALSE, ...) {

  object$delta <- as.mcmc(do.call(rbind, object$delta))
  object$delta.grp1 <- as.mcmc(do.call(rbind, object$delta.grp1))
  if(object$levels >= 3)
    object$delta.grp2 <- as.mcmc(do.call(rbind, object$delta.grp2))
  if(object$levels == 4)
    object$delta.grp3 <- as.mcmc(do.call(rbind, object$delta.grp3))

  class(object) <- "ictregBayesHier"

  coef(object, ranef = ranef, ... = ...)

}

coef.ictregBayesHier <- function(object, ranef = FALSE, ...) {

  M <- length(object$treat.labels)
  nPar <- length(object$coef.names)
  nPar.level.2 <- length(object$coef.names.level.2)
  if (object$levels>=3)
    nPar.level.3 <- length(object$coef.names.level.3)
  if (object$levels==4)
    nPar.level.4 <- length(object$coef.names.level.4)
  
  delta.coef <- delta.coef.level.2 <- delta.coef.level.3 <- delta.coef.level.4 <- list()
  for (m in 1:M) {
    delta.coef[[object$treat.labels[[m]]]] <- apply(object$delta[,( (m-1) * nPar + 1) : (m*nPar), drop = FALSE], 2, mean)
    names(delta.coef[[object$treat.labels[[m]]]]) <- object$coef.names
    delta.coef.level.2[[object$treat.labels[[m]]]] <- apply(object$delta.grp1[,( (m-1) * nPar.level.2 + 1) : (m*nPar.level.2), drop = FALSE], 2, mean)
    names(delta.coef.level.2[[object$treat.labels[[m]]]]) <- object$coef.names.level.2
    if (object$levels>=3) {
      delta.coef.level.3[[object$treat.labels[[m]]]] <- apply(object$delta.grp2[,( (m-1) * nPar.level.3 + 1) : (m*nPar.level.3), drop = FALSE], 2, mean)
      names(delta.coef.level.3[[object$treat.labels[[m]]]]) <- object$coef.names.level.3
    }
    if (object$levels==4) {
      delta.coef.level.4[[object$treat.labels[[m]]]] <- apply(object$delta.grp3[,( (m-1) * nPar.level.4 + 1) : (m*nPar.level.4), drop = FALSE], 2, mean)
      names(delta.coef.level.4[[object$treat.labels[[m]]]]) <- object$coef.names.level.4
    }
  }
  
  psi.coef <- apply(object$delta[,( M * nPar + 1) : ((M+1)*nPar), drop = FALSE], 2, mean)
  names(psi.coef) <- object$coef.names

  psi.coef.level.2 <- apply(object$delta.grp1[,( M * nPar.level.2 + 1) : ((M+1)*nPar.level.2), drop = FALSE], 2, mean)
  names(psi.coef.level.2) <- object$coef.names.level.2
  
  return.object <- list(delta = delta.coef, psi = psi.coef, delta.level.2 = delta.coef.level.2, psi.level.2 = psi.coef.level.2)

  if (object$levels>=3) {
    return.object$delta.level.3 <- delta.coef.level.3
    return.object$psi.level.3 <- apply(object$delta.grp2[,( M * nPar.level.3 + 1) : ((M+1)*nPar.level.3), drop = FALSE], 2, mean)
    names(return.object$psi.level.3) <- object$coef.names.level.3
  }
  
  if (object$levels==4) {
    return.object$delta.level.4 <- delta.coef.level.4
    return.object$psi.level.4 <- apply(object$delta.grp3[,( M * nPar.level.4 + 1) : ((M+1)*nPar.level.4), drop = FALSE], 2, mean)
    names(return.object$psi.level.4) <- object$coef.names.level.4
  }
  
  if (ranef == TRUE) {

    n.grp2 <- ncol(object$alpha.grp1)/3
    if(object$levels>=3) n.grp3 <- ncol(object$alpha.grp2)/3
    if(object$levels==4) n.grp4 <- ncol(object$alpha.grp3)/3
    
    gamma.coef.level.2 <- gamma.coef.level.3 <- gamma.coef.level.4 <- list()
    for(m in 1:M) {
      gamma.coef.level.2[[object$treat.labels[[m]]]] <- apply(object$alpha.grp1[,( (m-1) * n.grp2 + 1) : (m*n.grp2), drop = FALSE], 2, mean)
      if(object$levels>=3) gamma.coef.level.3[[object$treat.labels[[m]]]] <- apply(object$alpha.grp2[,( (m-1) * n.grp3 + 1) : (m*n.grp3), drop = FALSE], 2, mean)
      if(object$levels==4) gamma.coef.level.4[[object$treat.labels[[m]]]] <- apply(object$alpha.grp3[,( (m-1) * n.grp4 + 1) : (m*n.grp4), drop = FALSE], 2, mean)
    }
    
    return.object$ranef.gamma.level.2 <- gamma.coef.level.2
    return.object$ranef.zeta.level.2 <- apply(object$alpha.grp1[,( M * n.grp2 + 1) : ((M+1)*n.grp2), drop = FALSE], 2, mean)
    
    if(object$levels>=3) {
      return.object$ranef.gamma.level.3 <- gamma.coef.level.3
      return.object$ranef.zeta.level.3 <- apply(object$alpha.grp2[,( M * n.grp3 + 1) : ((M+1)*n.grp3), drop = FALSE], 2, mean)
    }
    if(object$levels==4) {
      return.object$ranef.gamma.level.4 <- gamma.coef.level.4
      return.object$ranef.zeta.level.4 <- apply(object$alpha.grp3[,( M * n.grp4 + 1) : ((M+1)*n.grp4), drop = FALSE], 2, mean)
    }
    
  }
  
  return.object

}

sd.ictregBayesHier.list <- function(object, ranef = FALSE, ...) {

  object$delta <- as.mcmc(do.call(rbind, object$delta))
  object$delta.grp1 <- as.mcmc(do.call(rbind, object$delta.grp1))
  if(object$levels >= 3)
    object$delta.grp2 <- as.mcmc(do.call(rbind, object$delta.grp2))
  if(object$levels == 4)
    object$delta.grp3 <- as.mcmc(do.call(rbind, object$delta.grp3))

  class(object) <- "ictregBayesHier"

  sd.ictregBayesHier(object, ranef = ranef, ... = ...)

}

sd.ictregBayesHier <- function(object, ranef = FALSE, ...) {

  M <- length(object$treat.labels)
  nPar <- length(object$coef.names)
  nPar.level.2 <- length(object$coef.names.level.2)
  if (object$levels>=3)
    nPar.level.3 <- length(object$coef.names.level.3)
  if (object$levels==4)
    nPar.level.4 <- length(object$coef.names.level.4)
  
  delta.coef <- delta.coef.level.2 <- delta.coef.level.3 <- delta.coef.level.4 <- list()
  for (m in 1:M) {
    delta.coef[[object$treat.labels[[m]]]] <- apply(object$delta[,( (m-1) * nPar + 1) : (m*nPar), drop = FALSE], 2, sd)
    names(delta.coef[[object$treat.labels[[m]]]]) <- object$coef.names
    delta.coef.level.2[[object$treat.labels[[m]]]] <- apply(object$delta.grp1[,( (m-1) * nPar.level.2 + 1) : (m*nPar.level.2), drop = FALSE], 2, sd)
    names(delta.coef.level.2[[object$treat.labels[[m]]]]) <- object$coef.names.level.2
    if (object$levels>=3) {
      delta.coef.level.3[[object$treat.labels[[m]]]] <- apply(object$delta.grp2[,( (m-1) * nPar.level.3 + 1) : (m*nPar.level.3), drop = FALSE], 2, sd)
      names(delta.coef.level.3[[object$treat.labels[[m]]]]) <- object$coef.names.level.3
    }
    if (object$levels==4) {
      delta.coef.level.4[[object$treat.labels[[m]]]] <- apply(object$delta.grp3[,( (m-1) * nPar.level.4 + 1) : (m*nPar.level.4), drop = FALSE], 2, sd)
      names(delta.coef.level.4[[object$treat.labels[[m]]]]) <- object$coef.names.level.4
    }
  }
  
  psi.coef <- apply(object$delta[,( M * nPar + 1) : ((M+1)*nPar), drop = FALSE], 2, sd)
  names(psi.coef) <- object$coef.names

  psi.coef.level.2 <- apply(object$delta.grp1[,( M * nPar.level.2 + 1) : ((M+1)*nPar.level.2), drop = FALSE], 2, sd)
  names(psi.coef.level.2) <- object$coef.names.level.2
  
  return.object <- list(delta = delta.coef, psi = psi.coef, delta.level.2 = delta.coef.level.2, psi.level.2 = psi.coef.level.2)

  if (object$levels>=3) {
    return.object$delta.level.3 <- delta.coef.level.3
    return.object$psi.level.3 <- apply(object$delta.grp2[,( M * nPar.level.3 + 1) : ((M+1)*nPar.level.3), drop = FALSE], 2, sd)
    names(return.object$psi.level.3) <- object$coef.names.level.3
  }
  
  if (object$levels==4) {
    return.object$delta.level.4 <- delta.coef.level.4
    return.object$psi.level.4 <- apply(object$delta.grp3[,( M * nPar.level.4 + 1) : ((M+1)*nPar.level.4), drop = FALSE], 2, sd)
    names(return.object$psi.level.4) <- object$coef.names.level.4
  }
  
  if (ranef == TRUE) {

    n.grp2 <- ncol(object$alpha.grp1)/3
    if(object$levels>=3) n.grp3 <- ncol(object$alpha.grp2)/3
    if(object$levels==4) n.grp4 <- ncol(object$alpha.grp3)/3
    
    gamma.coef.level.2 <- gamma.coef.level.3 <- gamma.coef.level.4 <- list()
    for(m in 1:M) {
      gamma.coef.level.2[[object$treat.labels[[m]]]] <- apply(object$alpha.grp1[,( (m-1) * n.grp2 + 1) : (m*n.grp2), drop = FALSE], 2, sd)
      if(object$levels>=3) gamma.coef.level.3[[object$treat.labels[[m]]]] <- apply(object$alpha.grp2[,( (m-1) * n.grp3 + 1) : (m*n.grp3), drop = FALSE], 2, sd)
      if(object$levels==4) gamma.coef.level.4[[object$treat.labels[[m]]]] <- apply(object$alpha.grp3[,( (m-1) * n.grp4 + 1) : (m*n.grp4), drop = FALSE], 2, sd)
    }
    
    return.object$ranef.gamma.level.2 <- gamma.coef.level.2
    return.object$ranef.zeta.level.2 <- apply(object$alpha.grp1[,( M * n.grp2 + 1) : ((M+1)*n.grp2), drop = FALSE], 2, sd)
    
    if(object$levels>=3) {
      return.object$ranef.gamma.level.3 <- gamma.coef.level.3
      return.object$ranef.zeta.level.3 <- apply(object$alpha.grp2[,( M * n.grp3 + 1) : ((M+1)*n.grp3), drop = FALSE], 2, sd)
    }
    if(object$levels==4) {
      return.object$ranef.gamma.level.4 <- gamma.coef.level.4
      return.object$ranef.zeta.level.4 <- apply(object$alpha.grp3[,( M * n.grp4 + 1) : ((M+1)*n.grp4), drop = FALSE], 2, sd)
    }
    
  }
  
  return.object

}

summary.ictregBayesHier <- function(object, ...) {
  structure(object, class = c("summary.ictregBayesHier", class(object)))
}

print.summary.ictregBayesHier <- function(x, ...) {
  
  cat("\nItem Count Technique Bayesian Hierarchical Regression \n\nCall: ")
  
  dput(x$call)

 ## if (x$multi == TRUE) {
    for (k in 1:length(x$treat.labels)) {
      cat(paste("\nSensitive item (", x$treat.labels[k], ")", "\n", sep = ""))
      print(matrix(c(round(cbind(coef(x)$delta[[k]], sd.ictregBayesHier(x)$delta[[k]]),5)),
                   nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
                   dimnames = list(x$coef.names, c("Est.", "S.E."))))
      cat("\nMetropolis acceptance ratio:", round(x$delta.accept[[k]], 3), "\n")    
    }
##  } else {
    
 ##   cat("\nSensitive item \n")
 ##   print(matrix(c(round(cbind(coef(x)$delta, sd.ictregBayes(x)$delta),5)), nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
 ##                dimnames = list(x$coef.names, c("Est.", "S.E."))))
 ##   cat("\nMetropolis acceptance ratio:", round(x$delta.accept, 3), "\n")
    
 ## }
    
  cat("\nControl items \n")
  print(matrix(c(round(cbind(coef(x)$psi, sd.ictregBayesHier(x)$psi),5)), nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
                  dimnames = list(x$coef.names, c("Est.", "S.E."))))
  ##cat("\nMetropolis acceptance ratio:", round(x$psi.accept, 3), "\n")

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

summary.ictregBayesHier.list <- function(object, ...) {
  structure(object, class = c("summary.ictregBayesHier.list", class(object)))
}

print.summary.ictregBayesHier.list <- function(x, ...) {
  
  cat("\nItem Count Technique Bayesian Hierarchical Regression \n\nCall: ")
  
  dput(x$call)

  cat("\nSummary from",length(x$delta),"chains\n")

  mcmc.list <- list()
  for (i in 1:length(x$delta))
    mcmc.list[[i]] <- as.mcmc(x$delta[[i]])
  
  mcmc.list <- as.mcmc.list(mcmc.list)
  gelmanrubin <- round(gelman.diag(mcmc.list)$psrf[,1],4)
  names(gelmanrubin) <- rep(x$coef.names, length(x$treat.labels)+1)

  ## if (x$multi == TRUE) {
  for (k in 1:length(x$treat.labels)) {
    cat(paste("\nSensitive item (", x$treat.labels[k], ")", "\n", sep = ""))
    print(matrix(c(round(cbind(coef(x)$delta[[k]], sd.ictregBayesHier.list(x)$delta[[k]]),5)),
                 nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
                 dimnames = list(x$coef.names, c("Est.", "S.E."))))
    cat("\nMetropolis acceptance ratio:", round(x$delta.accept[[k]], 3), "\n")
    
    cat("\nGelman-Rubin statistics:\n")
        
    print(gelmanrubin[( (k-1) * length(x$coef.names) + 1) : (k*length(x$coef.names))])
    
  }
##  } else {
 ##   cat("\nSensitive item \n")
 ##   print(matrix(c(round(cbind(coef(x)$delta, sd.ictregBayes(x)$delta),5)), nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
 ##                dimnames = list(x$coef.names, c("Est.", "S.E."))))
 ##   cat("\nMetropolis acceptance ratio:", round(x$delta.accept, 3), "\n")
    
 ## }
    
  cat("\nControl items \n")
  print(matrix(c(round(cbind(coef(x)$psi, sd.ictregBayesHier.list(x)$psi),5)), nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
                  dimnames = list(x$coef.names, c("Est.", "S.E."))))
  ##cat("\nMetropolis acceptance ratio:", round(x$psi.accept, 3), "\n")

  cat("\nGelman-Rubin statistics:\n")

  print(gelmanrubin[( (length(x$treat.labels)) * length(x$coef.names) + 1) : ((length(x$treat.labels)+1)*length(x$coef.names))])

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

print.ictregBayesHier <- print.ictregBayesHier.list <- function(x, ...) { 
  
  cat("\nItem Count Technique Bayesian Hierarchical Regression \n\nCall: ")
  
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

predict.ictregBayesHier.list <- function(object, ...) {

  object$delta <- as.mcmc(do.call(rbind, object$delta))
  object$delta.grp1 <- as.mcmc(do.call(rbind, object$delta.grp1))
  if(object$levels >= 3)
    object$delta.grp2 <- as.mcmc(do.call(rbind, object$delta.grp2))
  if(object$levels == 4)
    object$delta.grp3 <- as.mcmc(do.call(rbind, object$delta.grp3))

  class(object) <- "ictregBayesHier"

  predict(object, ... = ...)

}

predict.ictregBayesHier <- function(object, newdata, se.fit = FALSE,
                                interval = c("none","confidence"), level = .95, sensitive.item, ...){

  n.draws <- nrow(object$delta)
  
  if(missing(interval)) interval <- "none"

  if(missing(sensitive.item)) {
    sensitive.item <- 1
    warning("Using the first sensitive item for predictions. Change with the sensitive.item option.")
  }

  nPar <- length(object$coef.names)
   
  logistic <- function(object) exp(object)/(1+exp(object))
  
  if(missing(newdata)) {
    if (object$levels == 2) xvar <- list(object$data.level.1, object$data.level.2[object$group.level.2,])
    if (object$levels >= 3) xvar <- list(object$data.level.1, object$data.level.2[object$group.level.2,],
          object$data.level.3[object$group.level.3[object$group.level.2],])
    if (object$levels >= 4) xvar <- list(object$data.level.1, object$data.level.2[object$group.level.2,],
          object$data.level.3[object$group.level.3[object$group.level.2],],
          object$data.level.4[object$group.level.4[object$group.level.3[object$group.level.2]],])
  } else { 
    if(nrow(newdata[[1]])==0)
      stop("No data in the provided data frame.")
    ##xvar <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))), newdata)
  }

  M <- length(object$treat.labels)
  nPar <- length(object$coef.names)
  nPar.level.2 <- length(object$coef.names.level.2)
  if (object$levels>=3)
    nPar.level.3 <- length(object$coef.names.level.3)
  if (object$levels==4)
    nPar.level.4 <- length(object$coef.names.level.4)

  m <- sensitive.item
  
  delta.coef <- object$delta[,( (m-1) * nPar + 1) : (m*nPar), drop = FALSE]
  delta.coef.level.2 <- object$delta.grp1[,( (m-1) * nPar.level.2 + 1) : (m*nPar.level.2), drop = FALSE]
  if (object$levels>=3)
    delta.coef.level.3 <- object$delta.grp2[,( (m-1) * nPar.level.3 + 1) : (m*nPar.level.3), drop = FALSE]
  if (object$levels==4)
    delta.coef.level.4 <- object$delta.grp3[,( (m-1) * nPar.level.4 + 1) : (m*nPar.level.4), drop = FALSE]
  
  pred.list.mean <- rep(NA, n.draws)
  
  for (d in 1:n.draws) {
        
    if(object$levels==2)
      pred.list <- logistic(as.matrix(xvar[[1]]) %*% as.matrix(delta.coef[d, ]) + as.matrix(xvar[[2]]) %*% as.matrix(delta.coef.level.2[d, ]) )
    if(object$levels==3)
      pred.list <- logistic(as.matrix(xvar[[1]]) %*% as.matrix(delta.coef[d, ]) + as.matrix(xvar[[2]]) %*% as.matrix(delta.coef.level.2[d, ])
                            + as.matrix(xvar[[3]]) %*% as.matrix(delta.coef.level.3[d, ]) )
    if(object$levels==4)
      pred.list <- logistic(as.matrix(xvar[[1]]) %*% as.matrix(delta.coef[d, ]) + as.matrix(xvar[[2]]) %*% as.matrix(delta.coef.level.2[d, ])
                            + as.matrix(xvar[[3]]) %*% as.matrix(delta.coef.level.3[d, ]) + as.matrix(xvar[[4]]) %*% as.matrix(delta.coef.level.4[d, ]))
    
    pred.list.mean[d] <- mean(pred.list)   
    
  }

  critical.value <- qt(1-(1-level)/2, df = nrow(xvar[[1]]))
  
  est.list <- mean(pred.list.mean)
  se.list <- sd(pred.list.mean)
  ci.upper.list <- est.list + critical.value * se.list
  ci.lower.list <- est.list - critical.value * se.list
  
  return.object <- list(fit = est.list)
  
  if ( interval == "confidence") {
    return.object <- as.data.frame(rbind(c(est.list, ci.lower.list, ci.upper.list)))
    names(return.object) <- c("fit","lwr","upr")
  }
  
  if (se.fit == T)
    return.object$se.fit <- c(se.list)
  
  class(return.object) <- "predict.ictreg"
  
  return.object
  
}



ictregBayesMulti2Level.fit <- function(Y, treat, X, V, J, grp,
                                       ##ceiling, floor,
                                       n.draws, burnin, thin,
                                       verbose, delta.start, delta.grp.start, sigma.start, delta.mu0, 
                                       delta.A0, delta.grp.mu0, delta.grp.A0, sigma.df, sigma.scale,
                                       delta.tune, alpha.tune) {

  n <- length(Y)
  k <- ncol(X)
  n.covV <- ncol(V)
  n.grp <- length(table(grp))

  ## treatment variable is 0, 1, 2, ...
  levels.treat <- 1:length(unique(treat))
  tmax <- length(levels.treat)

  ## starting values, prior, and tuning parameters for sensitive item
  ## parameters: either the same starting values and prior for all
  ## sensitive items or a list with the names identical to the levels
  ## of the treatment factor variable

  ## delta includes (fixed effects) parameters for both control and
  ## sensitive item models at the individual level
  if (is.list(delta.start)) {
    delta.start.all <- NULL
    for (i in 1:tmax) {
      delta.start.all <- c(delta.start.all, delta.start[[levels.treat[i]]])
    }
  } else {
    delta.start.all <- rep(delta.start, tmax)
  }

  if (is.list(delta.mu0)) {
    delta.mu0.all <- NULL
    for (i in 1:tmax) {
      delta.mu0.all <- c(delta.mu0.all, delta.mu0[[levels.treat[i]]])
    }
  } else {
    delta.mu0.all <- rep(delta.mu0, tmax)
  }
  
  if (is.list(delta.A0)) {
    delta.A0.all <- NULL
    for (i in 1:tmax) {
      delta.A0.all <- c(delta.A0.all, as.double(delta.A0[[levels.treat[i]]]))
    }
  } else {
    delta.A0.all <- rep(as.double(delta.A0), tmax)
  }

  if (is.list(delta.tune)) {
    delta.tune.all <- NULL
    for (i in 1:tmax) {
      delta.tune.all <- c(delta.tune.all, as.double(delta.tune[[levels.treat[i]]]))
    }
  } else {
    delta.tune.all <- rep(as.double(delta.tune), tmax)
  }

  if (is.list(alpha.tune)) {
    alpha.tune.all <- NULL
    for (i in 1:tmax) {
      alpha.tune.all <- c(alpha.tune.all, as.double(alpha.tune[[levels.treat[i]]]))
    }
  } else {
    alpha.tune.all <- rep(as.double(alpha.tune), tmax)
  }


  ## delta.grp includes (fixed effects) parameters for both control
  ## and sensitive item models at the group level
  if (is.list(delta.grp.start)) {
    delta.grp.start.all <- NULL
    for (i in 1:tmax) {
      delta.grp.start.all <- c(delta.grp.start.all, delta.grp.start[[levels.treat[i]]])
    }
  } else {
    delta.grp.start.all <- rep(delta.grp.start, tmax)
  }

  if (is.list(delta.grp.mu0)) {
    delta.grp.mu0.all <- NULL
    for (i in 1:tmax) {
      delta.grp.mu0.all <- c(delta.grp.mu0.all, delta.grp.mu0[[levels.treat[i]]])
    }
  } else {
    delta.grp.mu0.all <- rep(delta.grp.mu0, tmax)
  }
  
  if (is.list(delta.grp.A0)) {
    delta.grp.A0.all <- NULL
    for (i in 1:tmax) {
      delta.grp.A0.all <- c(delta.grp.A0.all, as.double(delta.grp.A0[[levels.treat[i]]]))
    }
  } else {
    delta.grp.A0.all <- rep(as.double(delta.grp.A0), tmax)
  }

  ## sigma includes the variance parameters for random effects in both
  ## control and sensitive item models
  if (length(sigma.start) == tmax) {
    sigma.start.all <- sigma.start.all <- sigma.start
  } else {
    sigma.start.all <- rep(sigma.start[1], tmax)
  }
 
  if (length(sigma.scale) == tmax) {
    sigma.scale.all <- sigma.scale
  } else {
    sigma.scale.all <- rep(sigma.scale[1], tmax)
  } 

  if (length(sigma.df) == tmax) {
    sigma.df.all <- sigma.df
  } else {
    sigma.df.all <- rep(sigma.df[1], tmax)
  }

  ## fixed effects, sigma, Phi, random effects, acceptance ratios
  keep <- thin + 1
  alldraws <- floor((n.draws - burnin) / keep)
  n.par <- tmax * (k + n.grp + 2 + n.covV)    

  res <- .C("ictregBinomMulti2Level", as.integer(Y), as.integer(J),
            as.integer(n), as.integer(n.draws), as.integer(treat),
            as.integer(tmax-1), as.double(X), as.double(V),
            as.double(delta.start.all), as.double(delta.grp.start.all),
            as.integer(k), as.integer(n.covV), as.double(delta.mu0.all),
            as.double(delta.grp.mu0.all), as.double(delta.A0.all),
            as.double(delta.grp.A0.all), as.double(delta.tune.all),
            ##as.integer(ceiling), as.integer(floor),
            0, 0,
            as.integer(grp-1),
            as.integer(n.grp), as.double(alpha.tune.all),
            as.double(sigma.start.all), as.integer(sigma.df.all),
            as.double(sigma.scale.all), as.integer(burnin),
            as.integer(keep), as.integer(verbose), allresults =
            double(n.par*alldraws), PACKAGE = "list")$allresults

  res <- matrix(res, byrow = TRUE, ncol = n.par)

  return(list(delta = res[, 1:(k*tmax), drop = FALSE],
              alpha = res[, (k*tmax + 1):((k + n.grp) * tmax), drop = FALSE],
              delta.grp = res[, ((k + n.grp) * tmax + 1):((k + n.grp + n.covV) * tmax), drop = FALSE],
              sigma = res[, ((k + n.grp + n.covV) * tmax + 1):((k + n.grp + n.covV + 1) * tmax), drop = FALSE],
              delta.accept = res[alldraws, ((k + n.grp + n.covV + 1) * tmax + 1):n.par]))
}



ictregBayesMulti3Level.fit <- function(Y, treat, X, V, W, J, grp1, grp2, ##ceiling, floor,
                                       n.draws, burnin,
                                       thin, verbose, delta.start, delta.grp1.start, delta.grp2.start,
                                       sigma.grp1.start, sigma.grp2.start, delta.mu0, delta.A0,
                                       delta.grp1.mu0, delta.grp1.A0, delta.grp2.mu0, delta.grp2.A0,
                                       sigma.grp1.df, sigma.grp1.scale, sigma.grp2.df, sigma.grp2.scale,
                                       delta.tune, alpha.tune) {

  n <- length(Y)
  k <- ncol(X)
  n.covV <- ncol(V)
  n.covW <- ncol(W)
  n.grp1 <- length(table(grp1))
  n.grp2 <- length(table(grp2))

  ## treatment variable is 0, 1, 2, ...
  levels.treat <- 1:length(unique(treat))
  tmax <- length(levels.treat)

  ## starting values, prior, and tuning parameters for sensitive item
  ## parameters: either the same starting values and prior for all
  ## sensitive items or a list with the names identical to the levels
  ## of the treatment factor variable

  ## delta includes (fixed effects) parameters for both control and
  ## sensitive item models at the individual level
  if (is.list(delta.start)) {
    delta.start.all <- NULL
    for (i in 1:tmax) {
      delta.start.all <- c(delta.start.all, delta.start[[levels.treat[i]]])
    }
  } else {
    delta.start.all <- rep(delta.start, tmax)
  }

  if (is.list(delta.mu0)) {
    delta.mu0.all <- NULL
    for (i in 1:tmax) {
      delta.mu0.all <- c(delta.mu0.all, delta.mu0[[levels.treat[i]]])
    }
  } else {
    delta.mu0.all <- rep(delta.mu0, tmax)
  }
  
  if (is.list(delta.A0)) {
    delta.A0.all <- NULL
    for (i in 1:tmax) {
      delta.A0.all <- c(delta.A0.all, as.double(delta.A0[[levels.treat[i]]]))
    }
  } else {
    delta.A0.all <- rep(as.double(delta.A0), tmax)
  }

  if (is.list(delta.tune)) {
    delta.tune.all <- NULL
    for (i in 1:tmax) {
      delta.tune.all <- c(delta.tune.all, as.double(delta.tune[[levels.treat[i]]]))
    }
  } else {
    delta.tune.all <- rep(as.double(delta.tune), tmax)
  }

  if (is.list(alpha.tune)) {
    alpha.tune.all <- NULL
    for (i in 1:tmax) {
      alpha.tune.all <- c(alpha.tune.all, as.double(alpha.tune[[levels.treat[i]]]))
    }
  } else {
    alpha.tune.all <- rep(as.double(alpha.tune), tmax)
  }

  ## delta.grp1 includes (fixed effects) parameters for both control
  ## and sensitive item models at the group 1 level
  if (is.list(delta.grp1.start)) {
    delta.grp1.start.all <- NULL
    for (i in 1:tmax) {
      delta.grp1.start.all <- c(delta.grp1.start.all, delta.grp1.start[[levels.treat[i]]])
    }
  } else {
    delta.grp1.start.all <- rep(delta.grp1.start, tmax)
  }

  if (is.list(delta.grp1.mu0)) {
    delta.grp1.mu0.all <- NULL
    for (i in 1:tmax) {
      delta.grp1.mu0.all <- c(delta.grp1.mu0.all, delta.grp1.mu0[[levels.treat[i]]])
    }
  } else {
    delta.grp1.mu0.all <- rep(delta.grp1.mu0, tmax)
  }
  
  if (is.list(delta.grp1.A0)) {
    delta.grp1.A0.all <- NULL
    for (i in 1:tmax) {
      delta.grp1.A0.all <- c(delta.grp1.A0.all, as.double(delta.grp1.A0[[levels.treat[i]]]))
    }
  } else {
    delta.grp1.A0.all <- rep(as.double(delta.grp1.A0), tmax)
  }

  ## delta.grp2 includes (fixed effects) parameters for both control
  ## and sensitive item models at the group 2 level
  if (is.list(delta.grp2.start)) {
    delta.grp2.start.all <- NULL
    for (i in 1:tmax) {
      delta.grp2.start.all <- c(delta.grp2.start.all, delta.grp2.start[[levels.treat[i]]])
    }
  } else {
    delta.grp2.start.all <- rep(delta.grp2.start, tmax)
  }

  if (is.list(delta.grp2.mu0)) {
    delta.grp2.mu0.all <- NULL
    for (i in 1:tmax) {
      delta.grp2.mu0.all <- c(delta.grp2.mu0.all, delta.grp2.mu0[[levels.treat[i]]])
    }
  } else {
    delta.grp2.mu0.all <- rep(delta.grp2.mu0, tmax)
  }
  
  if (is.list(delta.grp2.A0)) {
    delta.grp2.A0.all <- NULL
    for (i in 1:tmax) {
      delta.grp2.A0.all <- c(delta.grp2.A0.all, as.double(delta.grp2.A0[[levels.treat[i]]]))
    }
  } else {
    delta.grp2.A0.all <- rep(as.double(delta.grp2.A0), tmax)
  }

  ## sigma.grp1 includes the variance parameters for random effects in both
  ## control and sensitive item models at group 1 level
  if (length(sigma.grp1.start) == tmax) {
    sigma.grp1.start.all <- sigma.grp1.start.all <- sigma.grp1.start
  } else {
    sigma.grp1.start.all <- rep(sigma.grp1.start[1], tmax)
  }
 
  if (length(sigma.grp1.scale) == tmax) {
    sigma.grp1.scale.all <- sigma.grp1.scale
  } else {
    sigma.grp1.scale.all <- rep(sigma.grp1.scale[1], tmax)
  } 

  if (length(sigma.grp1.df) == tmax) {
    sigma.grp1.df.all <- sigma.grp1.df
  } else {
    sigma.grp1.df.all <- rep(sigma.grp1.df[1], tmax)
  }
  
  ## sigma.grp2 includes the variance parameters for random effects in both
  ## control and sensitive item models
  if (length(sigma.grp2.start) == tmax) {
    sigma.grp2.start.all <- sigma.grp2.start.all <- sigma.grp2.start
  } else {
    sigma.grp2.start.all <- rep(sigma.grp2.start[1], tmax)
  }
 
  if (length(sigma.grp2.scale) == tmax) {
    sigma.grp2.scale.all <- sigma.grp2.scale
  } else {
    sigma.grp2.scale.all <- rep(sigma.grp2.scale[1], tmax)
  } 

  if (length(sigma.grp2.df) == tmax) {
    sigma.grp2.df.all <- sigma.grp2.df
  } else {
    sigma.grp2.df.all <- rep(sigma.grp2.df[1], tmax)
  }

  ## fixed effects, 2 sigmas, random effects, acceptance ratios
  keep <- thin + 1
  alldraws <- floor((n.draws - burnin) / keep)
  n.par <- tmax * (k + n.grp1 + n.grp2 + 3 + n.covV + n.covW)    

  res <- .C("ictregBinomMulti3Level", as.integer(Y), as.integer(J),
            as.integer(n), as.integer(n.draws), as.integer(treat),
            as.integer(tmax-1), as.double(X), as.double(V), as.double(W),
            as.double(delta.start.all), as.double(delta.grp1.start.all),
            as.double(delta.grp2.start.all), as.integer(k),
            as.integer(n.covV), as.integer(n.covW), as.double(delta.mu0.all),
            as.double(delta.grp1.mu0.all), as.double(delta.grp2.mu0.all),
            as.double(delta.A0.all), as.double(delta.grp1.A0.all),
            as.double(delta.grp2.A0.all), as.double(delta.tune.all),
           ## as.integer(ceiling), as.integer(floor),
            0, 0, as.integer(grp1-1),
            as.integer(n.grp1), as.integer(grp2-1), as.integer(n.grp2),
            as.double(alpha.tune.all), as.double(sigma.grp1.start.all),
            as.double(sigma.grp2.start.all), as.integer(sigma.grp1.df.all),
            as.double(sigma.grp1.scale.all), as.integer(sigma.grp2.df.all),
            as.double(sigma.grp2.scale.all), as.integer(burnin),
            as.integer(keep), as.integer(verbose), allresults =
            double(n.par*alldraws), PACKAGE = "list")$allresults

  res <- matrix(res, byrow = TRUE, ncol = n.par)

  return(list(delta = res[, 1:(k*tmax), drop = FALSE],
              alpha.grp1 = res[, (k*tmax + 1):((k + n.grp1) * tmax), drop = FALSE],
              delta.grp1 = res[, ((k + n.grp1) * tmax + 1):((k + n.grp1 + n.covV) * tmax), drop = FALSE],
              alpha.grp2 = res[, ((k + n.grp1 + n.covV) * tmax + 1):((k + n.grp1 + n.covV + n.grp2) * tmax), drop = FALSE],
              delta.grp2 = res[, ((k + n.grp1 + n.covV + n.grp2) * tmax + 1):((k + n.grp1 + n.covV + n.grp2 + n.covW) * tmax), drop = FALSE],
              sigma.grp1 = res[, ((k + n.grp1 + n.covV + n.grp2 + n.covW) * tmax + 1):((k + n.grp1 + n.covV + n.grp2 + n.covW + 1) * tmax), drop = FALSE],
              sigma.grp2 = res[, ((k + n.grp1 + n.covV + n.grp2 + n.covW + 1) * tmax + 1):((k + n.grp1 + n.covV + n.grp2 + n.covW + 2) * tmax), drop = FALSE],
              delta.accept = res[alldraws, ((k + n.grp1 + n.covV + n.grp2 + n.covW + 2) * tmax + 1):n.par]))
}



ictregBayesMulti4Level.fit <- function(Y, treat, X, V, W, Z, J, grp1, grp2, grp3,
                                       ##ceiling, floor,
                                       n.draws, burnin, thin, verbose, delta.start, delta.grp1.start,
                                       delta.grp2.start, delta.grp3.start, sigma.grp1.start,
                                       sigma.grp2.start, sigma.grp3.start, delta.mu0, delta.A0,
                                       delta.grp1.mu0, delta.grp1.A0, delta.grp2.mu0, delta.grp2.A0,
                                       delta.grp3.mu0, delta.grp3.A0, sigma.grp1.df, sigma.grp1.scale,
                                       sigma.grp2.df, sigma.grp2.scale, sigma.grp3.df, sigma.grp3.scale,
                                       delta.tune, alpha.tune) {

  n <- length(Y)
  k <- ncol(X)
  n.covV <- ncol(V)
  n.covW <- ncol(W)
  n.covZ <- ncol(Z)
  n.grp1 <- length(table(grp1))
  n.grp2 <- length(table(grp2))
  n.grp3 <- length(table(grp3))

  ## treatment variable is 0, 1, 2, ...
  levels.treat <- 1:length(unique(treat))
  tmax <- length(levels.treat)

  ## starting values, prior, and tuning parameters for sensitive item
  ## parameters: either the same starting values and prior for all
  ## sensitive items or a list with the names identical to the levels
  ## of the treatment factor variable

  ## delta includes (fixed effects) parameters for both control and
  ## sensitive item models at the individual level
  if (is.list(delta.start)) {
    delta.start.all <- NULL
    for (i in 1:tmax) {
      delta.start.all <- c(delta.start.all, delta.start[[levels.treat[i]]])
    }
  } else {
    delta.start.all <- rep(delta.start, tmax)
  }

  if (is.list(delta.mu0)) {
    delta.mu0.all <- NULL
    for (i in 1:tmax) {
      delta.mu0.all <- c(delta.mu0.all, delta.mu0[[levels.treat[i]]])
    }
  } else {
    delta.mu0.all <- rep(delta.mu0, tmax)
  }
  
  if (is.list(delta.A0)) {
    delta.A0.all <- NULL
    for (i in 1:tmax) {
      delta.A0.all <- c(delta.A0.all, as.double(delta.A0[[levels.treat[i]]]))
    }
  } else {
    delta.A0.all <- rep(as.double(delta.A0), tmax)
  }

  if (is.list(delta.tune)) {
    delta.tune.all <- NULL
    for (i in 1:tmax) {
      delta.tune.all <- c(delta.tune.all, as.double(delta.tune[[levels.treat[i]]]))
    }
  } else {
    delta.tune.all <- rep(as.double(delta.tune), tmax)
  }

  if (is.list(alpha.tune)) {
    alpha.tune.all <- NULL
    for (i in 1:tmax) {
      alpha.tune.all <- c(alpha.tune.all, as.double(alpha.tune[[levels.treat[i]]]))
    }
  } else {
    alpha.tune.all <- rep(as.double(alpha.tune), tmax)
  }


  ## delta.grp1 includes (fixed effects) parameters for both control
  ## and sensitive item models at the group 1 level
  if (is.list(delta.grp1.start)) {
    delta.grp1.start.all <- NULL
    for (i in 1:tmax) {
      delta.grp1.start.all <- c(delta.grp1.start.all, delta.grp1.start[[levels.treat[i]]])
    }
  } else {
    delta.grp1.start.all <- rep(delta.grp1.start, tmax)
  }

  if (is.list(delta.grp1.mu0)) {
    delta.grp1.mu0.all <- NULL
    for (i in 1:tmax) {
      delta.grp1.mu0.all <- c(delta.grp1.mu0.all, delta.grp1.mu0[[levels.treat[i]]])
    }
  } else {
    delta.grp1.mu0.all <- rep(delta.grp1.mu0, tmax)
  }
  
  if (is.list(delta.grp1.A0)) {
    delta.grp1.A0.all <- NULL
    for (i in 1:tmax) {
      delta.grp1.A0.all <- c(delta.grp1.A0.all, as.double(delta.grp1.A0[[levels.treat[i]]]))
    }
  } else {
    delta.grp1.A0.all <- rep(as.double(delta.grp1.A0), tmax)
  }

  ## delta.grp2 includes (fixed effects) parameters for both control
  ## and sensitive item models at the group 2 level
  if (is.list(delta.grp2.start)) {
    delta.grp2.start.all <- NULL
    for (i in 1:tmax) {
      delta.grp2.start.all <- c(delta.grp2.start.all, delta.grp2.start[[levels.treat[i]]])
    }
  } else {
    delta.grp2.start.all <- rep(delta.grp2.start, tmax)
  }

  if (is.list(delta.grp2.mu0)) {
    delta.grp2.mu0.all <- NULL
    for (i in 1:tmax) {
      delta.grp2.mu0.all <- c(delta.grp2.mu0.all, delta.grp2.mu0[[levels.treat[i]]])
    }
  } else {
    delta.grp2.mu0.all <- rep(delta.grp2.mu0, tmax)
  }
  
  if (is.list(delta.grp2.A0)) {
    delta.grp2.A0.all <- NULL
    for (i in 1:tmax) {
      delta.grp2.A0.all <- c(delta.grp2.A0.all, as.double(delta.grp2.A0[[levels.treat[i]]]))
    }
  } else {
    delta.grp2.A0.all <- rep(as.double(delta.grp2.A0), tmax)
  }

  ## delta.grp3 includes (fixed effects) parameters for both control
  ## and sensitive item models at the group 3 level
  if (is.list(delta.grp3.start)) {
    delta.grp3.start.all <- NULL
    for (i in 1:tmax) {
      delta.grp3.start.all <- c(delta.grp3.start.all, delta.grp3.start[[levels.treat[i]]])
    }
  } else {
    delta.grp3.start.all <- rep(delta.grp3.start, tmax)
  }

  if (is.list(delta.grp3.mu0)) {
    delta.grp3.mu0.all <- NULL
    for (i in 1:tmax) {
      delta.grp3.mu0.all <- c(delta.grp3.mu0.all, delta.grp3.mu0[[levels.treat[i]]])
    }
  } else {
    delta.grp3.mu0.all <- rep(delta.grp3.mu0, tmax)
  }
  
  if (is.list(delta.grp3.A0)) {
    delta.grp3.A0.all <- NULL
    for (i in 1:tmax) {
      delta.grp3.A0.all <- c(delta.grp3.A0.all, as.double(delta.grp3.A0[[levels.treat[i]]]))
    }
  } else {
    delta.grp3.A0.all <- rep(as.double(delta.grp3.A0), tmax)
  }
  
  ## sigma.grp1 includes the variance parameters for random effects in both
  ## control and sensitive item models at group 1 level
  if (length(sigma.grp1.start) == tmax) {
    sigma.grp1.start.all <- sigma.grp1.start.all <- sigma.grp1.start
  } else {
    sigma.grp1.start.all <- rep(sigma.grp1.start[1], tmax)
  }
 
  if (length(sigma.grp1.scale) == tmax) {
    sigma.grp1.scale.all <- sigma.grp1.scale
  } else {
    sigma.grp1.scale.all <- rep(sigma.grp1.scale[1], tmax)
  } 

  if (length(sigma.grp1.df) == tmax) {
    sigma.grp1.df.all <- sigma.grp1.df
  } else {
    sigma.grp1.df.all <- rep(sigma.grp1.df[1], tmax)
  }
  
  ## sigma.grp2 includes the variance parameters for random effects in both
  ## control and sensitive item models
  if (length(sigma.grp2.start) == tmax) {
    sigma.grp2.start.all <- sigma.grp2.start.all <- sigma.grp2.start
  } else {
    sigma.grp2.start.all <- rep(sigma.grp2.start[1], tmax)
  }
 
  if (length(sigma.grp2.scale) == tmax) {
    sigma.grp2.scale.all <- sigma.grp2.scale
  } else {
    sigma.grp2.scale.all <- rep(sigma.grp2.scale[1], tmax)
  } 

  if (length(sigma.grp2.df) == tmax) {
    sigma.grp2.df.all <- sigma.grp2.df
  } else {
    sigma.grp2.df.all <- rep(sigma.grp2.df[1], tmax)
  }

  ## sigma.grp3 includes the variance parameters for random effects in both
  ## control and sensitive item models
  if (length(sigma.grp3.start) == tmax) {
    sigma.grp3.start.all <- sigma.grp3.start.all <- sigma.grp3.start
  } else {
    sigma.grp3.start.all <- rep(sigma.grp3.start[1], tmax)
  }
 
  if (length(sigma.grp3.scale) == tmax) {
    sigma.grp3.scale.all <- sigma.grp3.scale
  } else {
    sigma.grp3.scale.all <- rep(sigma.grp3.scale[1], tmax)
  } 

  if (length(sigma.grp3.df) == tmax) {
    sigma.grp3.df.all <- sigma.grp3.df
  } else {
    sigma.grp3.df.all <- rep(sigma.grp3.df[1], tmax)
  }

  
  ## fixed effects, 2 sigmas, random effects, acceptance ratios
  keep <- thin + 1
  alldraws <- floor((n.draws - burnin) / keep)
  n.par <- tmax * (k + n.grp1 + n.grp2 + n.grp3 + 4 + n.covV + n.covW + n.covZ)    

  res <- .C("ictregBinomMulti4Level", as.integer(Y), as.integer(J),
            as.integer(n), as.integer(n.draws), as.integer(treat),
            as.integer(tmax-1), as.double(X), as.double(V), as.double(W),
            as.double(Z), as.double(delta.start.all), as.double(delta.grp1.start.all),
            as.double(delta.grp2.start.all), as.double(delta.grp3.start.all),
            as.integer(k), as.integer(n.covV), as.integer(n.covW), as.integer(n.covZ),
            as.double(delta.mu0.all), as.double(delta.grp1.mu0.all),
            as.double(delta.grp2.mu0.all), as.double(delta.grp3.mu0.all),
            as.double(delta.A0.all), as.double(delta.grp1.A0.all),
            as.double(delta.grp2.A0.all), as.double(delta.grp3.A0.all),
            as.double(delta.tune.all), ##as.integer(ceiling), as.integer(floor),
            0, 0,
            as.integer(grp1-1), as.integer(n.grp1), as.integer(grp2-1), as.integer(n.grp2),
            as.integer(grp3-1), as.integer(n.grp3), as.double(alpha.tune.all),
            as.double(sigma.grp1.start.all), as.double(sigma.grp2.start.all),
            as.double(sigma.grp3.start.all), as.integer(sigma.grp1.df.all),
            as.double(sigma.grp1.scale.all), as.integer(sigma.grp2.df.all),
            as.double(sigma.grp2.scale.all), as.integer(sigma.grp3.df.all),
            as.double(sigma.grp3.scale.all),as.integer(burnin),
            as.integer(keep), as.integer(verbose), allresults =
            double(n.par*alldraws), PACKAGE = "list")$allresults

  res <- matrix(res, byrow = TRUE, ncol = n.par)

  return(list(delta = res[, 1:(k*tmax), drop = FALSE],
              alpha.grp1 = res[, (k*tmax + 1):((k + n.grp1) * tmax), drop = FALSE],
              delta.grp1 = res[, ((k + n.grp1) * tmax + 1):((k + n.grp1 + n.covV) * tmax), drop = FALSE],
              alpha.grp2 = res[, ((k + n.grp1 + n.covV) * tmax + 1):((k + n.grp1 + n.covV + n.grp2) * tmax), drop = FALSE],
              delta.grp2 = res[, ((k + n.grp1 + n.covV + n.grp2) * tmax + 1):((k + n.grp1 + n.covV + n.grp2 + n.covW) * tmax), drop = FALSE],
              alpha.grp3 = res[, ((k + n.grp1 + n.covV + n.grp2 + n.covW) * tmax + 1):((k + n.grp1 + n.covV + n.grp2+ n.covW + n.grp3) * tmax), drop = FALSE],
              delta.grp3 = res[, ((k + n.grp1 + n.covV + n.grp2 + n.covW + n.grp3) * tmax + 1):((k + n.grp1 + n.covV + n.grp2 + n.covW + n.grp3 + n.covZ) * tmax), drop = FALSE],
              sigma.grp1 = res[, ((k + n.grp1 + n.covV + n.grp2 + n.covW + n.grp3 + n.covZ) * tmax + 1):((k + n.grp1 + n.covV + n.grp2 + n.covW + n.grp3 + n.covZ + 1) * tmax), drop = FALSE],
              sigma.grp2 = res[, ((k + n.grp1 + n.covV + n.grp2 + n.covW + n.grp3 + n.covZ + 1) * tmax + 1):((k + n.grp1 + n.covV + n.grp2 + n.covW + n.grp3 + n.covZ + 2) * tmax), drop = FALSE],
              sigma.grp3 = res[, ((k + n.grp1 + n.covV + n.grp2 + n.covW + n.grp3 + n.covZ + 2) * tmax + 1):((k + n.grp1 + n.covV + n.grp2 + n.covW + n.grp3 + n.covZ + 3) * tmax), drop = FALSE],
              delta.accept = res[alldraws, ((k + n.grp1 + n.covV + n.grp2 + n.covW + n.grp3 + n.covZ + 3) * tmax + 1):n.par]))
}
