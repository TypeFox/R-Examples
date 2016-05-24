# Effect generic and methods
# John Fox and Sanford Weisberg
# 12-21-2012 Allow for empty cells in factor interactions, S. Weisberg
# 2012-03-05: Added .merMod method for development version of lme4, J. Fox
# 2012-04-06: Added support for lme4.0, J. Fox
# 2013-07-15:  Changed default xlevels and default.levels
# 2013-10-15: Added Effect.default(). J. Fox
# 2013-10-22: fixed bug in Effect.lm() when na.action=na.exclude. J. Fox
# 2013-10-29: code to handle "valid" NAs in factors. J. Fox
# 2013-11-06: fixed bug in Effect.multinom() in construction of effect object
#             when there is only one focal predictor; caused as.data.frame.effpoly() to fail
# 2014-03-13: modified Effect.lm() to compute partial residuals. J. Fox
# 2014-05-06: fixed bug in Effect.gls() when cor or var structure depends on variables in the data set. J. Fox
# 2014-08-02: added vcov.=vcov argument to allow other methods of estimating var(coef.estimates)
# 2014-09-25: added KR argument to Effect.mer() and Effect.merMod(). J. Fox
# 2014-12-07: don't assume that pbkrtest is installed. J. Fox
# 2015-03-25: added "family" element to eff objects returned by Effect.lm(). J. Fox
# 2016-02-16: fixed problem in handling terms like polynomials for non-focal predictors. J. Fox
# 2016-03-01: recoded calculation of partial residuals. J. Fox

Effect <- function(focal.predictors, mod, ...){
  UseMethod("Effect", mod)
}

Effect.lm <- function (focal.predictors, mod, xlevels = list(), 
                       default.levels = NULL, given.values,
                       vcov. = vcov, se = TRUE, confidence.level = 0.95, 
                       transformation = list(link = family(mod)$linkfun, inverse = family(mod)$linkinv), 
                       typical = mean, offset = mean, partial.residuals=FALSE, quantiles=seq(0.2, 0.8, by=0.2),
                       x.var=NULL,  ...){
  data <- if (partial.residuals){
    all.vars <- all.vars(formula(mod))
    expand.model.frame(mod, all.vars)[, all.vars]
  }
  else NULL
  if (missing(given.values)) 
    given.values <- NULL
  else if (!all(which <- names(given.values) %in% names(coef(mod)))) 
    stop("given.values (", names(given.values[!which]), ") not in the model")
  off <- if (is.numeric(offset) && length(offset) == 1) offset
  else if (is.function(offset)) {
    mod.off <- model.offset(model.frame(mod))
    if (is.null(mod.off)) 0 else offset(mod.off)
  }
  else stop("offset must be a function or a number")
  formula.rhs <- formula(mod)[[3]]
  model.components <- Analyze.model(focal.predictors, mod, xlevels, default.levels, formula.rhs, 
                                    partial.residuals=partial.residuals, quantiles=quantiles, x.var=x.var, data=data, typical=typical)
  excluded.predictors <- model.components$excluded.predictors
  predict.data <- model.components$predict.data
  predict.data.all.rounded <- predict.data.all <- if (partial.residuals) na.omit(data[, all.vars(formula(mod))]) else NULL
  factor.levels <- model.components$factor.levels
  factor.cols <- model.components$factor.cols
  n.focal <- model.components$n.focal
  x <- model.components$x
  X.mod <- model.components$X.mod
  cnames <- model.components$cnames
  X <- model.components$X
  x.var <- model.components$x.var
  formula.rhs <- formula(mod)[c(1, 3)]
  Terms <- delete.response(terms(mod))
  mf <- model.frame(Terms, predict.data, xlev = factor.levels, na.action=NULL)
  mod.matrix <- model.matrix(formula.rhs, data = mf, contrasts.arg = mod$contrasts)
  if (is.null(x.var)) partial.residuals <- FALSE
  factors <- sapply(predict.data, is.factor)
  if (partial.residuals){
    for (predictor in focal.predictors[-x.var]){
      if (!factors[predictor]){
        values <- unique(predict.data[, predictor])
        predict.data.all.rounded[, predictor] <- values[apply(outer(predict.data.all[, predictor], values, function(x, y) (x - y)^2), 1, which.min)]
      }
    }
  }
  mod.matrix.all <- model.matrix(mod)
  wts <- weights(mod) 
  if (is.null(wts)) 
    wts <- rep(1, length(residuals(mod)))
  mod.matrix <- Fixup.model.matrix(mod, mod.matrix, mod.matrix.all, 
                                   X.mod, factor.cols, cnames, focal.predictors, excluded.predictors, 
                                   typical, given.values) #,
  # look for aliased coefficients and remove those columns from mod.matrix
  mod.matrix <- mod.matrix[, !is.na(mod$coefficients)]
  effect <- off + mod.matrix %*% mod$coefficients[!is.na(mod$coefficients)]
  if (partial.residuals){
    res <- na.omit(residuals(mod, type="working"))
    fitted <- na.omit(if (inherits(mod, "glm")) predict(mod, type="link") else predict(mod))
    partial.residuals.range <- range(fitted + res)
  }
  else {
    res <- partial.residuals.range <- NULL
  }
  result <- list(term = paste(focal.predictors, collapse="*"), 
                 formula = formula(mod), response = response.name(mod), 
                 variables = x, fit = effect, x = predict.data[, 1:n.focal, drop=FALSE], 
                 x.all=predict.data.all.rounded[, focal.predictors, drop=FALSE],
                 model.matrix = mod.matrix, 
                 data = X, 
                 discrepancy = 0, offset=off,
                 residuals=res, partial.residuals.range=partial.residuals.range,
                 x.var=x.var)
  # find empty cells, if any, and correct
  whichFact <- unlist(lapply(result$variables, function(x) x$is.factor))
  zeroes <- NULL
  if(sum(whichFact) > 1){
    nameFact <- names(whichFact)[whichFact]
    counts <- xtabs(as.formula( paste("~", paste(nameFact, collapse="+"))), 
                    model.frame(mod))
    zeroes <- which(counts == 0)  
  }
  if(length(zeroes) > 0){
    levs <- expand.grid(lapply(result$variables, function(x) x$levels)) 
    good <- rep(TRUE, dim(levs)[1])
    for(z in zeroes){
      good <- good &  
        apply(levs, 1, function(x) !all(x == levs[z, whichFact]))
    } 
    result$fit[!good] <- NA
  } 
  if (se) {
    if (any(family(mod)$family == c("binomial", "poisson"))) {
      z <- qnorm(1 - (1 - confidence.level)/2)
    }
    else {
      z <- qt(1 - (1 - confidence.level)/2, df = mod$df.residual)
    }
    V <- vcov.(mod)
    eff.vcov <- mod.matrix %*% V %*% t(mod.matrix)
    rownames(eff.vcov) <- colnames(eff.vcov) <- NULL
    var <- diag(eff.vcov)
    result$vcov <- eff.vcov        
    result$se <- sqrt(var)
    result$lower <- effect - z * result$se
    result$upper <- effect + z * result$se
    result$confidence.level <- confidence.level
    if(length(zeroes) > 0){
      result$se[!good] <- NA
      result$lower[!good] <- NA
      result$upper[!good] <- NA
    }
  }
  if (is.null(transformation$link) && is.null(transformation$inverse)) {
    transformation$link <- I
    transformation$inverse <- I
  }
  result$transformation <- transformation
  result$family <- family(mod)$family
  class(result) <- "eff"
  result
}

Effect.mer <- function(focal.predictors, mod, KR=FALSE, ...) {
  result <- Effect(focal.predictors, mer.to.glm(mod, KR=KR), ...)
  result$formula <- as.formula(formula(mod))
  result
}

Effect.merMod <- function(focal.predictors, mod, KR=FALSE, ...){
  Effect.mer(focal.predictors, mod, KR=KR, ...)
}

Effect.lme <- function(focal.predictors, mod, ...) {
  result <- Effect(focal.predictors, lme.to.glm(mod), ...)
  result$formula <- as.formula(formula(mod))
  result
}

Effect.gls <- function (focal.predictors, mod, xlevels = list(), default.levels = NULL, given.values,
                        vcov. = vcov, se = TRUE, confidence.level = 0.95, 
                        transformation = NULL, 
                        typical = mean, ...){
  if (missing(given.values)) 
    given.values <- NULL
  else if (!all(which <- names(given.values) %in% names(coef(mod)))) 
    stop("given.values (", names(given.values[!which]), ") not in the model")
  formula.rhs <- formula(mod)[[3]]
  .data <- eval(mod$call$data)
  mod.lm <- lm(as.formula(mod$call$model), data=.data, na.action=na.exclude)
  model.components <- Analyze.model(focal.predictors, mod.lm, xlevels, default.levels, formula.rhs, typical=typical)
  excluded.predictors <- model.components$excluded.predictors
  predict.data <- model.components$predict.data
  factor.levels <- model.components$factor.levels
  factor.cols <- model.components$factor.cols
  n.focal <- model.components$n.focal
  x <- model.components$x
  X.mod <- model.components$X.mod
  cnames <- model.components$cnames
  X <- model.components$X
  formula.rhs <- formula(mod)[c(1, 3)]
  nrow.X <- nrow(X)
  mf <- model.frame(formula.rhs, data=rbind(X[,names(predict.data),drop=FALSE], predict.data), 
                    xlev=factor.levels)
  mod.matrix.all <- model.matrix(formula.rhs, data=mf, contrasts.arg=mod$contrasts)
  mod.matrix <- mod.matrix.all[-(1:nrow.X),]
  mod.matrix <- Fixup.model.matrix(mod.lm, mod.matrix, model.matrix(mod.lm), 
                                   X.mod, factor.cols, cnames, focal.predictors, excluded.predictors, typical, given.values)
  fit.1 <- na.omit(predict(mod))
  mod.2 <- lm.fit(mod.matrix.all[1:nrow.X,], fit.1)
  class(mod.2) <- "lm"
  use <- !is.na(residuals(mod.lm))
  .data <- .data[use, ]
  .data$.y <- model.response.gls(mod)
  .data$.X <- mod.matrix.all[1:nrow.X, ]
  mod.3 <- update(mod, .y ~ .X - 1, data=.data)
  discrepancy <- 100*mean(abs(fitted(mod.2)- fit.1)/(1e-10 + mean(abs(fit.1))))
  if (discrepancy > 1e-3) warning(paste("There is a discrepancy of", round(discrepancy, 3),
                                        "percent \n     in the 'safe' predictions used to generate effect", paste(focal.predictors, collapse="*")))
  effect <- mod.matrix %*% mod$coefficients
  result <- list(term = paste(focal.predictors, collapse="*"), 
                 formula = formula(mod), response = response.name(mod), 
                 variables = x, fit = effect, x = predict.data[, 1:n.focal, drop=FALSE], model.matrix = mod.matrix, data = X, 
                 discrepancy = discrepancy, offset=0)
  if (se){
    df.res <- mod$dims[["N"]] - mod$dims[["p"]]
    z <- qt(1 - (1 - confidence.level)/2, df=df.res)
    mod.2$terms <- terms(mod)
    V <- vcov.(mod.3)
    eff.vcov <- mod.matrix %*% V %*% t(mod.matrix)
    rownames(eff.vcov) <- colnames(eff.vcov) <- NULL
    var <- diag(eff.vcov)
    result$vcov <- eff.vcov
    result$se <- sqrt(var)        
    result$lower <- effect - z*result$se
    result$upper <- effect + z*result$se
    result$confidence.level <- confidence.level
  }
  if (is.null(transformation$link) && is.null(transformation$inverse)){
    transformation$link <- I
    transformation$inverse <- I
  }
  result$transformation <- transformation
  class(result) <- "eff"
  result
}

Effect.multinom <- function(focal.predictors, mod, 
                            confidence.level=.95, xlevels=list(), default.levels=NULL,
                            given.values, vcov. = vcov, se=TRUE, typical=mean, ...){    
  if (length(mod$lev) < 3) stop("effects for multinomial logit model only available for response levels > 2")
  if (missing(given.values)) given.values <- NULL
  else if (!all(which <- colnames(given.values) %in% names(coef(mod)))) 
    stop("given.values (", colnames(given.values[!which]),") not in the model")
  formula.rhs <- formula(mod)[c(1, 3)]
  model.components <- Analyze.model(focal.predictors, mod, xlevels, default.levels, formula.rhs, typical=typical)
  excluded.predictors <- model.components$excluded.predictors
  predict.data <- model.components$predict.data
  factor.levels <- model.components$factor.levels
  factor.cols <- model.components$factor.cols
  #    n.focal <- model.components$n.focal
  x <- model.components$x
  X.mod <- model.components$X.mod
  cnames <- model.components$cnames
  X <- model.components$X
  formula.rhs <- formula(mod)[c(1, 3)]
  Terms <- delete.response(terms(mod))
  mf <- model.frame(Terms, predict.data, xlev = factor.levels)
  mod.matrix <- model.matrix(formula.rhs, data = mf, contrasts.arg = mod$contrasts)
  X0 <- Fixup.model.matrix(mod, mod.matrix, model.matrix(mod), 
                           X.mod, factor.cols, cnames, focal.predictors, excluded.predictors, typical, given.values)
  resp.names <- make.names(mod$lev, unique=TRUE)
  resp.names <- c(resp.names[-1], resp.names[1]) # make the last level the reference level
  B <- t(coef(mod))
  V <- vcov.(mod)
  m <- ncol(B) + 1
  p <- nrow(B)
  r <- p*(m - 1)	
  n <- nrow(X0)
  P <- Logit <- matrix(0, n, m)
  colnames(P) <-  paste("prob.", resp.names, sep="")
  colnames(Logit) <-  paste("logit.", resp.names, sep="")
  if (se){
    z <- qnorm(1 - (1 - confidence.level)/2)
    Lower.P <- Upper.P <- Lower.logit <- Upper.logit <- SE.P <- SE.logit <- matrix(0, n, m)
    colnames(Lower.logit) <-  paste("L.logit.", resp.names, sep="")
    colnames(Upper.logit) <-  paste("U.logit.", resp.names, sep="")
    colnames(Lower.P) <-  paste("L.prob.", resp.names, sep="")
    colnames(Upper.P) <-  paste("U.prob.", resp.names, sep="")
    colnames(SE.P) <-  paste("se.prob.", resp.names, sep="")
    colnames(SE.logit) <-  paste("se.logit.", resp.names, sep="")
  }
  for (i in 1:n){
    res <- eff.mul(X0[i,], B, se, m, p, r, V) # compute effects
    #        P[i,] <- prob <- res$p # fitted probabilities
    P[i,] <- res$p # fitted probabilities
    Logit[i,] <- logit <- res$logits # fitted logits
    if (se){
      #            SE.P[i,] <- se.p <- res$std.err.p # std. errors of fitted probs
      SE.P[i,] <- res$std.err.p # std. errors of fitted probs	
      SE.logit[i,] <- se.logit <- res$std.error.logits # std. errors of logits
      Lower.P[i,] <- logit2p(logit - z*se.logit)
      Upper.P[i,] <- logit2p(logit + z*se.logit)
      Lower.logit[i,] <- logit - z*se.logit
      Upper.logit[i,] <- logit + z*se.logit
    }
  }
  resp.levs <- c(m, 1:(m-1)) # restore the order of the levels
  P <- P[, resp.levs]
  Logit <- Logit[, resp.levs]
  if (se){
    Lower.P <- Lower.P[, resp.levs]
    Upper.P <- Upper.P[, resp.levs]
    Lower.logit <- Lower.logit[, resp.levs]
    Upper.logit <- Upper.logit[, resp.levs]
    SE.P <- SE.P[, resp.levs]
    SE.logit <- SE.logit[, resp.levs]
  }
  result <- list(term=paste(focal.predictors, collapse="*"), formula=formula(mod), response=response.name(mod),
                 y.levels=mod$lev, variables=x, x=predict.data[, focal.predictors, drop=FALSE],
                 model.matrix=X0, data=X, discrepancy=0, model="multinom",
                 prob=P, logit=Logit)
  if (se) result <- c(result, list(se.prob=SE.P, se.logit=SE.logit,
                                   lower.logit=Lower.logit, upper.logit=Upper.logit, 
                                   lower.prob=Lower.P, upper.prob=Upper.P,
                                   confidence.level=confidence.level))
  # find empty cells, if any, and correct
  whichFact <- unlist(lapply(result$variables, function(x) x$is.factor))
  zeroes <- NULL
  if(sum(whichFact) > 1){
    nameFact <- names(whichFact)[whichFact]
    counts <- xtabs(as.formula( paste("~", paste(nameFact, collapse="+"))), 
                    model.frame(mod))
    zeroes <- which(counts == 0)  
  }
  if(length(zeroes) > 0){
    levs <- expand.grid(lapply(result$variables, function(x) x$levels)) 
    good <- rep(TRUE, dim(levs)[1])
    for(z in zeroes){
      good <- good &  
        apply(levs, 1, function(x) !all(x == levs[z, whichFact]))
    }
    result$prob[!good, ] <- NA
    result$logit[!good, ] <- NA 
    if (se){
      result$se.prob[!good, ] <- NA
      result$se.logit[!good, ] <- NA
      result$lower.prob[!good, ] <- NA
      result$upper.prob[!good, ] <- NA
    }
  } 
  class(result) <-'effpoly'
  result
}

Effect.polr <- function(focal.predictors, mod, 
                        confidence.level=.95, xlevels=list(), default.levels=NULL,
                        given.values, vcov.=vcov, se=TRUE, typical=mean, latent=FALSE, ...){
  if (mod$method != "logistic") stop('method argument to polr must be "logistic"')    
  if (missing(given.values)) given.values <- NULL
  else if (!all(which <- names(given.values) %in% names(coef(mod)))) 
    stop("given.values (", names(given.values[!which]),") not in the model")
  formula.rhs <- formula(mod)[c(1, 3)]
  model.components <- Analyze.model(focal.predictors, mod, xlevels, default.levels, formula.rhs, typical=typical)
  excluded.predictors <- model.components$excluded.predictors
  predict.data <- model.components$predict.data
  factor.levels <- model.components$factor.levels
  factor.cols <- model.components$factor.cols
  #    n.focal <- model.components$n.focal
  x <- model.components$x
  X.mod <- model.components$X.mod
  cnames <- model.components$cnames
  X <- model.components$X
  Terms <- delete.response(terms(mod))
  mf <- model.frame(Terms, predict.data, xlev = factor.levels, na.action=NULL)
  mod.matrix <- model.matrix(formula.rhs, data = mf, contrasts.arg = mod$contrasts)
  X0 <- Fixup.model.matrix(mod, mod.matrix, model.matrix(mod), 
                           X.mod, factor.cols, cnames, focal.predictors, excluded.predictors, typical, given.values)
  resp.names <- make.names(mod$lev, unique=TRUE)
  X0 <- X0[,-1, drop=FALSE]
  b <- coef(mod)
  p <- length(b)  # corresponds to p - 1 in the text
  alpha <- - mod$zeta  # intercepts are negatives of thresholds
  z <- qnorm(1 - (1 - confidence.level)/2)
  result <- list(term=paste(focal.predictors, collapse="*"), formula=formula(mod), response=response.name(mod),
                 y.levels=mod$lev, variables=x, 
                 x=predict.data[, focal.predictors, drop=FALSE],
                 model.matrix=X0, data=X, discrepancy=0, model="polr")
  if (latent){
    res <- eff.latent(X0, b, vcov.(mod)[1:p, 1:p], se)
    result$fit <- res$fit
    if (se){
      result$se <- res$se
      result$lower <- result$fit - z*result$se
      result$upper <- result$fit + z*result$se
      result$confidence.level <- confidence.level
    }
    transformation <- list()
    transformation$link <- I
    transformation$inverse <- I
    result$transformation <- transformation
    result$thresholds <- -alpha
    class(result) <- c("efflatent", "eff")
    return(result)
  }
  m <- length(alpha) + 1
  r <- m + p - 1
  indices <- c((p+1):r, 1:p)
  V <- vcov.(mod)[indices, indices]
  for (j in 1:(m-1)){  # fix up the signs of the covariances
    V[j,] <- -V[j,]  #  for the intercepts
    V[,j] <- -V[,j]}	
  n <- nrow(X0)
  P <- Logit <- matrix(0, n, m)
  colnames(P) <-  paste("prob.", resp.names, sep="")
  colnames(Logit) <-  paste("logit.", resp.names, sep="")
  if (se){
    Lower.logit <- Upper.logit <- Lower.P <- Upper.P <- SE.P <- SE.Logit <- matrix(0, n, m)
    colnames(Lower.logit) <-  paste("L.logit.", resp.names, sep="")
    colnames(Upper.logit) <-  paste("U.logit.", resp.names, sep="")
    colnames(Lower.P) <-  paste("L.prob.", resp.names, sep="")
    colnames(Upper.P) <-  paste("U.prob.", resp.names, sep="")
    colnames(SE.P) <-  paste("se.prob.", resp.names, sep="")
    colnames(SE.Logit) <-  paste("se.logit.", resp.names, sep="")
  }
  for (i in 1:n){
    res <- eff.polr(X0[i,], b, alpha, V, m, r, se) # compute effects
    #        P[i,] <- prob <- res$p # fitted probabilities
    P[i,] <- res$p # fitted probabilities
    Logit[i,] <- logit <- res$logits # fitted logits
    if (se){
      #            SE.P[i,] <- se.p <- res$std.err.p # std. errors of fitted probs
      SE.P[i,] <- res$std.err.p # std. errors of fitted probs
      SE.Logit[i,] <- se.logit <- res$std.error.logits # std. errors of logits
      Lower.P[i,] <- logit2p(logit - z*se.logit)
      Upper.P[i,] <- logit2p(logit + z*se.logit)
      Lower.logit[i,] <- logit - z*se.logit
      Upper.logit[i,] <- logit + z*se.logit
    }
  }
  result$prob <- P
  result$logit <- Logit
  if (se) result <- c(result,
                      list(se.prob=SE.P, se.logit=SE.Logit,
                           lower.logit=Lower.logit, upper.logit=Upper.logit, 
                           lower.prob=Lower.P, upper.prob=Upper.P,
                           confidence.level=confidence.level))
  class(result) <-'effpoly'
  result
}

Effect.default <- function(focal.predictors, mod, xlevels = list(), default.levels = NULL, given.values,
                           vcov. = vcov, se = TRUE, confidence.level = 0.95, 
                           transformation = list(link = I, inverse = I), 
                           typical = mean, offset = mean, ...){
  if (missing(given.values)) 
    given.values <- NULL
  else if (!all(which <- names(given.values) %in% names(coef(mod)))) 
    stop("given.values (", names(given.values[!which]), ") not in the model")
  off <- if (is.numeric(offset) && length(offset) == 1) offset
  else if (is.function(offset)) {
    mod.off <- model.offset(model.frame(mod))
    if (is.null(mod.off)) 0 else offset(mod.off)
  }
  else stop("offset must be a function or a number")
  formula.rhs <- formula(mod)[[3]]
  model.components <- Analyze.model(focal.predictors, mod, xlevels, default.levels, formula.rhs, typical=typical)
  excluded.predictors <- model.components$excluded.predictors
  predict.data <- model.components$predict.data
  factor.levels <- model.components$factor.levels
  factor.cols <- model.components$factor.cols
  n.focal <- model.components$n.focal
  x <- model.components$x
  X.mod <- model.components$X.mod
  cnames <- model.components$cnames
  X <- model.components$X
  formula.rhs <- formula(mod)[c(1, 3)]
  Terms <- delete.response(terms(mod))
  mf <- model.frame(Terms, predict.data, xlev = factor.levels, na.action=NULL)
  mod.matrix <- model.matrix(formula.rhs, data = mf, contrasts.arg = mod$contrasts)
  mod.matrix <- Fixup.model.matrix(mod, mod.matrix, model.matrix(mod), 
                                   X.mod, factor.cols, cnames, focal.predictors, excluded.predictors, typical, given.values)
  mod.matrix <- mod.matrix[, !is.na(coef(mod))]
  effect <- off + mod.matrix %*% mod$coefficients[!is.na(coef(mod))]
  result <- list(term = paste(focal.predictors, collapse="*"), 
                 formula = formula(mod), response = response.name(mod), 
                 variables = x, fit = effect, x = predict.data[, 1:n.focal, drop=FALSE], model.matrix = mod.matrix, data = X, 
                 discrepancy = 0, offset=off)
  whichFact <- unlist(lapply(result$variables, function(x) x$is.factor))
  zeroes <- NULL
  if(sum(whichFact) > 1){
    nameFact <- names(whichFact)[whichFact]
    counts <- xtabs(as.formula( paste("~", paste(nameFact, collapse="+"))), 
                    model.frame(mod))
    zeroes <- which(counts == 0)  
  }
  if(length(zeroes) > 0){
    levs <- expand.grid(lapply(result$variables, function(x) x$levels)) 
    good <- rep(TRUE, dim(levs)[1])
    for(z in zeroes){
      good <- good &  
        apply(levs, 1, function(x) !all(x == levs[z, whichFact]))
    } 
    result$fit[!good] <- NA
  } 
  if (se) {
    z <- qnorm(1 - (1 - confidence.level)/2)
    V <- vcov.(mod)
    eff.vcov <- mod.matrix %*% V %*% t(mod.matrix)
    rownames(eff.vcov) <- colnames(eff.vcov) <- NULL
    var <- diag(eff.vcov)
    result$vcov <- eff.vcov    	
    result$se <- sqrt(var)
    result$lower <- effect - z * result$se
    result$upper <- effect + z * result$se
    result$confidence.level <- confidence.level
    if(length(zeroes) > 0){
      result$se[!good] <- NA
      result$lower[!good] <- NA
      result$upper[!good] <- NA
    }
  }
  result$transformation <- transformation
  class(result) <- "eff"
  result
}