SL.glm.condExpX2givenW <- function#SL  Wrapper for Estimation of Cond. Expect. of X^2 Given W
### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional expectation of \eqn{X^2} given \eqn{W}.
(Y, X, newX, family, obsWeights, ...) {
  ##seealso<< learnCondExpX2givenW, predict.SL.glm.condExpX2givenW
  varNames <- names(X)
  if (length(varNames)>20) {
    varNames <- varNames[1:20]
    warning(paste("Using only", paste(varNames, collapse=", "), "in 'SL.glm.condExpX2givenW'"))
  }
  theFormula <- paste(varNames, collapse=" + ")
  if (length(varNames)<=10) {
    theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
    theFormula <- paste("Y ~", theFormula, "+", theFormula2, sep="")
  } else {
    theFormula <- paste("Y ~", theFormula, sep="")
  } 
  ## formula.glm.condExpX2givenW <- as.formula(Y~W+I(W^2));
  formula.glm.condExpX2givenW <- as.formula(theFormula);
  
  fit.glm <- glm(formula.glm.condExpX2givenW, data = X, family = family,
                 weights = obsWeights)
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- c("SL.glm.condExpX2givenW")
  out <- list(pred = pred, fit = fit)
  return(out)
### Returns a fitted object.
}


predict.SL.glm.condExpX2givenW <- function#SL  Wrapper for Estimation of Cond. Expect. of X^2 Given W
### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional expectation of \eqn{X^2} given \eqn{W}.
(object, newdata, ...) {
  ##seealso<< SL.glm.condExpX2givenW
  out <- predict(object = object$object, newdata = newdata, 
                 type = "response")
  return(out)
### Returns a prediction.
}


SL.glm.condExpXYgivenW <- function#SL  Wrapper for Estimation of Cond. Expect. of XY Given W
### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional expectation of \eqn{XY} given \eqn{W}.
(Y, X, newX, family, obsWeights, ...) {
  ##seealso<< learnCondExpXYgivenW, predict.SL.glm.condExpXYgivenW
  varNames <- names(X)
  if (length(varNames)>20) {
    varNames <- varNames[1:20]
    warning(paste("Using only", paste(varNames, collapse=", "), "in 'SL.glm.condExpXYgivenW'", collapse=""))
  }
  theFormula <- paste(varNames, collapse=" + ")
  if (length(varNames)<=10) {
    theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
    theFormula <- paste("Y ~", theFormula, "+", theFormula2, sep="")
  } else {
    theFormula <- paste("Y ~", theFormula, sep="")
  }
  ## formula.glm.condExpXYgivenW <- as.formula(Y~W+I(W^2));
  formula.glm.condExpXYgivenW <- as.formula(theFormula);
  
  fit.glm <- glm(formula.glm.condExpXYgivenW, data = X, family = family,
                 weights = obsWeights)
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- c("SL.glm.condExpXYgivenW")
  out <- list(pred = pred, fit = fit)
  return(out)
### Returns a fitted object.
}


predict.SL.glm.condExpXYgivenW <- function#SL  Wrapper for Estimation of Cond. Expect. of XY Given W
### Prediction algorithm wrapper for SuperLearner.
(object, newdata, ...) {
  ##seealso<< SL.glm.condExpXYgivenW
  out <- predict(object = object$object, newdata = newdata, 
                 type = "response")
  return(out)
### Returns a prediction.
}
## environment(predict.SL.glm.condExpXYgivenW) <- asNamespace("SuperLearner")


SL.glm.g <- function#SL  Wrapper for Estimation of Cond. Prob. of X=0 Given W
### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional probability of \eqn{X=0} given \eqn{W}.
(Y, X, newX, family, obsWeights, ...) {
  ##seealso<< learnG, predict.SL.glm.g
  varNames <- names(X)
  if (length(varNames)>20) {
    varNames <- varNames[1:20]
    warning(paste("Using only", paste(varNames, collapse=", "), "in 'SL.glm.g'"))
  }
  theFormula <- paste(varNames, collapse=" + ")
  if (length(varNames)<=10) {
    theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
    theFormula <- paste("Y ~", theFormula, "+", theFormula2, sep=" ")
  } else {
    theFormula <- paste("Y ~", theFormula, sep="")
  }
  ## formula.glm.g <- as.formula(Y~W+I(W^2));
  formula.glm.g <- as.formula(theFormula);



  fit.glm <- glm(formula.glm.g, data = X, family = family,
                 weights = obsWeights)
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- c("SL.glm.g")
  out <- list(pred = pred, fit = fit)
  return(out)
### Returns a fitted object.
}
## environment(SL.glm.g) <- asNamespace("SuperLearner")


predict.SL.glm.g <- function#SL  Wrapper for Estimation of Cond. Prob. of X=0 Given W
### Prediction algorithm wrapper for SuperLearner.
(object, newdata, ...) {
  ##seealso<< SL.glm.g
    out <- predict(object = object$object, newdata = newdata, 
        type = "response")
    return(out)
### Returns a prediction.
}
## environment(predict.SL.glm.g) <- asNamespace("SuperLearner")


SL.glm.theta <- function#SL  Wrapper for Estimation of Cond. Prob. of X=0 Given W
### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional expectation of \eqn{Y} given \eqn{(X,W)}.
(Y, X, newX, family, obsWeights, ...) {
  ##seealso<< learnTheta, predict.SL.glm.theta
  varNames <- names(X)
  if (length(varNames)>20) {
    varNames <- varNames[1:20]
    warning(paste("Using only 'X' and", paste(varNames, collapse=", "), "in 'SL.glm.theta'"))
  }
  theFormula <- paste(varNames, collapse=" + ")
  if (length(varNames)<=10) {
    theFormula <- paste("Y ~ X+", theFormula, "+ X*(",
                        theFormula, ")", sep="")
  } else {
    theFormula <- paste("Y ~ X+", theFormula, sep="")
  }
  ## formula.glm.theta <- as.formula(Y~X*W);
  formula.glm.theta <- as.formula(theFormula);

  fit.glm <- glm(formula.glm.theta, data = X, family = family, 
                 weights = obsWeights)
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- c("SL.glm.theta")
  out <- list(pred = pred, fit = fit)
  return(out)
### Returns a fitted object.
}
## environment(SL.glm.theta) <- asNamespace("SuperLearner")



predict.SL.glm.theta <- function#SL  Wrapper for Estimation of Cond. Expect. of Y Given (X,W)
### Prediction algorithm wrapper for SuperLearner.
(object, newdata, ...) {
  ##seealso<< SL.glm.theta
    out <- predict(object = object$object, newdata = newdata, 
        type = "response")
    return(out)
### Returns a prediction.
}

