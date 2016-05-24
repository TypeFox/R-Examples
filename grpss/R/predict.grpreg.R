#' Model predictions for a fitted \code{grpreg} object
#' @description Similar to usual predict methods and \code{predict.grprep} in \code{grpreg} package.
#' @param object A fitted "\code{grpreg}" object from \code{\link{grpss}},
#' or \code{\link[grpreg]{grpreg}} function.
#' @param newdata Optionally, a matrix or data frame where to predict. If omits, the fitted
#' predictors are used.
#' @param lambda Value of the regularization parameter \code{lambda} at which predictions are
#' requested. See details for the default.
#' @param type The type of prediction: "\code{response}" gives the fitted values; \code{"class"}
#' returns the predicted class for the binomial outcome; "\code{probability}" returns the
#' predicted probabilities for the logistic regression.
#' @param ... Not used.
#' @details This function gives the predictions at \code{newdata} or all predictors if the
#' argument \code{newdata} is not supplied. The default \code{lambda} for "\code{grpreg}"
#' object is the one at which we obtain the minimum loss value, i.e., negative log-likelihood
#' value. Typically, \code{type = "response"} is
#' used for linear or poisson regression, and \code{type = "class"} or
#' \code{type = "probability"} is used for logistic regression.
#' @return The predicted values depending on the type.
#' @author Debin Qiu, Jeongyoun Ahn
#' @seealso \code{\link{grpss}}
#' @examples
#' library(MASS)
#' set.seed(23)
#' n <- 30 # sample size
#' p <- 3  # number of predictors in each group
#' J <- 50  # group size
#' group <- rep(1:J,each = 3)  # group indices
#' X <- mvrnorm(n,seq(0,5,length.out = p*J),diag(p*J))
#' beta <- runif(12,-2,5)
#' mu <- X%*%matrix(c(beta,rep(0,p*J-12)),ncol = 1)
#'
#' # linear regression with family = "gaussian"
#' y <-  mu + rnorm(n)
#'
#' ## without cross-validation
#' gss12 <- grpss(X,y,ncut = 10,group,select = TRUE)
#' predict(gss12) # fitted values
#' predict(gss12,lambda = 0.2) # fitted values at lambda = 0.2
#'
#' # logistic regression with family = "binomial"
#' set.seed(23)
#' y1 <- rbinom(n,1,1/(1 + exp(-mu)))
#' gss21 <- grpss(X,y1,group, criterion = "gDC",select = TRUE,
#'                family = "binomial")
#' predict(gss21)
#'
#' @export

predict.grpreg <- function(object, newdata, lambda = NULL,
                           type = c("response","class","probability"),...) {
  type <- match.arg(type)
  if (class(object) == "grpss")
    stop("No prediction method available for class 'grpss'
         without doing grouped variable selection")
  callArg <- object$call
  if (is.null(callArg$formula)) {
    y <- eval(callArg$y)
    X <- eval(callArg$X)
  }
  else {
    data <- eval(callArg$data)
    yresp <- as.character(callArg$formula)[2]
    y <- data[,yresp]
    data[,yresp] <- NULL
    X <- data
  }
  group <- eval(object$call$group)
  grp.scr <- object$group.screen
  XX <- as.matrix(cbind(rep(1,length(y)),X[,group %in% grp.scr]))
  if (class(object) == "grpreg") {
    index <- ifelse(is.null(lambda), which.min(object$loss),
                    which.min(abs(object$lambda - lambda)))
    lambda <- ifelse(is.null(lambda),object$lambda[index],lambda)
    if (lambda > max(object$lambda) || lambda < min(object$lambda))
      stop(paste("please specify 'lambda' between",min(object$lambda),
                 "and",max(object$lambda)))
    beta <- object$beta[,index]
  }
  else {
    index <- which.min(abs(object$lambda - object$lambda.min))
    beta <- object$fit$beta[,index]
  }
  if (missing(newdata))
    yhat <- XX%*%matrix(beta,ncol = 1)
  else {
    if (class(newdata) != "matrix") {
      temp <- try(newdata <- as.matrix(newdata), silent = TRUE)
      if (class(temp)[1] == "try-error")
        stop("'newdata' must be a matrix or can be coerced to a matrix")
    }
    if (any(c(NCOL(newdata),NROW(newdata)) == 1))
      newdata <- matrix(c(1,newdata),nrow = 1)
    else
      newdata <- cbind(rep(1,nrow(newdata)),newdata)
    beta.all <- numeric(ncol(X))
    beta.all[group %in% grp.scr] = beta[-1]
    beta.all <- c(beta[1],beta.all)
    yhat <- newdata%*%matrix(beta.all,ncol = 1)
  }
  colnames(yhat) <- "response"
  family <- if(class(object) == "grpreg") object$family else object$fit$family
  if (family == "binomial") {
    prob <- exp(yhat)/(1 + exp(yhat))
    yhat <- matrix(as.numeric(yhat > 0.5), ncol = 1,dimnames = list(NULL,"class"))
    if (type == "probability")
      yhat <- matrix(prob,ncol = 1,dimnames = list(NULL,"probability"))
  }
  if (family == "poisson")
    yhat <- matrix(exp(yhat),ncol = 1,dimnames = list(NULL,"mean"))
  return(yhat)
}
