##----------------------------
## methods for mlogit objects |
##----------------------------
##    * fitted                |
##    * residuals             |
##    * df.residual           |
##    * terms                 |
##    * model.matrix          |
##    * model.response        |
##    * update                |
##    * print                 |
##    * vcov                  |
##    * logLik                |
##    * summary               |
##    * print.summary         |
##    * index                 |
##    * predict               |
##    * coef                  |
##----------------------------

residuals.mlogit <- function(object, outcome = TRUE, ...){
  if (!outcome){
    result <- object$residuals
  }
  else{
    J <- ncol(object$residuals)
    y <- matrix(model.response(object$model), ncol = J, byrow = T)
    result <- apply(y * object$residuals, 1, sum)
  }
  result
}

df.residual.mlogit <- function(object, ...){
  n <- length(residuals(object))
  K <- length(coef(object))
  n-K
}

terms.mlogit <- function(x, ...){
  terms(x$formula)
}

model.matrix.mlogit <- function(object, ...){
  model.matrix(object$formula, object$model)
}

model.response.mlogit <- function(object, ...){
  y.name <- paste(deparse(object$formula[[2]]))
  object$model[[y.name]]
}

update.mlogit <- function (object, new, ...){
  call <- object$call
  if (is.null(call))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(new))
    call$formula <- update(formula(object), new)
  if(length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    ## do these individually to allow NULL to remove entries.
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  eval(call, parent.frame())
}

print.mlogit <- function (x, digits = max(3, getOption("digits") - 2),
                          width = getOption("width"), ...){
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

vcov.mlogit <- function(object, what = c('coefficient', 'errors', 'rpar'),
                        type = c('cov', 'cor', 'sd'), reflevel = NULL, ...){
  what <- match.arg(what)
  type <- match.arg(type)
  if (what == 'coefficient'){
    fixed <- attr(object$coefficients, "fixed")
    result <- solve(-object$hessian[!fixed, !fixed])
  }
  if (what == 'errors'){
    if (!is.null(object$omega)){
      if (is.null(reflevel)){
        if (is.list(object$omega)) result <- object$omega[[1]]
        else result <- object$omega
      }
      else result <- object$omega[[reflevel]]
    }
    result <- switch(type,
                     cov = result,
                     cor = result / tcrossprod(sqrt(diag(result))),
                     sd = sqrt(diag(result))
                     )
  }
  if (what == 'rpar'){
    if (is.null(object$rpar)) stop('no random parameters')
    if (is.null(attr(object$rpar, "covariance"))){
      result <- stdev(object)
      if (type != 'sd'){
        V  <- matrix(0, length(result), length(result), dimnames = list(names(result), names(result)))
        if (type == 'cor') diag(V) <- 1
        if (type == 'cov') diag(V) <- result^2
        result <- V
      }
    }
    else{
      thecov <- attr(object$rpar, "covariance")
      result <- switch(type,
                       cov = thecov,
                       cor = thecov / tcrossprod(sqrt(diag(thecov))),
                       sd = sqrt(diag(thecov))
                       )
    }
  }
  result
}

logLik.mlogit <- function(object,...){
  object$logLik
}

summary.mlogit <- function (object,...){
  fixed <- attr(object$coefficients, "fixed")
  b <- coef(object)[!fixed]
  std.err <- sqrt(diag(vcov(object)))
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "t-value", "Pr(>|t|)")
  object$CoefTable <- CoefTable
  if (has.intercept(object$formula)){
    object$lratio <- lratio(object)
    object$mfR2 <- mfR2(object)
  }
  if (!is.null(object$rpar)){
    rpar <- object$rpar
    object$summary.rpar <- t(sapply(rpar, summary))
  }

  class(object) <- c("summary.mlogit", "mlogit")
  return(object)
}

print.summary.mlogit <- function(x, digits = max(3, getOption("digits") - 2),
                                 width = getOption("width"), ...){
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  cat("Frequencies of alternatives:")
  print(prop.table(x$freq), digits = digits)
  cat("\n")
  print(x$est.stat)
  cat("\nCoefficients :\n")
  printCoefmat(x$CoefTable, digits = digits)
  cat("\n")
  cat(paste("Log-Likelihood: ", signif(x$logLik, digits), "\n", sep = ""))
  if (has.intercept(x$formula)){
    cat("McFadden R^2: ", signif(x$mfR2, digits), "\n")
    cat("Likelihood ratio test : ", names(x$lratio$statistic),
        " = ", signif(x$lratio$statistic, digits),
        " (p.value = ", format.pval(x$lratio$p.value, digits = digits), ")\n", sep = "")
  }
  if (!is.null(x$summary.rpar)){
    cat("\nrandom coefficients\n")
    print(x$summary.rpar)
  }
  invisible(x)
}

index.mlogit <- function(x, ...){
  index(model.frame(x))
}

predict.mlogit <- function(object, newdata = NULL, returnData = FALSE, ...){
  # if no newdata is provided, use the mean of the model.frame
  if (is.null(newdata)) newdata <- mean(model.frame(object))
  # if newdata is not a mlogit.data, it is coerced below
  if (!inherits(newdata, "mlogit.data")){
    rownames(newdata) <- NULL
    lev <- colnames(object$probabilities)
    J <- length(lev)
    choice.name <- attr(model.frame(object), "choice")
    if (nrow(newdata) %% J)
      stop("the number of rows of the data.frame should be a multiple of the number of alternatives")
    attr(newdata, "index") <- data.frame(chid = rep(1:(nrow(newdata)%/%J), each = J), alt = rep(lev, J))
    if (is.null(newdata[['choice.name']])){
      newdata[[choice.name]] <- FALSE
      newdata[[choice.name]][1] <- TRUE # probit and hev requires that one (arbitrary) choice is TRUE
    }
  }
  # if the updated model requires the use of mlogit.data, suppress all
  # the relevant arguments
  m <- match(c("choice", "shape", "varying", "sep",
               "alt.var", "chid.var", "alt.levels",
               "opposite", "drop.index", "id", "ranked"),
             names(object$call), 0L)
  if (sum(m) > 0) object$call <- object$call[ - m]
  # update the model and get the probabilities
  newobject <- update(object, start = coef(object), data = newdata, iterlim = 0, print.level = 0)
  result <- newobject$probabilities
  if (nrow(result) == 1){
    result <- as.numeric(result)
    names(result) <- colnames(object$probabilities)
  }
  if (returnData) attr(result, "data") <- newdata
  result
}

fitted.mlogit <- function(object, outcome = TRUE, ...){
  if (outcome) result <- object$fitted
  else result <- object$probabilities
  result
}

coef.mlogit <- function(object, ...){
  result <- object$coefficients
  attr(result, "fixed") <- NULL
  result
}

mean.mlogit.data <- function(x, ...){
  alt <- index(x)$alt
  J <- length(levels(alt))
  result <- data.frame(lapply(x,
                              function(x){
                                if (is.numeric(x)) result <- as.numeric(tapply(x, alt, mean))
                                else{
                                  if (is.logical(x)){
                                    z <- tapply(x, alt, sum)
                                    result <- z == max(z)
                                  }
                                  if(is.character(x)){
                                    x <- factor(x, levels = unique(x))
                                  }
                                  if (is.factor(x)){
                                    result <- factor(names(which.max(table(x))), levels = levels(x))
                                  }
                                }
                                result
                              }
                              )
                       )
  attr(result, "index") <- data.frame(alt = factor(levels(alt), levels =  levels(alt)), chid = rep(1, J))
  rownames(result) <- rownames(attr(result, "index")) <- paste(rep(1, J), levels(alt), sep = ".")
  class(result) <- c("mlogit.data", "data.frame")
  result
}


effects.mlogit <- function(object, covariate = NULL,
                           type = c("aa", "ar", "rr", "ra"),
                           data = NULL, ...){
  type <- match.arg(type)
  if (is.null(data)){
    P <- predict(object, returnData = TRUE)
    data <- attr(P, "data")
    attr(P, "data") <- NULL
  }
  else P <- predict(object, data)
  newdata <- data
  J <- length(P)
  alt.levels <- names(P)
  pVar <- substr(type, 1, 1)
  xVar <- substr(type, 2, 2)
  cov.list <- lapply(attr(formula(object), "rhs"), as.character)
  rhs <- sapply(cov.list, function(x) length(na.omit(match(x, covariate))) > 0)
  rhs <- (1:length(cov.list))[rhs]
  eps <- 1E-5

  if (rhs %in% c(1, 3)){
    if (rhs == 3){
      theCoef <- paste(alt.levels, covariate, sep = ":")
      theCoef <- coef(object)[theCoef]
    }
    else theCoef <- coef(object)[covariate]
    me <- c()
    for (l in 1:J){
      newdata[l, covariate] <- data[l, covariate] + eps
      newP <- predict(object, newdata)
      me <- rbind(me, (newP - P) / eps)
      newdata <- data
    }
    if (pVar == "r") me <- t(t(me) / P)
    if (xVar == "r") me <- me * matrix(rep(data[[covariate]], J), J)
    dimnames(me) <- list(alt.levels, alt.levels)
  }
  if (rhs == 2){
    newdata[, covariate] <- data[, covariate] + eps
    newP <- predict(object, newdata)
    me <- (newP - P) / eps
    if (pVar == "r") me <- me / P
    if (xVar == "r") me <- me * data[[covariate]]
    names(me) <- alt.levels
  }
  me
}
