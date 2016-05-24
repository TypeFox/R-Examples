Qa <- function(m1, m2) {
  k1 <- attr(logLik(m1),"df")
  k2 <- attr(logLik(m2),"df")
  exp(k1-k2)
}

Qc <- function(m1, m2) {
  n <- length(residuals(m1))
  k1 <- attr(logLik(m1),"df")
  k2 <- attr(logLik(m2),"df")
  exp(k1*(n/(n-k1-1))-k2*(n/(n-k2-1)))
}

Qb <- function(m1, m2) {
  n <- length(residuals(m1))
  Qa(m1, m2)^(log(n)/2)
}

L <- function(m1, m2) ((sd(residuals(m1))^2)/(sd(residuals(m2))^2))^(length(residuals(m1))/2)

#' Likelihood Ratio Test
#' 
#' Computes a likelihood ratio statistic which reflects the relative likelihood of the
#' data given two competing models.
#' 
#' @param object an object. See below for details.
#' @param ... further object specifications passed to methods. See below for details.
#' @param name a function for extracting a suitable name/description from a fitted model object. By default the name is queried by calling formula.
#' 
#' @section Details:
#' lr.glover performs comparisons of models via likelihood ratio tests. The default method 
#' consecutively compares the fitted model object object with the models passed in .... 
#' Subsequently, a likelihood ratio test for each two consecutive models is carried out.
#' 
#' @return An object of class "anova" which contains the log-likelihood, degrees of 
#' freedom, the difference in degrees of freedom, likelihood ratio, and AIC/BIC 
#' corrected likelihood ratios.
#' 
#' @export
#' 
#' @references Glover, S. & Dixon, P. (2004). Likelihood ratios: A simple and flexible statistic for empirical psychologists. Psychonomic Bulletin & Review, 11(5), 791-806.
#' 
#' @examples
#' m1 <- lm(mpg~.,mtcars)
#' m2 <- step(m1,~.,trace=0)
#' m3 <- step(m1,~.+.^2,trace=0)
#' lr.glover(m1,m2,m3)
#' 
lr.glover <- function(object, ..., name = NULL) {
  logLik0 <- if ("stats4" %in% loadedNamespaces()) 
    stats4::logLik
  else logLik
  update0 <- if ("stats4" %in% loadedNamespaces()) 
    stats4::update
  else update
  nobs0 <- function(x, ...) {
    nobs1 <- if ("stats4" %in% loadedNamespaces()) 
      stats4::nobs
    else nobs
    nobs2 <- function(x, ...) NROW(residuals(x, ...))
    rval <- try(nobs1(x, ...), silent = TRUE)
    if (inherits(rval, "try-error") | is.null(rval)) 
      rval <- nobs2(x, ...)
    return(rval)
  }
  cls <- class(object)[1]
  tlab <- function(x) attr(terms(x), "term.labels")
  if (is.null(name)) 
    name <- function(x) {
      rval <- try(formula(x), silent = TRUE)
      if (inherits(rval, "try-error") | is.null(rval)) 
        rval <- try(x$call, silent = TRUE)
      if (inherits(rval, "try-error") | is.null(rval)) 
        return(NULL)
      else return(paste(deparse(rval), collapse = "\n"))
    }
  modelUpdate <- function(fm, update) {
    if (is.numeric(update)) {
      if (any(update < 1)) {
        warning("for numeric model specifications all values have to be >=1")
        update <- abs(update)[abs(update) > 0]
      }
      if (any(update > length(tlab(fm)))) {
        warning(paste("more terms specified than existent in the model:", 
                      paste(as.character(update[update > length(tlab(fm))]), 
                            collapse = ", ")))
        update <- update[update <= length(tlab(fm))]
      }
      update <- tlab(fm)[update]
    }
    if (is.character(update)) {
      if (!all(update %in% tlab(fm))) {
        warning(paste("terms specified that are not in the model:", 
                      paste(dQuote(update[!(update %in% tlab(fm))]), 
                            collapse = ", ")))
        update <- update[update %in% tlab(fm)]
      }
      if (length(update) < 1) 
        stop("empty model specification")
      update <- as.formula(paste(". ~ . -", paste(update, 
                                                  collapse = " - ")))
    }
    if (inherits(update, "formula")) 
      update <- update0(fm, update)
    if (!inherits(update, cls)) 
      warning(paste("original model was of class \"", cls, 
                    "\", updated model is of class \"", class(update)[1], 
                    "\"", sep = ""))
    return(update)
  }
  objects <- list(object, ...)
  nmodels <- length(objects)
  if (nmodels < 2) {
    objects <- c(objects, . ~ 1)
    nmodels <- 2
  }
  no.update <- sapply(objects, function(obj) inherits(obj, 
                                                      cls))
  for (i in 2:nmodels) objects[[i]] <- modelUpdate(objects[[i - 
                                                              1]], objects[[i]])
  ns <- sapply(objects, nobs0)
  if (any(ns != ns[1])) {
    for (i in 2:nmodels) {
      if (ns[1] != ns[i]) {
        if (no.update[i]) 
          stop("models were not all fitted to the same size of dataset")
        else {
          commonobs <- row.names(model.frame(objects[[i]])) %in% 
            row.names(model.frame(objects[[i - 1]]))
          objects[[i]] <- eval(substitute(update(objects[[i]], 
                                                 subset = commonobs), list(commonobs = commonobs)))
          if (nobs0(objects[[i]]) != ns[1]) 
            stop("models could not be fitted to the same size of dataset")
        }
      }
    }
  }
  rval <- matrix(rep(NA, 9 * nmodels), ncol = 9)
  colnames(rval) <- c("#Df","LogLik","AIC","BIC","Df","L","L*Qa","L*Qc","L*Qb")
  rownames(rval) <- 1:nmodels
  logL <- lapply(objects, logLik0)
  rval[, 1] <- as.numeric(sapply(logL, function(x) attr(x, "df")))
  rval[, 2] <- sapply(logL, as.numeric)
  rval[, 3] <- sapply(objects, AIC)
  rval[, 4] <- sapply(objects, BIC)
  rval[2:nmodels, 5] <- rval[2:nmodels, 1] - rval[1:(nmodels - 1), 1]
  rval[2:nmodels, 6] <- sapply(2:nmodels, function(i) L(objects[[i-1]],objects[[i]]))
  rval[2:nmodels, 7] <- sapply(2:nmodels, function(i) Qa(objects[[i-1]],objects[[i]])*L(objects[[i-1]],objects[[i]]))
  rval[2:nmodels, 8] <- sapply(2:nmodels, function(i) Qc(objects[[i-1]],objects[[i]])*L(objects[[i-1]],objects[[i]]))
  rval[2:nmodels, 9] <- sapply(2:nmodels, function(i) Qb(objects[[i-1]],objects[[i]])*L(objects[[i-1]],objects[[i]]))
  variables <- lapply(objects, name)
  if (any(sapply(variables, is.null))) 
    variables <- lapply(match.call()[-1L], deparse)[1L:nmodels]
  title <- "Likelihood ratio test\n"
  topnote <- paste("Model ", format(1:nmodels), ": ", variables, 
                   sep = "", collapse = "\n")
  structure(as.data.frame(rval), heading = c(title, topnote), class = c("anova","data.frame"))
}