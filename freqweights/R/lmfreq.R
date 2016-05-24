## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2015-06-21 09:22 emilio on emilio-despacho>
## ============================================================

##' \code{lmfreq} is used to fit linear models with frequency tables
##' 
##' To fit linear models with data grouped in frequency tables. 
##' 
##' It computes the linear model of a frequency table. See
##' \code{\link[stats]{lm}} for further details.
##'
##' Any variables in the formula are removed from the data set.
##'
##' The dot function are for programming purpose. It does not check the data.
##' 
##' @aliases lmfreq AIC.lmfreq extractAIC.lmfreq logLik.lmfreq print.lmfreq
##' summary.lmfreq print.summary.lmfreq predict.lmfreq
##' @param formula an object of class \code{formula} 
##' @param data a data frame that must contain all variables present in
##' \code{formula} and \code{freq}
##' @param freq a character string specifying the variable of  frequency weights
##' @param tfq a \code{tablefreq} object
##' @return It returns an object of class \code{lmfreq}, very similar to \code{\link[stats]{lm}}
##' @export
##' @seealso \code{\link{tablefreq}}
##' @examples
##'
##' ## Benchmark
##' if(require(hflights)){
##'   formula <-  ArrDelay ~ DepDelay   
##'   print(system.time(a <- lm(formula, data=hflights)))  ## ~0.4 seconds 
##'   print(system.time(b <- lmfreq(formula, data=hflights))) ## ~0.12 seconds. 4x faster
##' }
##' 
##' l0 <- lm(Sepal.Length ~ Sepal.Width,iris)
##' summary(l0)
##' 
##' tfq <- tablefreq(iris[,1:2])
##' lf <- lmfreq(Sepal.Length ~ Sepal.Width,tfq, freq="freq")
##' summary(lf)
##' 
##' all.equal(coef(lf),coef(l0))
##' all.equal(AIC(lf),AIC(l0))
##' 
##' newdata <- data.frame(Sepal.Width=c(1,NA,7))
##' predict(lf, newdata)
##' 
##' if(require(MASS)){
##'    stepAIC(lf)
##' }
##' 
##' system.time(lmfreq(Sepal.Length ~ Sepal.Width,tfq, freq="freq"))
##' system.time(.lmfreq(Sepal.Length ~ Sepal.Width,tfq)) # Fast
##'
##' library(dplyr)
##' igrouped <- iris %>% group_by(Species)
##' models <- igrouped %>% do(model=lmfreq(Sepal.Length ~ Sepal.Width, .))
##' coefs <- models %>%
##'   do(cbind(as.data.frame(rbind(coef(.$model))),
##'            Species=.$Species))
##' coefs
##' 
##' \dontrun{
##' ## If data is too granular, benchmark is worst
##' n <- 10^6
##' data <- data.frame(y=rnorm(n),x=rnorm(n))
##' system.time(lm(y~x,data)) ## ~5 seconds
##' system.time(lmfreq(y~x,data)) ## ~ 15 seconds
##' system.time(tfq <- tablefreq(data)) ## ~ 5 seconds
##' nrow(tfq) # same number of rows than original data
##' system.time(.lmfreq(y~x,tfq)) ## ~ 10 seconds
##' }
lmfreq<- function(formula, data, freq=NULL) {
  cl <- match.call()
  m <- .lmfreq(formula,
                     tablefreq(data, vars=c(all.vars(formula)), freq=freq))
  m$freq <- freq
  m$call <- cl
  m
}

##' @rdname lmfreq
##' @export
.lmfreq<- function(formula, tfq) {
  cl <- match.call()
   ## x <- as.data.frame(tablefreq(data, vars=c(all.vars(formula)), freq=freq))
##  x <- as.data.frame(tablefreq(data, freq=freq))
  m <- do.call("lm", list(formula=formula,
                          data=tfq,
                          weights=unlist(tfq[,ncol(tfq)])))
  m <- m[c("coefficients", "rank", "assign",
           "xlevels", "terms", "qr",
           "residuals","weights",
           "model")]
  m$freq <- formula(paste("~", colnames(tfq)[ncol(tfq)]))
  m$call <- cl
  class(m) <- c("lmfreq")
  m
}


## ============================================================
##
## ============================================================

##' @param object a  \code{lmfreq} object
##' @param ... See Details 
##' @rdname lmfreq
##' @method logLik lmfreq
##' @export
logLik.lmfreq <- function (object, ...)
{
  if (!inherits(object, "lmfreq"))
    stop("logLik.lmfreq: error")
  REML <- FALSE
  res <- object$residuals
  p <- object$rank
  ## N <- length(res)
  N <- if (!is.null(w <- object$weights)) sum(w != 0) else nrow(object$residuals)
  if (is.null(w <- object$weights)) {
    w <- rep.int(1, N)
  }        else {
    excl <- w == 0
    if (any(excl)) {
      res <- res[!excl]
      ##N <- length(res)
      w <- w[!excl]
    }
  }
  ## Now sum the w
  N <- sum(w)
  N0 <- N
  if (REML)
    N <- N - p
  ## val <- 0.5 * (sum(log(w)) - N * (log(2 * pi) + 1 - log(N) +
  ##     log(sum(w * res^2))))
  val <- 0.5 * ( 0 - N * (log(2 * pi) + 1 - log(N) +
                          log(sum(w * res^2))))
  attr(val, "nall") <- N0
  attr(val, "nobs") <- N
  attr(val, "df") <- p + 1
  class(val) <- "logLik"
  val
}


## ============================================================
##
## ============================================================

## ============================================================
##
## ============================================================


##' @param fit a \code{lmfreq} object 
##' @param scale not used
##' @param k penalty parameter
##' @method extractAIC lmfreq
##' @rdname lmfreq
##' @export
extractAIC.lmfreq <- function (fit, scale = 0, k = 2,...) {
  if (!inherits(fit, "lmfreq"))
    stop("logLik.lmfreq: error")
  res <- fit$residuals
  edf <- fit$rank
  ## N <- length(res)
  N <- if (!is.null(w <- fit$weights)) sum(w != 0) else nrow(fit$residuals)
  if (is.null(w <- fit$weights)) {
    w <- rep.int(1, N)
  }        else {
    excl <- w == 0
    if (any(excl)) {
      res <- res[!excl]
      ##N <- length(res)
      w <- w[!excl]
    }
  }
  ## Now sum the w
  N <- sum(w)
  dev <-  N * (log( sum(w * res^2)) - log( N))
  c(edf, dev + 2 * edf)
}





##' @method AIC lmfreq
##' @rdname lmfreq 
##' @export
AIC.lmfreq <- function(object, ..., k=2){
  ## See http://stat.ethz.ch/R-manual/R-patched/library/stats/html/extractAIC.html
  m <- object
  if(!inherits(m, "lmfreq")){
    stop("AIC.lmfreq: Only use lmfreq objects")
  }
  ## print(paste("AIC.freq ",usefreq," ", is.logical(usefreq)))
  N <- if (!is.null(w <- m$weights)) sum(w != 0) else nrow(m$residuals)
  if(is.null(w <- m$weights)) {
    w <- rep.int(1, N)
  }
  p <- m$rank
  res <- m$residuals
  if( FALSE) {
    aic <- -2 * ( 0.5 * (sum(log(w)) -
                         N * (log(2 * pi) +
                              1 - log(N) +
                              log(sum(w * res^2))))) + k *( p+1)
  } else {
    N <- sum(w)
    aic <- -2 * ( 0.5 * (0  -
                         N * (log(2 * pi) +
                              1 - log(N) +
                              log(sum(w * res^2))))) + k *( p+1)

  }
  aic
}


##' @rdname lmfreq
##' @method nobs lmfreq
##' @importFrom stats nobs
##' @export
nobs.lmfreq <- function(object, ...) {
##  print(paste("nobs.lmfreq ",sum(object$weights)))
  sum(object$weights)
}



##' @rdname lmfreq
##' @method summary lmfreq
##' @export
summary.lmfreq <- function (object, ...){
  beta <- coef(object)
  N <- sum(object$weights)
  p <-  object$rank
  sigma2 <- sum(object$weights * object$residuals^2)/(N-p)
  if(sigma2>0) {
    vcov <- sigma2 * chol2inv(object$qr$qr)
    se <- sqrt(diag(vcov))
  } else {
    se <-  0
  }
  ## mat <- cbind(Coef = beta, `(95%` = beta - 2 * se, `CI)` = beta +
  ##     2 * se, SE = se, p = 2 * pnorm(abs(beta/se), lower.tail = FALSE))
  mat <- cbind(Coef = beta, `(95%` = beta - 2 * se, `CI)` = beta +
               2 * se, SE = se, p = 2 * pt(abs(beta/se), N-p, lower.tail = FALSE))
  rownames(mat) <- names(coef(object))
  rval <- list(obj = object, call=object$call, mat = mat,
               AIC=AIC(object), N=N, sigma=sqrt(sigma2), rdf= N-p)
  class(rval) <- "summary.lmfreq"
  rval
}


##' @param x a \code{lmfreq} object
##' @rdname lmfreq
##' @method print lmfreq
##' @export
print.lmfreq <- function (x, ...)
{
  ##print(x$obj)
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
}


##' @param digits digits
##' @rdname lmfreq
##' @export
print.summary.lmfreq <- function (x, digits = getOption("digits") - 3, ...)
{
  ##print(x$obj)
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  cat("\nCoefficients:\n")
  print(round(x$mat, digits))
  cat("\nN: ",x$N)
  cat("\nResidual standard error:", format(signif(x$sigma,
                                                  digits)), "on", x$rdf, "degrees of freedom")
  cat("\nAIC: ",  round(x$AIC, digits),"\n")
  invisible(x)
}


## ============================================================
##
## ============================================================


##' @rdname lmfreq
##' @method predict lmfreq
##' @export
predict.lmfreq <- function(object, ...){
  dots <- list(...)
  if(length(dots)) {
    
    ## pos <- which("newdata" == names(dots))
    ## if(length(pos) == 0) {
    ##   pos <- 1
    ## } else if(length(pos)>1) {
    ##   pos <- pos[1]
    ## }
    class(object) <- "lm"
    return(predict(object, dots[[1]]))
  } else {
    freq <- object$weights
    class(object) <- "lm"
    x <- cbind(predict(object), freq)
    colnames(x) <- c("yhat","freq")
    return(x)
  }
  }

