## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2014-04-28 Mon 13:48 emilio on emilio-despacho>
## ============================================================

##' Estimates the coefficients of a linear model
##' 
##' Estimates the coefficients of a linear model following the guidelines of
##' \code{\link[biglm]{biglm}}
##'
##' Any variables in the formula are removed from the data set.
##' 
##' It only computes the coefficients of the linear model.
##' 
##' @aliases biglmfreq  coef.biglmfreq predict.biglmfreq update.biglmfreq
##' @param formula a model formula
##' @param data data frame that must contain all variables in \code{formula}
##' and \code{freq} 
##' @param freq a string of the variable specifying frequency weights
##' @return A \code{biglmfreq} object. 
##' @seealso \code{\link[biglm]{biglm}},  \code{\link{make.readchunk}}
##' @import biglm
##' @export
##' @examples
##' mt <- biglmfreq(Sepal.Length ~ Sepal.Width, iris)
##' coef(mt)
##'
##' chunk1 <- iris[1:30,]
##' chunk2 <- iris[-c(1:30),]
##' mf1 <- biglmfreq(Sepal.Length ~ Sepal.Width, chunk1)
##' mf2 <- update(mf1, chunk2)
##' 
##' predict(mf2, iris)
biglmfreq <- function(formula, data, freq = NULL) {
  cl <- match.call()
  tfq <- tablefreq(data,vars=all.vars(formula), freq=freq)
  weights <- formula(paste("~", attr(tfq,"colweights")))
  m <- biglm(formula=formula, data=tfq, weights=weights)
##  m <- biglm(formula=as.formula(formula), data=as.data.frame(tfq), weights=weights)
  ## yhat <- predict(m, data)
  ## residuals <- data[,all.vars(formula)[1] ] - yhat
  m$N <- evaldp(tfq, summarise, paste("N = sum(", attr(tfq,"colweights"),")"))[1,1]
  m$freq <- freq
  m$call <- cl
  class(m) <- "biglmfreq"
  m
}




##' @param object a \code{biglmfreq} object
##' @param ... See Details
##' @export 
##' @import biglm
##' @rdname biglmfreq
##' @export
coef.biglmfreq <- function(object, ...) {
  class(object) <- "biglm"
   coef(object)
   }


##' @details \code{\dots} should be a data frame when \code{predict}. See Examples
##' @method predict biglmfreq
##' @import biglm
##' @rdname biglmfreq
##' @export
predict.biglmfreq <- function(object, ...) {
  newdata <- list(...)[[1]]
    class(object) <- "biglm"
  predict(object,newdata)
}


##' @param x a \code{biglmfreq} object
##' @method print biglmfreq
##' @import biglm
##' @rdname biglmfreq
##' @export
print.biglmfreq <- function (x, ...)
{
  ##print(x$obj)
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
}


##' @details \code{\dots} should be a data frame when \code{update}. See Examples
##' @method update biglmfreq
##' @import biglm
##' @rdname biglmfreq
##' @export
update.biglmfreq <- function (object, ...) {
## Fromo update.biglm
  freq <- object$freq
 N <-  object$N
  class(object) <- "biglm"
  moredata <- list(...)[[1]]
  if(is.null(moredata)) return(object)
  ##print(class(object))
  tfq <- tablefreq(moredata, vars=all.vars(object$terms), freq=object$freq)
  m <- update(object, tfq)
  ## yhat <- predict(m, moredata)
  ##residuals <- moredata[,all.vars(m$terms)[1] ] - yhat
  ##w  <- model.frame(m$weights, moredata)[[1]]
  ##m$ss <-  m$ss + sum(residuals^2 * w)
  m$N <- N + sum(tfq[,ncol(tfq)])
  m$freq <- freq
  class(m) <- c("biglmfreq")
  m
}


