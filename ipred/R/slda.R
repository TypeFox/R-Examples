# $Id: slda.R,v 1.9 2005/06/29 08:50:28 hothorn Exp $

# stabilized linear discriminant analysis according to Laeuter & Kropf

slda <- function(y, ...) UseMethod("slda")

slda.default <- function(y, ...) 
  stop(paste("Do not know how to handle objects of class", class(data)))

slda.formula <- function(formula, data, subset, na.action=na.rpart, ...) {
   cl <- match.call()
   if(missing(formula)
       || (length(formula) != 3)
       || (length(attr(terms(formula[-2], data = data), "term.labels")) < 1)
       || (length(attr(terms(formula[-3], data = data), "term.labels")) != 1))
        stop("formula missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    grouping <- model.extract(m, "response")
    x <- model.matrix(Terms, m)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    if(xint > 0) x <- x[, -xint, drop=FALSE]
    RET <- slda(y=grouping, X=x, ...)
    RET$terms <- Terms
    RET$call <- match.call()
    RET$contrasts <- attr(x, "contrasts")
    RET$xlevels <- xlev
    attr(RET, "na.message") <- attr(m, "na.message")
    if(!is.null(attr(m, "na.action"))) RET$na.action <- attr(m, "na.action")
    RET
}


slda.factor <- function(y, X, q=NULL, ...) {

  p <- ncol(X)
  # substract global mean 
  Xnull <- X - apply(X, 2, mean)
  if (!is.null(q)) {
    if (floor(q) != q) stop("q is not an integer")
    if (q > p) {
      q <- p
      warning("q is greater ncol(X), using q = ncol(X)")
    }
    if (q < 1) {
      q <- 1
      warning("q is less than 1, using q = 1")
    }
  }

  # this is S_0 in Kropf (2000)
  Snull <- cov(Xnull)
  ewp <- svd(solve(diag(diag(Snull), ncol = ncol(Snull)))%*%Snull)
  if (!is.complex(ewp$d)) {
    # determine q by the number of eigenvalues > 1
    if (is.null(q)) q <- sum(ewp$d > 1)
    D <- ewp$v[,1:q]
    if (q == 1) D <- as.matrix(D)
    # Xstab is still spherically distributed (Fang & Zhang, Laeuter, Kropf &
    # Glimm)!
  } else {
    D <- diag(p)
  }
  Xstab <- as.matrix(X) %*% D
  colnames(Xstab) <- paste("V", 1:ncol(Xstab), sep="")
  mylda <- lda(Xstab, grouping = y, ...)
  RET <- list(scores = D, mylda = mylda)
  class(RET) <- "slda"
  RET
}

predict.slda <- function(object, newdata, ...) {
    if(!inherits(object, "slda")) stop("object not of class slda")
    if(!is.null(Terms <- object$terms)) { #
    # formula fit (only)
        if(missing(newdata)) newdata <- model.frame(object)
        else {
            newdata <- model.frame(as.formula(delete.response(Terms)),
                                   newdata, na.action=function(x) x,
                                   xlev = object$xlevels)
        }
        x <- model.matrix(delete.response(Terms), newdata,
                          contrasts = object$contrasts)
        xint <- match("(Intercept)", colnames(x), nomatch=0)
        if(xint > 0) x <- x[, -xint, drop=FALSE]
    } else { 
      stop("object has no terms element") 
    }
    if(ncol(x) != nrow(object$scores)) stop("wrong number of variables")
#   <FIXME>: check for variable names!
#    if(length(colnames(x)) > 0 &&
#      any(colnames(x) != dimnames(object$means)[[2]]))
#         warning("Variable names in newdata do not match those in object")
#   </FIXME>
    X <- x %*% object$scores
    if (inherits(object$mylda, "lda"))
      return(predict(object$mylda, newdata=as.data.frame(X), ...))
    else
      stop(paste("Do not know how to predict from objects of class", class(object$mylda)))

}
