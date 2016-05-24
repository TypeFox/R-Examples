# $Id: ipredknn.R,v 1.5 2005/06/29 08:50:28 hothorn Exp $

# k-NN compatible with the fit(formula) - predict(object) framework

ipredknn <- function(formula, data, subset, na.action, k=5, ...) {
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
    m$k <- NULL
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")   
    y <- model.extract(m, "response")
    x <- model.matrix(Terms, m)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    if(xint > 0) x <- x[, -xint, drop=FALSE]
    RET <- list(learn=list(y=y, X=x))
    RET$k <- k
    RET$terms <- Terms
    RET$call <- match.call()
    RET$contrasts <- attr(x, "contrasts")
    RET$xlevels <- xlev
    attr(RET, "na.message") <- attr(m, "na.message")
    if(!is.null(attr(m, "na.action"))) RET$na.action <- attr(m, "na.action")
    class(RET) <- "ipredknn"
    RET
}

predict.ipredknn <- function(object, newdata, type=c("prob", "class"), ...) {
    type <- match.arg(type)
    if(!inherits(object, "ipredknn")) stop("object not of class ipredknn")
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
#    <FIXME>: check for variable names
#    if(length(colnames(x)) > 0 &&
#      any(colnames(x) != dimnames(object$means)[[2]]))
#         warning("Variable names in newdata do not match those in object")
#   </FIXME>
  RET <- knn(object$learn$X, x, 
             object$learn$y, k=object$k, prob=TRUE)
  if (type=="prob") return(attr(RET, "prob"))
  else return(RET)
}
