# $Id: bagging.R,v 1.19 2005/06/29 08:50:28 hothorn Exp $

bagging <- function(formula, data, ...) UseMethod("bagging", data)

bagging.default <- function(formula, data, ...)
  stop(paste("Do not know how to handle objects of class", class(data)))

bagging.data.frame <-
function(formula, data, subset, na.action=na.rpart, ...)
{
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
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    # just extract the data.frame, no handling of contrasts or NA's here.
    # this is done by rpart or the user supplied methods

    DATA <- list(y = mf[,response], X = mf[,-response, drop = FALSE]) 
    names(DATA) <- c("y", "X")
    y <- do.call("ipredbagg", c(DATA, list(...)))
    y$call <- cl
    return(y)
}

