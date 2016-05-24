confusion.fda <-
function (object, data, ...) 
{
    if (missing(data)) 
        return(object$confusion)
    Terms <- terms(object)
    attr(Terms, "intercept") <- 0
    m <- model.frame(Terms, data)
    x <- model.matrix(Terms, m)
    g <- model.extract(m, "response")
    confusion.default(predict(object, x, ...), g)
}

