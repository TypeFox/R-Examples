predict1 <- 
# predict wrapper to deal with awkward predict methods
function (object, ...)
{
    type <- if (inherits(object, "nnet"))
        "class"
    else if (inherits(object, "rpart")) 
        "vector"
    else "response"    
    predict(object, ..., n.trees = object$n.trees, type = type)
}