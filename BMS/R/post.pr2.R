post.pr2 <-
function (object, exact = FALSE) 
{
    if (!(is.bma(object) | is(object, "lm"))) 
        stop("Required input is an object of class 'bma' or 'lm'/'zlm'.")
    od = deviance(object, exact = exact)
    oy = model.frame(object)[, 1, drop = TRUE]
    return(1 - (od/crossprod(oy - mean(oy))[[1]]))
}
