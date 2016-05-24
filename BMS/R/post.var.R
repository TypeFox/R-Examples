post.var <-
function (object, exact = FALSE) 
{
    if (!(is.bma(object) | is(object, "lm"))) 
        stop("Required input is an object of class 'bma' or 'lm'/'zlm'.")
    od = deviance(object, exact = exact)
    oy = model.frame(object)[, 1, drop = TRUE]
    ret = od/length(oy)
    attr(ret, "nobs") = length(oy)
    return(ret)
}
