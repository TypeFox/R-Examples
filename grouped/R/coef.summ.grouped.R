"coef.summ.grouped" <-
function(object, ...){
    if(!inherits(object, "summ.grouped"))
        stop("Use only with 'summ.grouped' objects.\n")
    round(object$coef, 3)
}

