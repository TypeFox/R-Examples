estimates.bma <-
function (object, exact = FALSE, order.by.pip = TRUE, include.constant = FALSE, 
    incl.possign = TRUE, std.coefs = FALSE, condi.coef = FALSE) 
{
    bmao = object
    rm(object)
    if (!is.bma(bmao)) {
        stop("you need to provide a BMA object")
        return()
    }
    if (exact) {
        if (bmao$topmod$nbmodels == 0) 
            stop("exact=TRUE needs at least one 'top model': Run estimation again and set nmodel>0")
    }
    bmaest = .post.estimates(bmao$info$b1mo, bmao$info$b2mo, 
        bmao$info$cumsumweights, bmao$info$inccount, bmao$topmod, 
        bmao$X.data, bmao$reg.names, bmao$info$pos.sign, exact, 
        order.by.pip, include.constant, incl.possign, std.coefs, 
        condi.coef)
    return(bmaest)
}
