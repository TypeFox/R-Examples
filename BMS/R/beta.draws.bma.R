beta.draws.bma <-
function (bmao, stdev = FALSE) 
{
    if (!is.bma(bmao)) {
        stop("you need to provide a BMA object")
        return()
    }
    resmat = .post.beta.draws(bmao$topmod, bmao$reg.names, FALSE)
    if (stdev) {
        mom2 = .post.beta.draws(bmao$topmod, bmao$reg.names, 
            TRUE)
        resmat = sqrt(mom2 - resmat^2)
    }
    return(resmat)
}
