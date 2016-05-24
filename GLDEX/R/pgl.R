pgl <-
function (q, lambda1 = 0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
    param = "fmkl", inverse.eps = 1e-08, max.iterations = 500) 
{
# tidy up the parameters
    lambdas <- .gl.parameter.tidy(lambda1, lambda2, lambda3, 
        lambda4, param)

    if (!fun.check.gld(lambdas, param = param)) {
       return(rep(NA,length(q)))
    }

    if (tolower(param)=="fkml" | tolower(param)=="fmkl" ){

    param<-"fmkl"

    }

    u <- rep(0,length(q))
    result <- .C("pgl", param,  as.double(lambdas[1]), as.double(lambdas[2]), 
    as.double(lambdas[3]), as.double(lambdas[4]), 
    inverse.eps, as.integer(max.iterations), as.double(q), as.double(u), as.integer(length(q)),
    as.double(.Machine$double.eps))

    if (!(is.numeric(result[[2]]))) {
        stop("Values for quantiles outside range. This shouldn't happen")
    }
    else {
    # u is probably not the right length, might have trailing 0s because it's
    # actually the same size as q but might not contain all of the elements?
        u <- result[[9]]
    }
    #print(u)

   u
}
