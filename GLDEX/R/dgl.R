dgl <-
function (x, lambda1 = 0, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
    param = "fmkl", inverse.eps =1e-08, max.iterations = 500) 
{
    lambdas <- .gl.parameter.tidy(lambda1, lambda2, lambda3, 
        lambda4, param)
    if (!fun.check.gld(lambdas, param = param)) {
         return(rep(NA,length(x)))
    }

    u <- rep(0,length(x))
    result <- .C("dgl", param, as.double(lambdas[1]), as.double(lambdas[2]), 
    as.double(lambdas[3]), as.double(lambdas[4]), 
    inverse.eps, as.integer(max.iterations), as.double(x), as.double(u), 
    as.integer(length(x)), as.double((.Machine$double.eps)))

    result[[9]]
}
