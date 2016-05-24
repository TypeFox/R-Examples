fun.check.gld <-
function(lambda1, lambda2, lambda3,lambda4,param){

    lambdas <- .gl.parameter.tidy(lambda1, lambda2, lambda3, 
        lambda4, param)    
    param <- switch(param, FKML= , fkml = , freimer = , frm = , FMKL = , fmkl = {
        ret <- .C("check_gld",as.double(lambdas[1]),as.double(lambdas[2]),
        as.double(lambdas[3]),as.double(lambdas[4]),as.integer(1),as.integer(0)) }
    , ramberg = , ram = , RS = , rs = {
        ret <- .C("check_gld",as.double(lambdas[1]),as.double(lambdas[2]),
        as.double(lambdas[3]),as.double(lambdas[4]),as.integer(2),as.integer(0)) })
    return(as.logical(ret[[6]]))
}
