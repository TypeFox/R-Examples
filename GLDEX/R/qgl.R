qgl <-
function (p, lambda1, lambda2 = NULL, lambda3 = NULL, lambda4 = NULL, 
    param = "fmkl") 
{
    lambdas <- .gl.parameter.tidy(lambda1, lambda2, lambda3, 
        lambda4, param)
    if (!fun.check.gld(lambdas, param = param)) {
       return(rep(NA,length(p)))
    }
    result <- switch(param, FKML= , fkml = , freimer = , frm = , FMKL = , fmkl = .qgl.fmkl(p, 
        lambdas), ramberg = , ram = , RS = , rs = .qgl.rs(p, 
        lambdas), stop("Error: Parameterisation must be fmkl or rs"))
    result
}
