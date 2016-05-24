gl.check.lambda.alt1 <-
function (l1, l2 = NULL, l3 = NULL, l4 = NULL, param = "fmkl", 
    vect = FALSE) 
{
    if (vect == TRUE && length(l1) == 1) {
        stop("You need to provide a vector of parameters")
    }
    lambdas <- .gl.parameter.tidy(lambda1 = l1, lambda2 = l2, 
        lambda3 = l3, lambda4 = l4, param = param)
    l4 = lambdas[4]
    l3 = lambdas[3]
    l2 = lambdas[2]
    l1 = lambdas[1]
    param <- switch(param, FKML= , fkml = , freimer = , frm = , FMKL = , fmkl = {
        ret <- l2 > 0
    }, ramberg = , ram = , RS = , rs = {
        ret <- rep(0, length(l1))
        con1 <- (l3 < -1) * (l4 > 1)
        con2 <- (l3 > 1) * (l4 < -1)
        con3 <- (l4 > 1) * (l3 > -1) * (l3 < 0) * (((1 - l3)^(1 - 
            l3) * (l4 - 1)^(l4 - 1))/((l4 - l3)^(l4 - l3)) < 
            -l3/l4)
        con4 <- (l3 < 0) * (l4 <= 0)
        con5 <- (l3 == 0) * (l4 < 0)
        con6 <- (l3 > 1) * (l4 > -1) * (l4 < 0) * (((1 - l4)^(1 - 
            l4) * (l3 - 1)^(l3 - 1))/((l3 - l4)^(l3 - l4)) < 
            -l4/l3)
        con6[which.na(con6)] <- 0
        con3[which.na(con3)] <- 0
        ret[(l2 < 0)] <- ((con1 + con2 + con3 + con4 + con5 + 
            con6) > 0)[(l2 < 0)]
        con7 <- (l3 > 0) * (l4 >= 0)
        con8 <- (l3 == 0) * (l4 > 0)
        ret[l2 > 0] <- (con7 + con8)[(l2 > 0)]
    })
    ret <- as.logical(ret * ((is.finite(l1) * is.finite(l2) * 
        is.finite(l3) * is.finite(l4)) == 1))
    return(ret)
}
