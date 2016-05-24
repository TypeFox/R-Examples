`fun.fit.gl.v.lma` <-
function (a, b, c, d, data, param) 
{
    init.sol <- c(a, b, c, d)
    optim.result <- optim(init.sol, optim.fun.lm, data = data, 
        param = param, control = list(maxit = 2e+05))
    return(optim.result$par)
}

