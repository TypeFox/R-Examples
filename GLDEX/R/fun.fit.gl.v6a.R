`fun.fit.gl.v6a` <-
function (a, b, c, d, data, param,len=1000,type=7) 
{
    init.sol <- c(a, b, c, d)
    optim.result <- optim(init.sol, optim.fun6, data = data, 
        param = param, len=len,type=type,control = list(maxit = 2e+05))
    return(optim.result$par)
}

