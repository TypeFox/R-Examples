"fun.fit.gl.v3a" <-
function(a, b, c, d, data, param)
{

init.sol <- c(a, b, c, d)
optim.result <- optim(init.sol, optim.fun3, data = data, param = param, control = list(maxit = 200000))
return(optim.result$par)
}

