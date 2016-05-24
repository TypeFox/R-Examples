laplace=function (logpost, mode, ...) 
{
    options(warn=-1)
    fit=optim(mode, logpost, gr = NULL, ..., hessian=TRUE,
      control=list(fnscale=-1))
    options(warn=0)
    mode=fit$par
    h=-solve(fit$hessian)
    p=length(mode)
    int = p/2 * log(2 * pi) + 0.5 * log(det(h)) +
        logpost(mode, ...)
    stuff = list(mode = mode, var = h, int = int, 
         converge=fit$convergence==0)   
    return(stuff)
}