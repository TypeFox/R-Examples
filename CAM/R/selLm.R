selLm <-
function(X,pars = list(cutOffPVal = 0.001),output = FALSE,k)
{
    if(output)
    {
        cat("Performing variable selection for variable", k, ": \n")
    }
    result <- list()
    p <- dim(as.matrix(X))
    if(p[2] > 1)
    {
        selVec <- rep(FALSE, p[2])
        x <- X[,-k]; y <- X[,k] 
        modfitLm <- lm(y ~ x)
        pValVec <- coef(summary(modfitLm))[-1,4] 
        selVec[-k] <- as.vector((pValVec < pars$cutOffPVal))
    } else
    {
        selVec <- list()
    }
    return(selVec)
}
