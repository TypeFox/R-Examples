bootreg <-
function(formula, data, nboot=1000){
    obj <- lm(formula, data=data)
    res <- resid(obj)
    hat <- fitted(obj)
    m <- length(coef(obj))
    bootmat <- matrix(0, nrow=nboot, ncol=m)
    rhs <- deparse(delete.response(terms(formula)))
    bootform <- as.formula(paste("yBoot", rhs))
    for(i in 1:nboot){
        resboot <- sample(res, replace=T)
        yBoot <- hat+resboot
        bootmat[i,] <- coef(lm(bootform, data=data))
    }
    colnames(bootmat) <- names(coef(obj))
    bootmat
}
