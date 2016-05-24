simreg <-
function(formula, data, nsim=1000){
        obj <- lm(formula, data=data)
        res <- resid(obj)
        sd <- summary(obj)$sigma
        hat <- fitted(obj)
        m <- length(coef(obj))
        regmat <- matrix(0, nrow=nsim, ncol=m)
        nvar <- length(hat)
        simform <- as.formula(paste("ySim", deparse(delete.response(terms(formula)))))
        for(i in 1:nsim){
            resSim <- rnorm(nvar, 0, sd)
            ySim <- hat+resSim
            regmat[i,] <- coef(lm(simform, data=data))
        }
        colnames(regmat) <- names(coef(obj))
        regmat
    }
