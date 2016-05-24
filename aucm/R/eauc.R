eauc = function (formula, dat, t0=NULL, t1=NULL){
    
    fit.rlogit=rlogit(formula, dat)
    if (fit.rlogit$convergence) {
        if (length(coef(fit.rlogit))!=3) stop ("formula can only have two covariates: "%+%formula)
        a = coef(fit.rlogit)[3]/coef(fit.rlogit)[2]
    } else {
        fit.glm = glm(formula, dat, family="binomial")
        if (length(coef(fit.glm))!=3) stop ("formula can only have two covariates: "%+%formula)
        a = coef(fit.glm)[3]/coef(fit.glm)[2]
    }                
    beta = cbind(1, seq(a-5, a+5, length=1000))
    
    grid.auc (formula, dat, beta, loss=TRUE, t0=t0, t1=t1)
}
