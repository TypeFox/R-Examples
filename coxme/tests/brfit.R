library(coxme)
options(na.action=na.omit)

# This data set comes from a family study, where there were doubts about
#  the correctness of the lmekin fit. (Courtesy M deAndrade)
# It's a good test of several things: there are missing values, and the
#  data is not in family order
load('brdat.rda')

library(kinship2)
library(nlme)  
bped <- with(brdat, pedigree(id, father, mother, sex, famid=family))
kmat <- 2*kinship(bped)

fit1 <- lmekin(height ~ sex + (1|id), data= brdat, varlist=kmat)

# Compute the loglik of a fit
blik <- function(fit, kin=kmat) {
    keptrows <- as.numeric(names(fit$resid))  #not tossed away
    indx <- match(dimnames(kin)[[1]], brdat$id[keptrows], nomatch=0)
    
    k2 <- kin[indx>0, indx>0]   # retained subjects
    n <- nrow(k2)                # number who remain
    vmat <- fit$sigma^2 * (diag(n) + unlist(VarCorr(fit))*k2)
    
    r2 <- fit$resid[indx]   #residuals, in kmat order
    smat <- gchol(as(vmat, "bdsmatrix"))  #cholesky decomposition
    
    rsum <- sum(r2 * solve(smat, r2, full=TRUE))
    csum <- n*log(2*pi) + sum(log(diag(smat)))

    # Check using Matrix routines
    rsum2  <- sum(r2 * solve(vmat, r2))
    csum2  <- n*log(2*pi) + 2*sum(log(diag(chol(vmat))))
    if (!all.equal(c(rsum, csum), c(rsum2, csum2))) cat("Matrix error")
   -(rsum + csum)/2
}

# Test on a model I can verify with lme
tfit1 <- lmekin(height ~ sex + (1|family), brdat)
tfit2 <- lme(height ~ sex, random= ~1|family, brdat, method="ML",
             na.action=na.omit)
tkmat <- with(brdat, bdsBlock(id, family))
blik(tfit1, as(tkmat, 'dsCMatrix'))
                                 
