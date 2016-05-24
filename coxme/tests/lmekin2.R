library(coxme)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)
#
# A set of tests using user-defined matrices
#  The lung data is handy, though models of time that ignore status aren't
#  very sensible.
#
fit1 <- lmekin(time ~ age + (1|ph.ecog), lung, vfixed=.5)

# Solve using the extended normal equations
#  The loglik is sum(y - hat y)^2 + b'P b/2 where the penalty matrix
#  is all zeros for the fixed effects, and A-inverse for the random
#  where b is the coefficient vector and A is the variance of the random effects
#  In our case A-inverse/2 = identity matrix.  
# If the variance is known this equation is easy to solve.
efit <- function(p, pmat) {
    tlung <- na.omit(lung[,c("time", "age", "ph.ecog")])
    cmat <- cbind(0,0, rbind(0,0, pmat))
    x <- with(tlung, cbind(1, age, ph.ecog==0, ph.ecog==1,
                          ph.ecog==2, ph.ecog==3))
    as.vector(solve(t(x) %*% x + p*cmat, t(x) %*% tlung$time))
}
efit1 <- efit(2, diag(4))

aeq(fixef(fit1), efit1[1:2])
aeq(unlist(ranef(fit1)), efit1[3:6])

# Now put in a more complex pmat (linear trend constraint)
# Since coxme is using the inverse variance as the penalty, we
#  need to give the inverse penalty as the "variance".
# (This is an lmekin test with non-zero off diagonal elements 
#  in the variance for the random effect, which is key to working
#  with kinship matrices.)
pmat <- matrix(c(2,-1,0,0, -1,2,-1,0, 0, -1, 2, -1, 0,0,-1,2),4)
dimnames(pmat) <- list(0:3, 0:3)

fit2 <- lmekin(time ~ age + (1|ph.ecog), lung, vfixed=.5,
               varlist= coxmeMlist(solve(pmat), rescale=F))
efit2 <- efit(2, pmat)
aeq(fixef(fit2), efit2[1:2])
aeq(unlist(ranef(fit2)), efit2[3:6])

