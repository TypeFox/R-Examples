obj.nloptr <- function (r, softinfo, epiparameters, caseinfo, c.out, constraints)
{
#
# Compute objective function value and gradients. If integrateToOne is FALSE,
# compute the integral and add it to the objective. If it's TRUE, that integral
# is done explicity elsewhere, inside integrate.to.one(). 
#

N <- epiparameters$Ndiscr
m0 <- epiparameters$m0
mN <- epiparameters$mN
p <- epiparameters$order 
Delta <- (mN-m0)/N
samplesize <- caseinfo$samplesize


if (!any (names (softinfo) == "integrateToOne") || softinfo$integrateToOne == FALSE)  {
#
# For 20 points of evaluation
#
    w <- c(0.152753387, 0.152753387, 0.149172986, 0.149172986, 0.142096109,
           0.142096109, 0.131688638, 0.131688638, 0.118194532, 0.118194532,
           0.10193012, 0.10193012, 0.083276742, 0.083276742,0.062672048,
           0.062672048, 0.04060143, 0.04060143, 0.017614007, 0.017614007)

    xivec <- c(-0.076526521,  0.076526521, -0.227785851,  0.227785851, -0.373706089,
               0.373706089, -0.510867002,  0.510867002, -0.636053681,  0.636053681,
              -0.746331906,  0.746331906, -0.839116972,  0.839116972, -0.912234428,
               0.912234428, -0.963971927,  0.963971927, -0.993128599,  0.993128599)
    
    intval <- 0
    intvalvec <- numeric ( (p+2) * N + 1)
    for (k in 1:N) {
        dex <- (Delta/2)*xivec + ((k-1)*Delta + m0 + k*Delta + m0)/2
        intval <- intval + (Delta/2) * sum(w * 
                     epispline(epiparameters, r, dex, exponentiate = TRUE))
        intvalvec <- intvalvec + (Delta/2) * grad_expepispline (epiparameters, r, dex) %*% w
    }

    f_val <- sum (c.out * r/samplesize) + intval
    f_grad <- c.out / samplesize + intvalvec

} # end "if integrate to one is "no" "
else
{    

    f_val <- sum (c.out * r/samplesize)
    f_grad <- c.out / samplesize
}

if (1 > 0) {
    if (f_val > 1e10) {
        cat ("Objective too big; returning 1e10\n")
        f_val <- 1e10
        f_grad[!is.finite(f_grad) | f_grad > 1e10] <- 1e10

    }
    if (f_val < -1e10) {
        cat ("Objective too negative; returning -1e10\n")
        f_val <- -1e10
        f_grad[!is.finite(f_grad) | f_grad < -1e10] <- -1e10
    }
}

return (list (objective = f_val, gradient=f_grad))

}
