pri.mle <-
function(nrisk, ndeath, age = c(0, 1, seq(5, 100, 5)), lo = c(1e-08, 
    1e-07, 1e-07, 1e-07, 1e-07, 15, 1e-07, 1), hi = c(1, 1, 1, 
    0.5, 15, 55, 0.1, 1.5), senum = 15, theta.test = c(0.06008, 
    0.31087, 0.34431, 0.00698, 1.98569, 26.71071, 0.00022, 1.088), 
    opt.meth = "Nelder-Mead") 
{
    x.fit <- age
    mp8.mle <- function(theta) {
        p.hat <- mod8p(theta = theta, x = x.fit)
        ll <- ll.binom(x = ndeath, n = nrisk, p = p.hat)
        return(ll)
    }
    if (opt.meth != "L-BFGS-B") {
        opt.out <- optim(par = theta.test, fn = mp8.mle, method = opt.meth, 
            hessian = TRUE, control = list(fnscale = -1, maxit = 10000))
        mle <- opt.out$par
    }
    if (opt.meth == "L-BFGS-B") {
        opt.out <- optim(par = theta.test, fn = mp8.mle, method = "L-BFGS-B", 
            hessian = TRUE, lower = lo, upper = hi, control = list(fnscale = -1, 
                maxit = 10000))
        mle <- opt.out$par
    }
    hess <- opt.out$hessian
    var.cov <- solve(-hess)
    d.vc <- diag(var.cov)
    se.out <- sqrt(abs(d.vc))
    mle1.low <- NULL
    if (mle[1] - senum * se.out[1] < lo[1]) {
        mle1.low <- lo[1]
    }
    else mle1.low <- mle[1] - senum * se.out[1]
    mle2.low <- NULL
    if (mle[2] - senum * se.out[2] < lo[2]) {
        mle2.low <- lo[2]
    }
    else mle2.low <- mle[2] - senum * se.out[2]
    mle3.low <- NULL
    if (mle[3] - senum * se.out[3] < lo[3]) {
        mle3.low <- lo[3]
    }
    else mle3.low <- mle[3] - senum * se.out[3]
    mle4.low <- NULL
    if (mle[4] - senum * se.out[4] < lo[4]) {
        mle4.low <- lo[4]
    }
    else mle4.low <- mle[4] - senum * se.out[4]
    mle5.low <- NULL
    if (mle[5] - senum * se.out[5] < lo[5]) {
        mle5.low <- lo[5]
    }
    else mle5.low <- mle[5] - senum * se.out[5]
    mle6.low <- NULL
    if (mle[6] - senum * se.out[6] < lo[6]) {
        mle6.low <- lo[6]
    }
    else mle6.low <- mle[6] - senum * se.out[6]
    mle7.low <- NULL
    if (mle[7] - senum * se.out[7] < lo[7]) {
        mle7.low <- lo[7]
    }
    else mle7.low <- mle[7] - senum * se.out[7]
    mle8.low <- NULL
    if (mle[8] - senum * se.out[8] < lo[8]) {
        mle8.low <- lo[8]
    }
    else mle8.low <- mle[8] - senum * se.out[8]
    mle1.hi <- NULL
    if (mle[1] + senum * se.out[1] > hi[1]) {
        mle1.hi <- hi[1]
    }
    else mle1.hi <- mle[1] + senum * se.out[1]
    mle2.hi <- NULL
    if (mle[2] + senum * se.out[2] > hi[2]) {
        mle2.hi <- hi[2]
    }
    else mle2.hi <- mle[2] + senum * se.out[2]
    mle3.hi <- NULL
    if (mle[3] + senum * se.out[3] > hi[3]) {
        mle3.hi <- hi[3]
    }
    else mle3.hi <- mle[3] + senum * se.out[3]
    mle4.hi <- NULL
    if (mle[4] + senum * se.out[4] > hi[4]) {
        mle4.hi <- hi[4]
    }
    else mle4.hi <- mle[4] + senum * se.out[4]
    mle5.hi <- NULL
    if (mle[5] + senum * se.out[5] > hi[5]) {
        mle5.hi <- hi[5]
    }
    else mle5.hi <- mle[5] + senum * se.out[5]
    mle6.hi <- NULL
    if (mle[6] + senum * se.out[6] > hi[6]) {
        mle6.hi <- hi[6]
    }
    else mle6.hi <- mle[6] + senum * se.out[6]
    mle7.hi <- NULL
    if (mle[7] + senum * se.out[7] > hi[7]) {
        mle7.hi <- hi[7]
    }
    else mle7.hi <- mle[7] + senum * se.out[7]
    mle8.hi <- NULL
    if (mle[8] + senum * se.out[8] > hi[8]) {
        mle8.hi <- hi[8]
    }
    else mle8.hi <- mle[8] + senum * se.out[8]

    pri.lo <- c(mle1.low, mle2.low, mle3.low, mle4.low, mle5.low, 
        mle6.low, mle7.low, mle8.low)
    pri.hi <- c(mle1.hi, mle2.hi, mle3.hi, mle4.hi, mle5.hi, 
        mle6.hi, mle7.hi, mle8.hi)
    q0 <- prior.form(pri.lo = pri.lo, pri.hi = pri.hi)
    return(list(q0 = q0, mle = mle, se.out = se.out, pri.lo=pri.lo, pri.hi=pri.hi))
}

