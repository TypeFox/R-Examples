### 'transfer' class
setClass("transfer", representation(para = "numeric",
                                    Y = "numeric")
         )


### Constructor
transfer = function (N = 10000, d = 0.5, deffect = 1, lambda = 120,
                     Q = 0.05, l0 = 0.8, u0 = 0.9, lstar0 = 0.1, ustar0 = 0.15,
                     lj = 0.45, uj = 0.7, lstarj = 0.05, ustarj=0.1, lR=0.5,
                     uR = 0.7, t = 1.5, r = 0.5) {

    Y = rep(NA,N)
    para = c(N, d, deffect, lambda, Q, l0, u0, lstar0, ustar0,
             lj, uj, lstarj, ustarj, lR, uR, t, r)
    names(para) = c("N", "d", "deffect", "lambda", "Q", "l0", "u0",
                    "lstar0", "ustar0", "lj", "uj", "lstarj", "ustarj",
                    "lR", "uR", "t", "r")

    for (i in 1:N) {
        di = rgamma(1, d)
        lambdai = rnorm(1, lambda, lambda/2)
        ##weight = rnorm(1, 1, 0.5)
        ##lambdai = lambda*weight
        ti = rnbinom(1, t, r)

        if(deffect == 1){
            x0 = rpois(1, abs(lambdai * exp(1-(di/d))))
        }
        else if (deffect == 0){
            x0 = rpois(1, abs(lambdai))
        }
        else {
            cat("Error in transfer(N, d, deffect, ...) : invalid distance effect", "\n")
        }

        q0 = rbinom(1, x0, Q)
        xj = max(x0 - q0 - rbinom(1, x0, runif(1, l0, u0)), 0)
        qj = max(q0 - rbinom(1, q0, runif(1, lstar0, ustar0)), 0)

        if (ti > 1) {
            for (j in 2:ti) {
                qj = max(qj - rbinom(1, qj, runif(1, lstarj, ustarj)), 0)
                xj = max(xj- rbinom(1, xj, runif(1, lj, uj)),0)
            }
        }
        yi = xj + qj
        Y[i] = max(yi - rbinom(1, yi, 1-runif(1, lR, uR)),0)
    }
    invisible(Y)
    new("transfer", para = para, Y = Y)
}

values = function(object) object@Y
para = function(object) object@para
