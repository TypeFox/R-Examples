plsrPIsboot = function (X, Y, Xnew, ncomp, B, alpha, weight, beta = .1) 
{
    pred = t(X)
    resp = t(Y)
    n = dim(pred)[1]
    p = dim(pred)[2]
    mpred = apply(pred, 2, mean)
    mresp = apply(resp, 2, mean)
    decenpred = scale(pred, scale = FALSE)
    decenresp = scale(resp, scale = FALSE)
    resi = matrix(NA, n, p)
    resib = resiY = matrix(NA, n, p)
    predictvalue = matrix(NA, B, p)
    resibb = matrix(NA, B, p)
    resi2 = matrix(, n, p)
    if (weight == FALSE){
        for (i in 1:p) {
             mod = plsr(pred, resp[, i], Xtest = pred, ncomp = ncomp, 
                        type = "simpls", unit.weights = FALSE)
             resi[, i] = (mod$Ypred[, , ] - resp[, i]) / sqrt(1 - (n - p - 1) / n)
        }
        for (j in 1:B) {
             for (i in 1:p) {
                  resib[, i] = sample(resi[, i], size = n, replace = TRUE)
                  resiY[, i] = resp[, i] + resib[, i]
             }
             mod = plsr(pred, resiY, Xtest = t(Xnew), ncomp = ncomp, 
                   type = "simpls", unit.weights = FALSE)
             predictvalue[j, ] = mod$Ypred[, , ]
         }
         for (i in 1:p) {
              mod = plsr(pred, resp[, i], Xtest = pred, ncomp = ncomp, 
                    type = "simpls", unit.weights = FALSE)
              resi2[, i] = (mod$Ypred[, , ] - resp[, i]) / sqrt(1 - (n - p - 1) / n)
         }
    }
    if (weight == TRUE){
        w = matrix(,n,1)
        for(i in 1:n){
            w[i,] = beta * (1 - beta)^(i - 1)
        }
        weight = diag(rev(w))
        newpred = newresp = matrix(,n,p)
        for(i in 1:p){
            newpred[,i] = (weight %*% decenpred)[,i] + mpred[i]
            newresp[,i] = (weight %*% decenresp)[,i] + mresp[i]
        }
        for (i in 1:p) {
            mod = plsr(newpred, newresp[,i], Xtest = newpred, ncomp = ncomp,
                        type = "simpls", unit.weights = FALSE)
            resi[, i] = (mod$Ypred[, , ] - newresp[,i])
        }
        for(j in 1:B){
            for (i in 1:p) {
                 resib[, i] = sample(resi[, i], size = n, replace = TRUE)
                 resiY[, i] = newresp[, i] + resib[, i]
            }
            mod = plsr(newpred, resiY, Xtest = t(Xnew), ncomp = ncomp,
                       type = "simpls", unit.weights = FALSE)
            predictvalue[j, ] = mod$Ypred[, , ]
        }
        resibb = matrix(NA, B, p)
        resi2 = matrix(, n, p)
        for (i in 1:p) {
             mod = plsr(pred, resp[,i], Xtest = pred, ncomp = ncomp, 
                        type = "simpls", unit.weights = FALSE)
             resi2[, i] = (mod$Ypred[, , ] - resp[,i]) / sqrt(1 - (n - p - 1) / n)
        }
    }
    for (i in 1:p) {
         resibb[, i] = sample(resi2[, i], size = B, replace = TRUE)
    }
    q1 = q2 = matrix(, p, 1)
    for (i in 1:p) {
         q1[i, ] = quantile(predictvalue[, i] + resibb[, i], alpha/2, 
                   na.rm = TRUE)
    }
    for (i in 1:p) {
         q2[i, ] = quantile(predictvalue[, i] + resibb[, i], 1 - alpha/2, 
                   na.rm = TRUE)
    }
    return(list(bootsamp = t(predictvalue + resibb), lb = q1, ub = q2))
}
