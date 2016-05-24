"mr.test" <- function(object1, object2, object, x, var.equal = TRUE, component = 1)
{
    dnFct <- deriv(~c+(d-c)/(1+(x/e)^b),c("b","c","d","e"), function(b,c,d,e,x){}, hessian = TRUE)
#    dnFct <- deriv(~d*exp(-e*x),c("d","e"), function(d,e,x){}, hessian = TRUE)  # exponential model
#    dnFct <- deriv(~d*exp(-e/x),c("d","e"), function(d,e,x){}, hessian = TRUE)  # exponential model
    daFct <- deriv(~c+(d-c)*exp(-exp(b*(log(x)-log(e)))),c("b","c","d","e"), function(b,c,d,e,x){}, 
    hessian = TRUE)

    beta1 <- coef(object1)
    beta2 <- coef(object2)
    beta <- coef(object)
    diffVec <- beta2 - beta
    
    lenx <- length(x)    
#    if (identical(var.equal, FALSE))
#    {
#        res1 <- residuals(object1)
#        res2 <- residuals(object2)
#    } else {
#        res1 <- rep(sqrt(summary(object1)$resVar), lenx) 
#        res2 <- res1 
#    }
#    res1 <- rep(1, lenx)
#    res2 <- res1
    
    Df1 <- attr(dnFct(beta1[1], beta1[2], beta1[3], beta1[4], x), "gradient")  # * res1
#    Df1 <- attr(dnFct(beta1[1], beta1[2], x), "gradient")  # exponential model
    
#    D2g1 <- attr(dnFct(theta1[1],theta1[2],theta1[3],theta1[4], x), "hessian")  # not needed
    Df1Df1 <- matrix(apply(apply(Df1, 1, function(x){x%*%t(x)}), 1, mean), 4, 4)
#    print(Df1)
#    print(Df1Df1)

    Df2 <- attr(daFct(beta2[1], beta2[2], beta2[3], beta2[4], x), "gradient")  # * res1
    Df2Df2 <- matrix(apply(apply(Df2, 1, function(x){x%*%t(x)}), 1, mean), 4, 4)
#    print(Df2)    
#    print(Df2Df2)

    ## Product of first derivatives from the two models
#    lenx <- length(x)
    Df1Df2.0 <- matrix(0, 4 * 4, lenx)
    for (i in 1:lenx)
    {
        Df1Df2.0[, i] <- as.vector(Df1[i, ]%*%t(Df2[i, ]))   
    }
    Df1Df2 <- matrix(apply(Df1Df2.0, 1, mean), 4, 4)

    ## Second derivatives for the alternative model
    fitDiff <- fitted(object2) - fitted(object1)
#    print(fitDiff)
    D2f2.0 <- attr(daFct(beta2[1], beta2[2], beta2[3], beta2[4], x), "hessian") 
    D2f2 <- matrix(0, 4, 4)
    for (i in 1:4)
    {
        D2f2[i, ] <- apply(D2f2.0[, , i] * fitDiff, 2, mean) 
    }
#    print(Df2Df2)
#    print(D2f2)

#    ## Using derivatives without residuals multiplied
#    Df2.2 <- attr(daFct(beta2[1], beta2[2], beta2[3], beta2[4], x), "gradient")
#    Df2Df2.2 <- matrix(apply(apply(Df2.2, 1, function(x){x%*%t(x)}), 1, mean), 4, 4)
    H <- Df2Df2 + D2f2
    Hinv <- solve(H)
    cov12 <- Df1Df2
    cov21 <- t(Df1Df2)
    var1 <- Df1Df1
    var2 <- Df2Df2

    Wfit <- Hinv %*% cov21 %*% solve(var1) %*% cov12 %*% Hinv
#    Wfit <- Hinv %*% cov21 %*% ginv(var1) %*% cov12 %*% Hinv  # exponential model
    Wobs <- Hinv %*% var2 %*% Hinv
    
#    varEst <- (Wobs - Wfit)/(lenx*summary(object1)$resVar)
    varEst <- (Wobs - Wfit) * (summary(object1)$resVar) / lenx
#    varEst <- (Wobs - Wfit)  / lenx

#    print(diffVec/sqrt(diag(varEst)))
#    print(varEst)
#    print(solve(varEst[1:4,1:4]))
     if (identical(var.equal, FALSE))
     {
         Hinv <- Hinv / lenx
         varEst0 <- - Hinv %*% t(Df2) + Hinv %*% cov21 %*% solve(var1) %*% t(Df1)
         
         p2 <- length(beta2)
         veMat <- matrix(NA, lenx, p2^2)
         res1 <- residuals(object1)
         for (i in 1:lenx)
         {
#             veMat[i, ] <- as.vector(outer(varEst0[, i], varEst0[, i])) * summary(object1)$resVar * lenx  # (res1[i]^2)
             veMat[i, ] <- as.vector(outer(varEst0[, i], varEst0[, i])) * (res1[i]^2) * lenx
         }
         varEst <- matrix(apply(veMat, 2, mean), p2, p2)
     }
    
#    sInd <- c(1, 2, 3, 4)
#    chiSq <- diffVec[sInd]%*%solve(varEst[sInd, sInd])%*%diffVec[sInd]
#    chiSq <- diffVec[sInd]%*%ginv(varEst[sInd, sInd])%*%diffVec[sInd]
#    pVal <- 1 - pchisq(chiSq, length(sInd))
    
    
    stdErr <- sqrt(varEst[component, component])
    diffVal <- diffVec[component]
    chiSq <- diffVal / stdErr
    pVal <- 2 * (1 - pnorm(abs(chiSq)))

    retVec <- c(chiSq, pVal, diffVal ,stdErr)
    names(retVec) <- c("Statistic", "p-value", "Difference", "SE")
    return(retVec)
}


