### Program to calculate L (multi) step ahead forecast of a single sample path
### based on a glarma model fit
### Assumes that the regressor and offset values at future
### times T + 1,... T + L are specified
### by the user either as a known value
### or as forecast that the user has prepared themselves.

forecast <- function(object, ...) UseMethod("forecast")

forecast.glarma <- function(object, n.ahead = 1, newdata = 0,
                            newoffset = 0, newm = 1, ...){

    L <- NROW(newdata)
    if (L == 0) return

    T <- dim(object$X)[1]
    r <- object$r
    beta <- object$delta[1:r]
    phiLags <- object$phiLags
    thetaLags <- object$thetaLags
    e <- object$residuals
    Z <- object$W - object$eta[, 1]
    eta1toL <- numeric(L)
    Wfore1toL <- numeric(L)
    mu1toL <- numeric(L)
    Yfore1toL <- numeric(L)
    for (lead in 1:L){
        Z.lead <- 0
        if(length(phiLags)>0){
            phi <- object$delta[r + 1:length(phiLags)]
            Z.lead <- Z.lead + phi%*%(Z[T + lead - phiLags] +
                                      e[T + lead - phiLags])
        }
        if(length(thetaLags)>0){
            theta <- object$delta[r + length(phiLags) + 1:length(thetaLags)]
            Z.lead <- Z.lead + theta%*%(e[T + lead - thetaLags])
        }
        eta.lead <- newdata[lead, ]%*%beta + newoffset[lead]
        W.lead <-  eta.lead + Z.lead

        if(object$type == "Bin") {
            pi.lead <- 1/(1 + exp(-W.lead))
            mu.lead <- newm[lead]*pi.lead
            sig2.lead <- mu.lead*(1 - pi.lead)
            Y.lead <- rbinom(1, newm[lead], 1/(1 + exp(-W.lead)))
            if(object$residType == "Pearson"){
                e.lead <- (Y.lead - mu.lead)/sig2.lead^0.5
            }
            if(object$residType == "Score"){
                e.lead <- (Y.lead - mu.lead)/sig2.lead
            }
            if(object$residType == "Identity"){
                e.lead <- (Y.lead - mu.lead)
            }
        }
        if(object$type == "Poi") {
            mu.lead <- exp(W.lead)
            Y.lead <- rpois(1, exp(W.lead))
            if(object$residType == "Pearson"){
                e.lead <- (Y.lead - mu.lead)/mu.lead^0.5
            }
            if(object$residType == "Score"){
                e.lead <- (Y.lead - mu.lead)/mu.lead
            }
            if(object$residType == "Identity"){
                e.lead <- (Y.lead - mu.lead)
            }
        }
        if(object$type == "NegBin") {
            alpha <- coef(object)$NB
            mu.lead <- exp(W.lead)
            sig2.lead <-  mu.lead + mu.lead^2/(alpha^2)
            Y.lead <- rnbinom(1, mu = mu.lead, size =  alpha)
            if(object$residType == "Pearson"){
                e.lead <- (Y.lead - mu.lead)/sig2.lead^0.5
            }
            if(object$residType == "Score"){
                e.lead <- (Y.lead - mu.lead)/sig2.lead
            }
            if(object$residType == "Identity"){
                e.lead <- (Y.lead - mu.lead)
            }
        }
        eta1toL[lead] <- eta.lead
        Wfore1toL[lead] <- W.lead
        mu1toL[lead] <- mu.lead
        Yfore1toL[lead] <- Y.lead
        Z <- c(Z, Z.lead)
        e <- c(e, e.lead)
    }
    out <- list(eta = eta1toL, W = Wfore1toL,
                mu = mu1toL, Y = Yfore1toL,
                n.ahead = n.ahead, newdata = newdata,
                newoffset = newoffset, newm = newm, model = object)
    class(out) <- "glarmaForecast"
    return(out)
}

