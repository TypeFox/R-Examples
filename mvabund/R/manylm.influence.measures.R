################################################################################
## MANYLM.INFLUENCE.MEASURES Influence measures for a manylm object.          ##
################################################################################

manylm.influence.measures <-
function (model) {

    is.influential <- function(infmat, n) {
        k   <- ncol(infmat) - 4
        if (n <= k) 
            stop("too few cases, n < k")
        absmat <- abs(infmat)
        result <- cbind(absmat[, 1:k] > 1, absmat[, k + 1] > 
            3 * sqrt(k/(n - k)), abs(1 - infmat[, k + 2]) > (3 * 
            k)/(n - k), pf(infmat[, k + 3], k, n - k) > 0.5, 
            infmat[, k + 4] > (3 * k)/n)
        dimnames(result) <- dimnames(infmat)
        result
    }

    infl <- manylm.influence(model)
    p   <- model$rank
    e   <- weighted.residuals(model)
    n.vars <- NCOL(e)
    n   <- nrow(e)
    s   <- sqrt(colSums(e^2, na.rm = TRUE)/df.residual(model))
    s   <- matrix( rep(s, each=n ) , nrow=n, ncol=n.vars )
    xxi <- chol2inv(model$qr$qr, model$qr$rank)
    si  <- infl$sigma
    h   <- infl$hat

    vn  <- variable.names(model)
    vn[vn == "(Intercept)"] <- "1_"
    
    dffits <- e * sqrt(h)/(si * (1 - h))
    cov.ratio <- (si/s)^(2 * p)/(1 - h)
    cooks.d <- ((e/(s * (1 - h)))^2 * h)/p
    
    coefs <- infl$coefficients
    sigma <- infl$sigma
    dfbs  <- infmat <- is.inf <- list()
        
    for(i in 1:n.vars) {
      dfbs[[i]] <- coefs[[i]]/outer(sigma[,i], sqrt(diag(xxi)))
      colnames(dfbs[[i]]) <- paste("dfb", abbreviate(vn), sep = ".")
      infmati <- cbind(dfbs[[i]], dffit = dffits[,i], cov.r = cov.ratio[,i],
      cook.d = cooks.d[,i], hat = h)
      infmati[is.infinite(infmati)] <- NaN
      infmat[[i]] <- infmati
      is.inf[[i]] <- is.influential(infmati, sum(h > 0))
    }
    names(infmat) <- names(is.inf)  <- colnames(e)

    ans <- list(infmat = infmat, is.inf = is.inf, call = model$call)
    class(ans) <- "infl.mvabund"
    ans

}

# setGeneric("influence.measures")   
# setMethod("influence.measures", "manylm", influence.measures.manylm)

