activeSetLogCon <- function(x, xgrid = NULL, print = FALSE, w = NA){

    prec <- 1e-10
    xn <- sort(x)
    
    if ((!identical(xgrid, NULL) & (!identical(w, NA)))){stop("If w != NA then xgrid must be NULL!\n")}  
    
    if (identical(w, NA)){
        tmp <- preProcess(x, xgrid = xgrid)
        x <- tmp$x
        w <- tmp$w
        sig <- tmp$sig
        }
        
    if (!identical(w, NA)){
        tmp <- cbind(x, w)
        tmp <- tmp[order(x), ]
        x <- tmp[, 1]
        w <- tmp[, 2]
        
        est.m <- sum(w * x)
        est.sd <- sum(w * (x - est.m) ^ 2)
        est.sd <- sqrt(est.sd * length(x) / (length(x) - 1))
        sig <- est.sd
    }

    n <- length(x)    
    phi <- LocalNormalize(x, 1:n * 0)
    IsKnot <- 1:n * 0
    IsKnot[c(1, n)] <- 1
    res1 <- LocalMLE(x, w, IsKnot, phi, prec)
    phi <- res1$phi
    L <- res1$L
    conv <- res1$conv
    H <- res1$H
    iter1 <- 1
    while ((iter1 < 500) & (max(H) > prec * mean(abs(H)))){
        IsKnot_old <- IsKnot
        iter1 <- iter1 + 1
        tmp <- max(H)
        k <- (1:n) * (H == tmp)
        k <- min(k[k > 0])
        IsKnot[k] <- 1
        res2 <- LocalMLE(x, w, IsKnot, phi, prec)
        phi_new <- res2$phi
        L <- res2$L
        conv_new <- res2$conv
        H <- res2$H
        while ((max(conv_new) > prec * max(abs(conv_new)))){
            JJ <- (1:n) * (conv_new > 0)
            JJ <- JJ[JJ > 0]
            tmp <- conv[JJ]/(conv[JJ] - conv_new[JJ])
            lambda <- min(tmp)
            KK <- (1:length(JJ)) * (tmp == lambda)
            KK <- KK[KK > 0]
            IsKnot[JJ[KK]] <- 0
            phi <- (1 - lambda) * phi + lambda * phi_new
            conv <- pmin(c(LocalConvexity(x, phi), 0))
            res3 <- LocalMLE(x, w, IsKnot, phi, prec)
            phi_new <- res3$phi
            L <- res3$L
            conv_new <- res3$conv
            H <- res3$H
        }
        phi <- phi_new
        conv <- conv_new
        if (sum(IsKnot != IsKnot_old) == 0){break}
        if (print == TRUE){
            print(paste("iter1 = ", iter1 - 1, " / L = ", round(L, 4), " / max(H) = ", round(max(H), 4), " / #knots = ", 
            sum(IsKnot), sep = ""))
        }
    }
    Fhat <- LocalF(x, phi)
    res <- list(xn = xn, x = x, w = w, phi = as.vector(phi), IsKnot = IsKnot, L = L, Fhat = as.vector(Fhat), 
        H = as.vector(H), n = length(xn), m = n, knots = x[IsKnot == 1], mode = x[phi == max(phi)], sig = sig)
    return(res)
}
