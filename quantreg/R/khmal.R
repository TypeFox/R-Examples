"rqProcess" <-
function (formula, data, taus, nullH = "location", ...) 
{
    z <- summary(f <- lm(formula, data = data, x = TRUE))
    xbar <- apply(f$x,2,mean)
    vars <- names(z$coef[-1, 1])
    p <- length(z$coef[, 1])
    n <- nrow(z$x)
    Jinv <- z$cov.unscaled
    pivot <- any(taus < 0) || any(taus > 1)
    if(!pivot){ #grid method
       if(abs(diff(range(diff(taus)))) > sqrt(.Machine$double.eps))
		stop("rqProcess must be evaluated on equally spaced taus")
       ntaus <- length(taus)
       coefs <- matrix(0, ntaus, p)
       Cov <- array(0, c(p, p, ntaus))
       cat("taus: ")
       for (i in 1:ntaus) {
           cat(taus[i], " ")
           z <- summary(rq(formula, data = data, tau = taus[i], method = "fn"), 
		covariance = TRUE, ...)
           coefs[i, ] <- z$coef[, 1]
           Cov[, , i] <- z$cov/(taus[i]*(1 - taus[i]))
           }
    qtaus <- coefs %*% xbar
    Vhat <- t(coefs)[-1,,drop=FALSE]
    vhat <- t(coefs)[-1,,drop=FALSE]
    J <- solve(Jinv)
    p <- nrow(J)
    if (nullH == "location-scale") {
	f <- lm(coefs[,-1] ~ coefs[,1]) 
        b <- matrix(f$coef,2,p-1)[2,]
        R <- matrix(f$resid,ntaus,p-1)
        for (j in 1:length(taus)) {
            V <- Cov[, , j] 
            v <- V[-1, -1] + V[1, 1] * outer(b,b) - 
                outer(V[-1, 1], b) - t(outer(V[-1, 1], b))
            v <- solve(v)
            v <- chol(0.5 * (v + t(v)))
            Vhat[,j] <- v %*% R[j,]
            for (i in 2:p) {
                v <- V[i, i] + V[1, 1] * b[i-1]^2 - 2 * V[i, 1] * b[i-1]
                vhat[i-1,j] <- R[j, i-1]/sqrt(v)
               }
           }
        }
    else if (nullH == "location") {
        b <- apply(coefs, 2, mean)
        R <- t(coefs) - b
        for (j in 1:length(taus)) {
            V <- Cov[, , j] 
            A <- solve(V[-1, -1,drop=FALSE])
            B <- chol(0.5 * (A + t(A)))
            Vhat[,j] <- B %*% R[-1, j,drop=FALSE]
            vhat[,j] <- R[-1, j,drop=FALSE] / (sqrt(diag(V))[-1])
            }
         }
       }
   else{
	z <- rq(formula,data = data, tau = -1)
	taus <-  z$sol[1,]
	ntaus <- length(taus)
	qtaus <- z$sol[2,]
	qden <- qdensity(taus, qtaus)
	A <- solve(Jinv[-1,-1,drop=FALSE])
	B <- z$sol[-(1:3),,drop=FALSE]
	if(nullH == "location")
	   R <- B[-1,,drop=FALSE] - c(((B[-1,-1] + B[-1,-ntaus])/2) %*% diff(taus))
	else if(nullH == "location-scale")
	   R <- t(lm(t(B[-1,,drop=FALSE]) ~ B[1,])$resid)
	Vhat <- (chol(A) %*% t(t(R) * qden$s) )[,!qden$trim,drop=FALSE]
	vhat <- (t(t(R) * qden$s)/(sqrt(diag(Jinv))[-1]))[,!qden$trim,drop=FALSE]
	taus <- taus[!qden$trim]
	qtaus <- qtaus[!qden$trim]
	}
    dimnames(Vhat) <- list(vars, NULL)
    dimnames(vhat) <- list(vars, NULL)
    x <- list(taus = taus, qtaus = qtaus, Vhat = Vhat, vhat = vhat)
    class(x) <- "rqProcess"
    x
}
"qdensity" <- function(u,q,alpha = .05) {
#Computes Siddiqui estimate of quantile density function, fhat(F^{-1}(tau)) 
#Based on quantile regression process, using trimming proportion alpha 
#linear interpolation calls approx()  april, 2006
#h <- 0.6 * bandwidth.rq(u, length(u),hs=FALSE) #local bandwidth 
h <-  bandwidth.rq(u, length(u)) #local bandwidth 
trim <- ((u - h) < alpha) | ((u + h) > 1 - alpha)
qlo <- approx(u, q, u - h)$y
qup <- approx(u, q, u + h)$y
s <- (2 * h)/(qup - qlo)
list(s=s, trim = trim)
}
"KhmaladzeTest"  <-
function (formula, data = NULL, taus = -1, nullH = "location", 
	trim = c(0.05, 0.95), ...) 
{
    f <- rqProcess(formula, data = data, taus=taus, nullH = nullH, ...)
    Vtilde <- khmaladzize(f$taus, f$qtaus, f$Vhat, nullH)
    vtilde <- khmaladzize(f$taus, f$qtaus, f$vhat, nullH)
    trim <- (f$taus >= trim[1] & f$taus <= trim[2])
    Tvtilde <- (vtilde - vtilde[, 2])/sqrt(max(f$taus) - min(f$taus))
    TVtilde <- apply(abs(Vtilde - Vtilde[, 2])/
	sqrt(max(f$taus) - min(f$taus)), 2, "sum")[trim]
    Tn <- max(TVtilde)
    THn <- apply(abs(Tvtilde[, trim,drop = FALSE]), 1, max)
    x <- list(nullH = nullH, Tn = Tn, THn = THn)
    class(x) <- "KhmaladzeTest"
    x
}
"khmaladzize" <- 
function (taus, qtaus, Z, nullH) {
    dtaus <- diff(taus)
    dtaus <- c(dtaus[1], dtaus)
    score <- akj(qtaus, qtaus, dtaus)
    L <- length(taus)
    gdot2 <- -score$psi
    gdot <- cbind(rep(1, L), gdot2)
    if (nullH == "location-scale") {
        gdot3 <- gdot2 * qtaus
        gdot <- cbind(gdot, gdot3)
    }
    kmin <- 0
    p <- nrow(Z)
    for (i in 1:p) {
        v <- Z[i, ]
        dv <- c(0, diff(v))
        x1 <- gdot * sqrt(dtaus)
        x1 <- x1[L:1, ]
        y1 <- rev(dv/sqrt(dtaus))
        bhat <- lm.fit.recursive(x1, y1, int = FALSE)
        bhat <- bhat[, L:1]
        dvhat <- diag(gdot %*% bhat) * dtaus
        vhat <- cumsum(dvhat)
        v <- v[kmin:L - kmin]
        vhat <- vhat[kmin:L - kmin]
        Z[i, ] <- v - vhat
    }
    return(Z)
}
