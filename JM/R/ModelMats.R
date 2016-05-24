ModelMats <-
function (time, ii, obs.times, survTimes) {
    if (method %in% c("weibull-AFT-GH", "weibull-PH-GH", 
            "spline-PH-GH", "spline-PH-Laplace")) {
        id.GK <- if (!LongFormat) {
            rep(ii, each = object$control$GKk)
        } else {
            rep(which(idT == ii), each = object$control$GKk)
        }
        wk <- gaussKronrod(object$control$GKk)$wk
        sk <- gaussKronrod(object$control$GKk)$sk
        if (!LongFormat) {
            P <- time / 2
            st <- P * (sk + 1)
        } else {
            time0 <- obs.times[[ii]]
            time1 <- c(time0[-1], time)
            P <- (time1 - time0) / 2
            P1 <- (time1 + time0) / 2
            st <- outer(P, sk) + P1
            st <- c(t(st))
        }
        data.id2 <- data.id[id.GK, ]
        data.id2[[timeVar]] <- pmax(st - lag, 0)
        out <- list(st = st, wk = rep(wk, length(P)), P = P)
        if (parameterization %in% c("value", "both")) {
            mfX <- model.frame(delete.response(TermsX), data = data.id2)
            mfZ <- model.frame(TermsZ, data = data.id2)
            out$Xs <- model.matrix(formYx, mfX)
            out$Zs <- model.matrix(formYz, mfZ)
            out$Ws.intF.vl <- WintF.vl[id.GK, , drop = FALSE]
        }
        if (parameterization %in% c("slope", "both")) {
            mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
            mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
            out$Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
            out$Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
            out$Ws.intF.sl <- WintF.sl[id.GK, , drop = FALSE]
        }
    }
    if (method == "piecewise-PH-GH") {
        wk <- gaussKronrod(7)$wk
        sk <- gaussKronrod(7)$sk
        nk <- length(sk)
        qs <- c(0, sort(object$control$knots), 
            max(survTimes, object$control$knots) + 1)
        ind <- findInterval(time, qs, rightmost.closed = TRUE)
        Tiq <- outer(time, qs, pmin)
        Lo <- Tiq[, 1:Q]
        Up <- Tiq[, 2:(Q+1)]
        T <- Up - Lo
        P <- T / 2
        if (!all(P < sqrt(.Machine$double.eps)))
            P[P < sqrt(.Machine$double.eps)] <- as.numeric(NA)
        P1 <- (Up + Lo) / 2
        st <- rep(P, each = nk) * rep(sk, Q) + rep(P1, each = nk)
        data.id2 <- data.id[rep(ii, each = nk*Q), ]
        data.id2[[timeVar]] <- pmax(st - lag, 0)
        data.id2 <- data.id2[!is.na(st), ]
        id.GK <- rep(ii, sum(!is.na(st)))
        out <- list(st = st, wk = wk, P = P, ind = ind)
        if (parameterization %in% c("value", "both")) {
            mfX <- model.frame(TermsX, data = data.id2)
            mfZ <- model.frame(TermsZ, data = data.id2)
            out$Xs <- model.matrix(formYx, mfX)
            out$Zs <- model.matrix(object$formYz, mfZ)
            out$Ws.intF.vl <- WintF.vl[id.GK, , drop = FALSE]
        }
        if (parameterization %in% c("slope", "both")) {
            mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
            mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
            out$Xs.deriv <- model.matrix(derivForm$fixed, mfX.deriv)
            out$Zs.deriv <- model.matrix(derivForm$random, mfZ.deriv)
            out$Ws.intF.sl <- WintF.sl[id.GK, , drop = FALSE]
        }
    }
    out
}
