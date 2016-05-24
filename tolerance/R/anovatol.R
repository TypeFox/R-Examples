anovatol.int <- function (lm.out, data, alpha = 0.05, P = 0.99, side = 1, method = c("HE", "HE2",
	    "WBE", "ELL", "KM", "EXACT", "OCT"), m = 50) 
{
    method <- match.arg(method)
    out <- anova(lm.out)
    dims <- dim(out)
    s <- sqrt(out[dims[1], 3])
    df <- out[, 1]
    xlev <- lm.out$xlevels
    resp <- names(attr(lm.out$terms, "dataClasses"))[1]
    resp.ind <- which(colnames(data) == resp)
    pred.ind <- c(1:ncol(data))[c(colnames(data) %in% names(xlev))]
    factors <- names(xlev)
    out.list <- list()
    bal <- NULL
    for (i in 1:length(xlev)) {
        ttt <- by(data[, resp.ind], data[, pred.ind[i]], mean)
        temp.means <- as.numeric(ttt)
        temp.eff <- as.numeric(by(data[, resp.ind], data[, pred.ind[i]], 
            length))
        K <- NULL
        for (j in 1:length(temp.eff)) {
            K <- c(K, K.factor(n = temp.eff[j], f=tail(df,1), alpha = alpha, 
                P = P, side = side, method = method, m = m))
        }
        temp.low <- temp.means - K * s
        temp.up <- temp.means + K * s
        temp.mat <- data.frame(cbind(temp.means, temp.eff, K, 
            temp.low, temp.up))
        rownames(temp.mat) <- names(ttt)
        if (side == 1) {
            colnames(temp.mat) <- c("mean", "n", "k", "1-sided.lower", 
                "1-sided.upper")
        }
        else colnames(temp.mat) <- c("mean", "n", "k", "2-sided.lower", 
            "2-sided.upper")
        out.list[[i]] <- temp.mat
        bal = c(bal, sum(abs(temp.eff - mean(temp.eff)) > 3))
    }
    bal <- sum(bal)
    if (bal > 0) 
        warning("This procedure should only be used for balanced (or nearly-balanced) designs.", 
            call. = FALSE)
    names(out.list) <- names(xlev)
    if (side == 1) {
        message("These are ", (1 - alpha) * 100, "%/", P * 100, "% ", 
            side, "-sided tolerance limits.")
    }
    else message("These are ", (1 - alpha) * 100, "%/", P * 100, 
        "% ", side, "-sided tolerance intervals.")
    attr(out.list, "comment") <- c(resp, alpha, P)
    out.list
}
