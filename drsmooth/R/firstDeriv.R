#' @name firstDeriv
#' @title First Derivative Function(s)
#' @param mod  The gam model.
#' @param n   Prediction increments.
#' @importFrom stats binomial terms
#' @keywords internal


firstDeriv <- function(mod, n) {

    if(isTRUE(all.equal(class(mod), "list")))
    
    # Future: sensitivity of first derivative should not be hard-coded to 0.05.
    mod <- mod$gam
    eps <- 1e-7
    alpha <- .05
    m.terms <- attr(terms(mod), "term.labels")
    newD <- sapply(stats::model.frame(mod)[, m.terms, drop = FALSE],
                       function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms

    X0 <- stats::predict(mod, data.frame(newD), type = "lpmatrix")
    newD <- newD + eps
    X1 <- stats::predict(mod, data.frame(newD), type = "lpmatrix")
    Xp <- (X1 - X0) / eps
    Xp.r <- NROW(Xp)
    Xp.c <- NCOL(Xp)
    ## dims of bs
    bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
    ## number of smooth terms
    t.labs <- attr(mod$terms, "term.labels")
    nt <- length(t.labs)
    ## list to hold the derivatives
    lD <- vector(mode = "list", length = nt)
    names(lD) <- t.labs
    for(i in seq_len(nt)) 
	{
        Xi <- Xp * 0
        want <- grep(t.labs[i], colnames(X1))
        Xi[, want] <- Xp[, want]
        df <- Xi %*% stats::coef(mod)
        df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
        lD[[i]] <- list(deriv = df, se.deriv = df.sd)
    	}
    class(lD) <- "firstDeriv"
    lD$gamModel <- mod
    lD$eps <- eps
    lD$eval <- newD - eps
    return(lD)
}


