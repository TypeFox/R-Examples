anova.rma <-
function (object, object2, btt, L, digits, ...) 
{
    if (!is.element("rma", class(object))) 
        stop("Argument 'object' must be an object of class \"rma\".")
    if (any(is.element(c("rma.mh", "rma.peto"), class(object)))) 
        stop("Function not applicable for objects of class \"rma.mh\" or \"rma.peto\".")
    if (is.element("rma.glmm", class(object))) 
        stop("Method not yet implemented for objects of class \"rma.glmm\". Sorry!")
    if (missing(digits)) 
        digits <- object$digits
    if (missing(object2)) {
        x <- object
        k <- x$k
        p <- x$p
        b <- x$b
        vb <- x$vb
        if (missing(L)) {
            if (missing(btt)) 
                btt <- NULL
            if (is.null(btt)) {
                if (p > 1) {
                  if (x$int.incl) {
                    btt <- seq.int(from = 2, to = p)
                  }
                  else {
                    btt <- seq_len(p)
                  }
                }
                else {
                  btt <- 1
                }
            }
            else {
                btt <- btt[(btt >= 1) & (btt <= p)]
                btt <- unique(round(btt))
                if (length(intersect(btt, seq_len(p))) == 0L) {
                  stop("Non-existent coefficients specified via 'btt'.")
                }
            }
            m <- length(btt)
            QM <- c(t(b)[btt] %*% chol2inv(chol(vb[btt, btt])) %*% 
                b[btt])
            if (x$knha) {
                QM <- QM/m
                QMp <- pf(QM, df1 = m, df2 = x$dfs, lower.tail = FALSE)
            }
            else {
                QMp <- pchisq(QM, df = m, lower.tail = FALSE)
            }
            res <- list(QM = QM, QMp = QMp, btt = btt, k = k, 
                p = p, m = m, knha = x$knha, dfs = x$dfs, digits = digits, 
                test = "Wald.b")
        }
        else {
            if (is.vector(L)) 
                L <- rbind(L)
            if (is.data.frame(L)) 
                L <- as.matrix(L)
            if (is.character(L)) 
                stop("Argument 'L' must be a numeric vector/matrix.")
            if (x$int.incl && ncol(L) == (p - 1)) 
                L <- cbind(1, L)
            if (ncol(L) != p) 
                stop("Length or number of columns of 'L' does not match number of model coefficients.")
            m <- nrow(L)
            Lb <- L %*% b
            vLb <- L %*% vb %*% t(L)
            se <- sqrt(diag(vLb))
            zval <- c(Lb/se)
            if (x$knha) {
                pval <- 2 * pt(abs(zval), df = x$dfs, lower.tail = FALSE)
            }
            else {
                pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
            }
            QM <- NULL
            QMp <- NULL
            if (rankMatrix(L) == m) {
                QM <- try(t(Lb) %*% chol2inv(chol(vLb)) %*% Lb, 
                  silent = TRUE)
                if (!inherits(QM, "try-error")) {
                  if (x$knha) {
                    QM <- QM/m
                    QMp <- pf(QM, df1 = m, df2 = x$dfs, lower.tail = FALSE)
                  }
                  else {
                    QMp <- pchisq(QM, df = m, lower.tail = FALSE)
                  }
                }
            }
            hyp <- rep("", m)
            for (j in 1:m) {
                Lj <- L[j, ]
                sel <- Lj != 0
                hyp[j] <- paste(paste(Lj[sel], rownames(b)[sel], 
                  sep = "*"), collapse = " + ")
                hyp[j] <- gsub("1*", "", hyp[j], fixed = TRUE)
                hyp[j] <- gsub("+ -", "- ", hyp[j], fixed = TRUE)
            }
            hyp <- paste0(hyp, " = 0")
            hyp <- data.frame(hyp)
            colnames(hyp) <- ""
            rownames(hyp) <- paste0(seq_len(m), ":")
            res <- list(QM = QM, QMp = QMp, hyp = hyp, Lb = Lb, 
                se = se, zval = zval, pval = pval, k = k, p = p, 
                m = m, knha = x$knha, dfs = x$dfs, digits = digits, 
                test = "Wald.L")
        }
    }
    else {
        if (!is.element("rma", class(object2))) 
            stop("Argument 'object2' must be an object of class \"rma\".")
        if (any(is.element(c("rma.mh", "rma.peto"), class(object2)))) 
            stop("Function not applicable for objects of class \"rma.mh\" or \"rma.peto\".")
        if (is.element("rma.glmm", class(object2))) 
            stop("Method not yet implemented for objects of class \"rma.glmm\". Sorry!")
        if (!identical(class(object), class(object2))) 
            stop("Class of 'object1' must be the same as class of 'object2'.")
        m.f <- object
        m.r <- object2
        if (is.element("rma.uni", class(object))) {
            if (!(identical(as.vector(m.f$yi), as.vector(m.r$yi)) && 
                identical(as.vector(m.f$vi), as.vector(m.r$vi)))) 
                stop("Observed outcomes and/or sampling variances not equal in the full and reduced model.")
        }
        else {
            if (!(identical(as.vector(m.f$yi), as.vector(m.r$yi)) && 
                identical(m.f$V, m.r$V))) 
                stop("Observed outcomes and/or sampling variances/covariances not equal in the full and reduced model.")
        }
        p.f <- m.f$parms
        p.r <- m.r$parms
        if (p.f == p.r) 
            stop("Models have the same number of parameters. LRT not meaningful.")
        if (p.f < p.r) {
            m.f <- object2
            m.r <- object
            p.s <- p.f
            p.f <- p.r
            p.r <- p.s
        }
        if (m.f$method == "FE" && m.r$method != "FE") 
            stop("Full model uses a fixed- and reduced model uses random/mixed-effects model.")
        p.diff <- p.f - p.r
        if (m.f$method == "REML") {
            LRT <- abs(m.r$fit.stats$REML[2] - m.f$fit.stats$REML[2])
            fit.stats.f <- m.f$fit.stats$REML
            fit.stats.r <- m.r$fit.stats$REML
            if (!identical(m.f$X, m.r$X)) 
                warning("Models with different fixed effects. REML comparisons are not meaningful.")
        }
        else {
            LRT <- abs(m.r$fit.stats$ML[2] - m.f$fit.stats$ML[2])
            fit.stats.f <- m.f$fit.stats$ML
            fit.stats.r <- m.r$fit.stats$ML
        }
        pval <- pchisq(LRT, df = p.diff, lower.tail = FALSE)
        if (is.element("rma.uni", class(object))) {
            if (m.f$method == "FE" || identical(m.r$tau2, 0)) {
                R2 <- NA
            }
            else {
                R2 <- round(100 * max(0, (m.r$tau2 - m.f$tau2)/m.r$tau2), 
                  2)
            }
        }
        else {
            R2 <- NA
        }
        if (is.element("rma.uni", class(object))) {
            tau2.f <- m.f$tau2
            tau2.r <- m.r$tau2
        }
        else {
            tau2.f <- NA
            tau2.r <- NA
        }
        res <- list(fit.stats.f = fit.stats.f, fit.stats.r = fit.stats.r, 
            p.f = p.f, p.r = p.r, LRT = LRT, pval = pval, QE.f = m.f$QE, 
            QE.r = m.r$QE, tau2.f = tau2.f, tau2.r = tau2.r, 
            R2 = R2, method = m.f$method, class.f = class(m.f), 
            digits = digits, test = "LRT")
    }
    class(res) <- c("anova.rma")
    return(res)
}
