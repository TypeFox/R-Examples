"lambdaPlot" <-
function (mod, lambda = seq(-1, 1, by = 0.1), stat = "F", global = TRUE, 
    cex = par("cex"), ...) 
{
    if (stat == "F") {
        org.fit <- mod
        y <- org.fit$model[, 1]
        resp <- names(org.fit$model)[1]
        print(resp)
        dat <- org.fit$model
        form <- as.formula(org.fit$call$formula)
        tt <- lm(form, data = dat)
        sav <- anova(tt)
        n <- nrow(sav)
        numdf <- sum(sav[-n, "Df"])
        dendf <- sav[n, "Df"]
        k <- length(lambda)
        if (global) {
            f.lambda <- matrix(NA, nrow = 1, ncol = k)
            dimnames(f.lambda) <- list("Model", paste("l", round(lambda, 
                2), sep = ""))
            for (j in seq(lambda)) {
                l <- lambda[j]
                if (l == 0) 
                  dat[, resp] <- log(y)
                else dat[, resp] <- (y^l - 1)/l
                tt <- lm(form, data = dat)
                f.lambda[1, j] <- (sum(anova(tt)[-n, "Sum Sq"])/numdf)/(anova(tt)[n, 
                  "Sum Sq"]/dendf)
            }
        }
        else {
            f.lambda <- matrix(NA, nrow = n, k)
            dimnames(f.lambda) <- list(dimnames(sav)[[1]], paste("l", 
                round(lambda, 2), sep = ""))
            for (j in seq(lambda)) {
                l <- lambda[j]
                if (l == 0) 
                  dat[, resp] <- log(y)
                else dat[, resp] <- (y^l - 1)/l
                tt <- lm(form, data = dat)
                f.lambda[, j] <- anova(tt)[, "F value"]
            }
            f.lambda <- f.lambda[-n, ]
        }
        Labels <- data.frame(term = dimnames(f.lambda)[[1]], 
            label = LETTERS[seq(nrow(f.lambda))])
        plot(lambda, lambda, xlim = range(lambda), ylim = range(f.lambda), 
            type = "n", xlab = "lambda", ylab = "F", ...)
        for (i in 1:nrow(f.lambda)) lines(lambda, f.lambda[i, 
            ])
        xlab <- lambda[k]
        ylab <- f.lambda[, k]
        lab <- paste(" ", as.character(Labels[, "label"]), sep = "")
        for (i in 1:nrow(f.lambda)) text(lambda[k], f.lambda[, 
            k], labels = lab, adj = 0)
        lab <- paste(as.character(Labels[, "label"]), " ", sep = "")
        for (i in 1:nrow(f.lambda)) text(lambda[1], f.lambda[, 
            1], labels = lab, adj = 1)
        print(Labels)
        invisible(list(lambda = lambda, f.lambda = f.lambda))
    }
    else if (stat == "t") {
        y <- mod$model[, 1]
        org.fit <- lm(y ~ ., qr = TRUE, data = mod$model[, -1])
        QR <- org.fit$qr
        n <- length(y)
        p <- length(coef(org.fit))
        idx <- 1:p
        rdf <- n - p
        coef.lambda <- matrix(NA, nrow = p, ncol = length(lambda))
        dimnames(coef.lambda) <- list(names(coef(org.fit)), paste("l", 
            round(lambda, 2), sep = ""))
        t.lambda <- se.lambda <- coef.lambda
        for (j in seq(lambda)) {
            l <- lambda[j]
            if (l == 0) 
                y.lambda <- log(y)
            else y.lambda <- (y^l - 1)/l
            resvar <- sum(qr.resid(QR, y.lambda)^2)/rdf
            coef.lambda[, j] <- qr.coef(QR, y.lambda)
            R <- chol2inv(QR$qr[idx, idx, drop = FALSE])
            se.lambda[, j] <- sqrt(diag(R) * resvar)
        }
        t.lambda <- coef.lambda/se.lambda
        Labels <- data.frame(term = names(coef(org.fit)), label = c(" ", 
            LETTERS[seq(nrow(t.lambda) - 1)]))
        plot(lambda, lambda, xlim = range(lambda), ylim = range(t.lambda[-1, 
            ]), type = "n", xlab = "lambda", ylab = "t", ...)
        for (i in 2:nrow(t.lambda)) lines(lambda, t.lambda[i, 
            ])
        xlab <- lambda[length(lambda)]
        ylab <- t.lambda[, ncol(t.lambda)]
        lab <- paste(" ", as.character(Labels[, "label"]), sep = "")
        for (i in 2:nrow(t.lambda)) text(lambda[length(lambda)], 
            t.lambda[, ncol(t.lambda)], labels = lab, adj = 0)
        lab <- paste(as.character(Labels[, "label"]), " ", sep = "")
        for (i in 2:nrow(t.lambda)) text(lambda[1], t.lambda[, 
            1], labels = lab, adj = 1)
        print(Labels)
        invisible(list(lambda = lambda, coef = coef.lambda, se = se.lambda))
    }
    else {
        warning("argument stat should be either \"F\" or \"t\"")
        invisible(NULL)
    }
}
