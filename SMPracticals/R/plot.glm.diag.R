"plot.glm.diag" <-
function (x, glmdiag = glm.diag(x), subset = NULL, 
    iden = FALSE, labels = NULL, ret = FALSE, ...) 
{
    if (is.null(glmdiag)) 
        glmdiag <- glm.diag(x)
    if (is.null(subset)) 
        subset <- c(1:length(glmdiag$h))
    else if (is.logical(subset)) 
        subset <- (1:length(subset))[subset]
    else if (is.numeric(subset) && all(subset < 0)) 
        subset <- (1:(length(subset) + length(glmdiag$h)))[subset]
    else if (is.character(subset)) {
        if (is.null(labels)) 
            labels <- subset
        subset <- seq(along = subset)
    }
    par(mfrow = c(2, 2))
    x1 <- predict(x)
    plot(x1, glmdiag$res, xlab = "Linear predictor", ylab = "Residuals")
    pars <- vector(4, mode = "list")
    pars[[1]] <- par("usr")
    y2 <- glmdiag$rd
    x2 <- qnorm(ppoints(length(y2)))[rank(y2)]
    plot(x2, y2, ylab = "Quantiles of standard normal", xlab = "Ordered deviance residuals")
    abline(0, 1, lty = 2)
    pars[[2]] <- par("usr")
    hh <- glmdiag$h/(1 - glmdiag$h)
    plot(hh, glmdiag$cook, xlab = "h/(1-h)", ylab = "Cook statistic")
    rx <- range(hh)
    ry <- range(glmdiag$cook)
    rank.fit <- x$rank
    nobs <- rank.fit + x$df.residual
    cooky <- 8/(nobs - 2 * rank.fit)
    hy <- (2 * rank.fit)/(nobs - 2 * rank.fit)
    if ((cooky >= ry[1]) && (cooky <= ry[2])) 
        abline(h = cooky, lty = 2)
    if ((hy >= rx[1]) && (hy <= rx[2])) 
        abline(v = hy, lty = 2)
    pars[[3]] <- par("usr")
    plot(subset, glmdiag$cook, xlab = "Case", ylab = "Cook statistic")
    if ((cooky >= ry[1]) && (cooky <= ry[2])) 
        abline(h = cooky, lty = 2)
    xx <- list(x1, x2, hh, subset)
    yy <- list(glmdiag$res, y2, glmdiag$cook, glmdiag$cook)
    pars[[4]] <- par("usr")
    if (is.null(labels)) 
        labels <- names(x1)
    while (iden) {
        cat("****************************************************\n")
        cat("Please Input a screen number (1,2,3 or 4)\n")
        cat("0 will terminate the function \n")
        num <- as.numeric(readline())
        if ((length(num) > 0) && ((num == 1) || (num == 2) || 
            (num == 3) || (num == 4))) {
            cat(paste("Interactive Identification for screen", 
                num, "\n"))
            cat("left button = Identify, center button = Exit\n")
            nm <- num + 1
            par(mfg = c(trunc(nm/2), 1 + nm%%2, 2, 2))
            par(usr = pars[[num]])
            identify(xx[[num]], yy[[num]], labels)
        }
        else iden <- FALSE
    }
    par(mfrow = c(1, 1))
    if (ret) 
        glmdiag
    else invisible()
}

