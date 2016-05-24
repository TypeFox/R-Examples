regtol.int <- function (reg, new.x = NULL, side = 1, alpha = 0.05, P = 0.99) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (class(reg) != "lm") {
        stop(paste("Input must be of class 'lm'.", "\n"))
    }
    if (is.vector(new.x)) 
        new.x <- matrix(new.x, ncol = 1)
    n <- length(reg$res)
    pars <- length(reg$coef)
    if(is.null(reg$weights)){
      x <- data.matrix(reg$model[, -1], rownames.force = FALSE)
    } else{
      x <- data.matrix(reg$model[, -c(1,ncol(reg$model))], rownames.force = FALSE)
    }
    new.length <- 0
    if (is.null(new.x) == FALSE) {
        new.length <- nrow(new.x)
        x <- rbind(x, new.x)
    }
    x <- data.frame(x)
    y <- c(reg$model[, 1], rep(NA, new.length))
    names(x) <- names(reg$coef[-1])
    est <- predict(reg, newdata = x, se.fit = TRUE)
    y.hat <- est$fit
    se.y <- est$se.fit
    a.out <- anova(reg)
    MSE <- a.out$"Mean Sq"[length(a.out$"Mean Sq")]
    df <- a.out$Df[length(a.out$Df)]
    n.star <- MSE/se.y^2
    if (side == 1) {
        z.p <- qnorm(P)
        delta <- sqrt(n.star) * z.p
        t.delta <- suppressWarnings(qt(1 - alpha, df = n - pars, 
            ncp = delta))
        K <- t.delta/sqrt(n.star)
        upper <- y.hat + sqrt(MSE) * K
        lower <- y.hat - sqrt(MSE) * K
        temp <- data.frame(cbind(alpha, P, y, y.hat, lower, upper))
        colnames(temp) <- c("alpha", "P", "y", "y.hat", "1-sided.lower", 
            "1-sided.upper")
    }
    else {
	K <- sqrt(df * qchisq(P, 1, 1/n.star)/qchisq(alpha, df))
        upper <- y.hat + sqrt(MSE) * K
        lower <- y.hat - sqrt(MSE) * K
        temp <- data.frame(cbind(alpha, P, y, y.hat, lower, upper))
        colnames(temp) <- c("alpha", "P", "y", "y.hat", "2-sided.lower", 
            "2-sided.upper")
    }
    index <- which(names(temp) == "y.hat")
    temp <- data.matrix(temp[order(temp[, index]), ], rownames.force = FALSE)
    temp <- data.frame(temp, check.names = FALSE)
    temp
}

