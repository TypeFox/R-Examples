llwr <-
function (theta, h, Bd, Pr, 
    particle.density, air, 
    critical.PR, h.FC, h.WP, 
    water.model = c("Silva", "Ross"), 
    Pr.model = c("Busscher", "noBd"), 
    pars.water = NULL, pars.Pr = NULL, 
    graph = TRUE, graph2 = TRUE, 
    xlab = expression(Bulk ~ Density ~ (Mg ~ m^{-3})), 
    ylab = expression(theta ~ (m^{3} ~ m^{-3})), 
    main = "Least Limiting Water Range", ...) 
{
    n <- length(theta)
    if (length(h) != n || length(Pr) != n) 
        stop("incompatible dimensions!")
    dat <- cbind(theta, h, Bd, Pr)
    if (!is.numeric(dat)) 
        stop("non-numeric data!")
    limits <- c(particle.density, air, critical.PR, h.FC, h.WP)
    if (length(limits) > 5) 
        stop("each limiting value must be a single value!")
    # water model ------------------------------------------
    water.model <- match.arg(water.model)
    if (is.null(pars.water)) {
        if (water.model == "Silva" & length(unique(Bd)) > 1L) {
            fit1. <- lm(log(theta) ~ Bd + log(h))
            a. <- coef(fit1.)
            fit1 <- nls(theta ~ exp(a + b * Bd) * h^c, 
                start = list(a = a.[1], b = a.[2], c = a.[3]))
            a <- coef(fit1)
        }
        else {
            fit1. <- lm(log(theta) ~ log(h))
            a. <- coef(fit1.)
            fit1 <- nls(theta ~ a * h^b, 
                start = list(a = exp(a.[1]), b = a.[2]))
            a <- c(log(coef(fit1)[1]), 0, coef(fit1)[2])
        }
        rsq1 <- Rsq(fit1)
    }
    else if (length(pars.water) == 3) {
        a <- pars.water
    }
    # PR model --------------------------------------------
    Pr.model <- match.arg(Pr.model)
    if (is.null(pars.Pr)) {
        if (Pr.model == "Busscher" & length(unique(Bd)) > 1L) {
            fit2 <- fitbusscher(Pr, theta, Bd)
            b <- coef(fit2)
        }
        else {
            fit2. <- lm(log(Pr) ~ log(theta))
            b. <- coef(fit2.)
            fit2 <- nls(Pr ~ b0 * theta^b1, 
                start = list(b0 = exp(b.[1]), b1 = b.[2]))
            b <- c(coef(fit2)[1], coef(fit2)[2], 0)
        }
        rsq2 <- Rsq(fit2)
    }
    else if (length(pars.Pr) == 3) {
        b <- pars.Pr
    }
    # limiting functions -----------------------------------
    Dp <- particle.density
    thetaA <- 1 - Bd/Dp - air
    PRc <- critical.PR
    thetaPR <- (PRc/(b[1] * Bd^b[3]))^(1/b[2])
    thetaFC <- exp(a[1] + a[2] * Bd) * h.FC^a[3]
    thetaWP <- exp(a[1] + a[2] * Bd) * h.WP^a[3]
    theta. <- cbind(thetaA, thetaPR, thetaFC, thetaWP)
    if (length(unique(Bd)) > 1L) {
        x. <- seq(range(Bd)[1], range(Bd)[2], length.out = 100)
        thetaA. <- 1 - x./Dp - air
        thetaPR. <- (PRc/(b[1] * x.^b[3]))^(1/b[2])
        mi <- which.min((thetaA. - thetaPR.)^2)
        x <- seq(range(Bd)[1], x.[mi], length.out = 100)
    }
    else {
        x <- Bd
    }
    # defining LLWR ----------------------------------------
    yUp. <- cbind(1 - x/Dp - air, 
        exp(a[1] + a[2] * x) * h.FC^a[3])
    yUp <- apply(yUp., 1, min)
    yLow. <- cbind((PRc/(b[1] * x^b[3]))^(1/b[2]), 
        exp(a[1] + a[2] * x) * h.WP^a[3])
    yLow <- apply(yLow., 1, max)
    iho <- as.vector(yUp - yLow)
    # graph ------------------------------------------------
    if (graph) {
        if (length(unique(Bd)) > 1L) {
            plot(range(Bd), range(theta.) * c(1, 1.1), pch = "", 
                xlab = xlab, ylab = ylab, main = main, ...)
            polygon(c(x[1], x, x[100]), c(yLow[1], yUp, yLow[100]), 
                col = gray(0.8, alpha = 1), border = FALSE)
            polygon(c(x[1], x, x[100]), c(yUp[1], yLow, yUp[100]), 
                col = gray(0.8, alpha = 1), border = FALSE)
            curve(1 - x/Dp - air, add = TRUE)
            curve(exp(a[1] + a[2] * x) * h.FC^a[3], add = TRUE)
            curve((PRc/(b[1] * x^b[3]))^(1/b[2]), add = TRUE, 
                col = "blue", lty = 2)
            curve(exp(a[1] + a[2] * x) * h.WP^a[3], add = TRUE, 
                col = "blue", lty = 2)
            points(Bd, thetaA)
            points(Bd, thetaFC, pch = 16)
            points(Bd, thetaPR, col = "blue")
            points(Bd, thetaWP, col = "blue", pch = 16)
            tex <- c(expression(theta[A]), expression(theta[FC]), 
                expression(theta[PR]), expression(theta[WP]))
            legend("topright", tex, lty = c(1, 1, 2, 2), col = c(1, 
                1, 4, 4), pch = c(1, 16, 1, 16), cex = 0.8, bg = "white")
            if (graph2) {
                dev.new(width = 3, height = 3)
                plot(x, yUp - yLow, type = "l", xlab = xlab, 
                  ylab = "LLWR", ...)
            }
        }
        else {
            plot(rep(Bd, 4), theta., ylab = ylab, xlab = xlab, 
                xaxt = "n", ...)
            axis(1, Bd, Bd)
            text(rep(Bd, 4) * 1.05, theta., c("A", "PR", "FC", 
                "WP"), cex = 0.8)
            polygon(x = c(-99, -99, 99, 99), y = c(yUp, yLow, 
                yLow, yUp), border = NA, col = adjustcolor("blue", 
                alpha.f = 0.1))
            arrows(Bd, yLow, Bd, yUp, length = 0.08, angle = 90, 
                code = 3, col = 4, ...)
        }
    }
    # output -------------------------------------------------
    if (length(unique(Bd)) > 1L) 
        area <- trapez(x, yUp) - trapez(x, yLow)
    out <- list(limiting.theta = drop(theta.), 
        pars.water = if (is.null(pars.water)) fit1 else a, 
        r.squared.water = if (is.null(pars.water)) rsq1, 
        pars.Pr = if (is.null(pars.Pr)) fit2 else b, 
        r.squared.Pr = if (is.null(pars.Pr)) rsq2, 
        area = if (length(unique(Bd)) > 1L) area, 
        LLWR = if (length(unique(Bd)) == 1L) iho)
    class(out) <- "llwr"
    return(out)
}

# -----------------------------------------------
# print method
print.llwr <-
function (x, ...) 
{
    cat("\n          Least Limiting Water Range \n")
    if (!is.null(x$area)) {
        cat("\n---------- \nLimiting theta (6 first rows):\n")
        print(head(x$limiting.theta))
    }
    else {
        cat("\n---------- \nLimiting theta:\n")
        print(drop(x$limiting.theta))
    }
    cat("\n---------- \nEstimates of the soil water content model: \n")
    if (inherits(x$pars.water, "nls")) {
        print(summary(x$pars.water))
        rsq1 <- x$r.squared.water$pseudo.R.squared
        cat("pseudo R-squared:", rsq1, "\n")
        adj.rsq1 <- x$r.squared.water$adj.R.squared
        cat("adjusted R-squared:", adj.rsq1, "\n")
    }
    else if (inherits(x$pars.water, "numeric")) {
        print(x$pars.water)
    }
    cat("\n---------- \nEstimates of the soil penetration resistance model: \n")
    if (inherits(x$pars.Pr, "nls")) {
        print(summary(x$pars.Pr))
        rsq2 <- x$r.squared.Pr$pseudo.R.squared
        cat("pseudo R-squared:", rsq2, "\n")
        adj.rsq2 <- x$r.squared.Pr$adj.R.squared
        cat("adjusted R-squared:", adj.rsq2, "\n")
    }
    else if (inherits(x$pars.Pr, "numeric")) {
        print(x$pars.Pr)
    }
    if (!is.null(x$area)) {
        cat("\n---------- \nShaded area:", x$area, "\n")
    }
    else {
        cat("\n---------- \nLLWR:", x$LLWR, "\n")
    }
    invisible(x)
}

# ---------------------------------------------------
# numerical integration: trapezoidal rule
trapez <- 
function(x, y)
{
   if (length(x) != length(y))
      stop("incompatible dimensions!")
   xy <- sortedXyData(x, y)
   x <- xy[["x"]]
   y <- xy[["y"]]
   n <- nrow(xy)
   out <- 0.5 * sum((x[2:n] - x[2:n - 1]) * 
      (y[2:n] + y[2:n - 1]))
   return(out)
}
