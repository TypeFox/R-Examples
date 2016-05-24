hemc <- 
function (x, theta_R, theta_S, alpha, n, b1, b2, graph = TRUE, 
    from = 1, to = 30, xlab = expression(Psi ~ (J~kg^{-1})), ylab = expression(d ~ 
        theta/d ~ Psi), ...) 
{
    if (!inherits(x, c("numeric", "integer"))) 
        stop("non-numeric argument!")
    lengths <- sapply(list(theta_R, theta_S, alpha, n, b1, b2), length)
    if (any(lengths != 2))
        stop("all the parameters must be a numeric vector of length two, in the following order: 'fast' and 'slow')!")

    # Required functions ('der' and 'peak')
    der <- function(x, theta_R, theta_S, alpha, n, b1, b2) {
        out <- abs((theta_S - theta_R) * (1 + (alpha * x)^n)^(1/n - 
            1) * (1/n - 1) * (alpha * x)^n * (n/(x * (1 + (alpha * 
            x)^n))) + 2 * b2 * x + b1)
        return(out)
    }
    peak <- function(ran, theta_R, theta_S, alpha, n, b1, b2) {
        x <- seq(ran[1], ran[2], length.out = 1000)
        y <- der(x, theta_R, theta_S, alpha, n, b1, b2)
        indMax <- which.max(y)
        return(x[indMax])
    }

    # Modal Suction
    p0Fast <- peak(range(x), theta_R[1], theta_S[1], alpha[1], 
        n[1], b1[1], b2[1])
    y0Fast <- der(p0Fast, theta_R[1], theta_S[1], alpha[1], 
       n[1], b1[1], b2[1])
    p0Slow <- peak(range(x), theta_R[2], theta_S[2], alpha[2], 
        n[2], b1[2], b2[2])
    y0Slow <- der(p0Slow, theta_R[2], theta_S[2], alpha[2], 
       n[2], b1[2], b2[2])

    # Graphic
    if (graph) {
        curve(der(x, theta_R[1], theta_S[1], alpha[1], n[1], b1[1], b2[1]), 
           from = from, to = to, lwd = 2, xlab = xlab, ylab = ylab,
           col = "blue", ...)
        curve(abs(2 * b2[1] * x + b1[1]), add = TRUE, col = "gray")
        lines(x = c(p0Fast, p0Fast), y = c(0, y0Fast), col = "gray", lty = 3)
        curve(der(x, theta_R[2], theta_S[2], alpha[2], n[2], b1[2], b2[2]), 
           add = TRUE, lwd = 2, lty = 2, col = "red", ...)
        curve(abs(2 * b2[2] * x + b1[2]), add = TRUE, col = "gray")
        lines(x = c(p0Slow, p0Slow), y = c(0, y0Slow), col = "gray", lty = 3)
        legend("topright", c("Fast", "Slow"), lty = 1:2, lwd = 2, 
           col = c("blue", "red"), cex = 0.8, bg = "lightyellow")
    }

    # VDP (area)
    f.area <- function(ran, theta_R, theta_S, alpha, n, b1, b2) {
        x <- seq(ran[1], ran[2], length.out = 1000)
        y <- abs((theta_S - theta_R) * (1 + (alpha * x)^n)^(1/n - 
            1) * (1/n - 1) * (alpha * x)^n * (n/(x * (1 + (alpha * 
            x)^n))))
        return(simpson(x, y))
    }
    areaFast <- f.area(range(x), theta_R[1], theta_S[1], alpha[1], 
       n[1], b1[1], b2[1])
    areaSlow <- f.area(range(x), theta_R[2], theta_S[2], alpha[2], 
       n[2], b1[2], b2[2])

    # Structural index
    SIFast <- areaFast / p0Fast
    SISlow <- areaSlow / p0Slow

    # Stability ratio
    SR <- SIFast / SISlow

    # Output
    mout <- matrix(c(p0Fast, areaFast, SIFast, p0Slow, areaSlow, SISlow), 
       nrow = 3, ncol = 2, dimnames = list(c("Modal Suction", "VDP", "Str Index"),
       c("Fast", "Slow")))
    out <- list(mout, StabilityRatio = SR)
    class(out) <- "hemc"
    return(out)
}

# ------------------------------------------------------
# print method
# ------------------------------------------------------
# print method
print.hemc <- function(x, digits = 4, ...)
{
    cat("\n          High-Energy-Moisture-Characteristics  \n\n")
    print(round(x[[1]], digits))
    cat("\nStability Ratio: ", x$StabilityRatio, "\n")
    invisible(x)
}

# ------------------------------------------------------
# Simpson's integration rule
# ------------------------------------------------------
simpson <- 
function(x, y)
{
   if (length(x) != length(y))
      stop("incompatible dimensions!")
   xy <- sortedXyData(x, y)
   x <- xy[["x"]]
   y <- xy[["y"]]
   n <- nrow(xy)
   if (n %% 2 != 0)
      stop("the length of both vectors must be an even number!")
   delta <- (x[n] - x[1]) / n
   mult <- c(1, rep(c(4, 2), (n - 2)/2), 1)
   out <- delta / 3 * sum(mult * y)
   return(out)
}
