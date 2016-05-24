Sindex <-
function (theta_R, theta_S, alpha, n, m = 1 - 1/n, vcov = NULL,
      nsim = 999, conf.level = 0.95, graph = TRUE, ...) 
{
    x <- NULL
    h_i <- 1/alpha * (1/m)^(1/n)
    theta_i <- theta_R + (theta_S - theta_R) * (1 + 1/m)^(-m)
    b1 <- -n * (theta_S - theta_R) * (1 + 1/m)^(-(1 + m))
    b0 <- theta_i - b1 * h_i
    if (graph) {
        curve(soilwater(x, theta_R, theta_S, alpha, n, m),
            xlab = "Matric potential", 
            ylab = "Soil water content",
            main = "Soil Water Retention Curve", ...)
        lines(x = c(h_i, h_i), y = c(-1e+09, theta_i), lty = 3)
        lines(x = c(-1e+09, h_i), y = c(theta_i, theta_i), lty = 3)
    }
    S <- abs(b1)
    if (S >= 0.05) {
        clas <- "Very good"
    }
    else if (S < 0.05 & S >= 0.035) {
        clas <- "Good"
    }
    else if (S < 0.035 & S > 0.02) {
        clas <- "Poor"
    }
    else {
        clas <- "Very poor"
    }
    if (!is.null(vcov)) {
       mu <- c(theta_R, theta_S, alpha, n)
       if (length(mu) != nrow(vcov))
          stop("'vcov' misspecified!")
       sim <- mvrnorm(nsim, mu, vcov)
       m.sim <- 1 - 1/sim[, 4]
       sim <- cbind(sim, m.sim)
       S.sim <- c()
       for(i in 1:nsim) {
          S.sim[i] <- abs(-sim[i, 4] * (sim[i, 2] - sim[i, 1]) *
             (1 + 1/sim[i, 5])^(-(1 + sim[i, 5])))
       }
       sig <- (1 - conf.level)/2
       sci <- as.vector(quantile(S.sim, p = c(sig, 1 - sig)))
    }    

    if (!is.null(vcov)) {
        out <- list(h_i = as.vector(h_i), theta_i = as.vector(theta_i),
           S.index = as.vector(S), PhysicalQuality = clas,
               simCI = sci, conf.level = conf.level)
    } else {
        out <- list(h_i = as.vector(h_i), theta_i = as.vector(theta_i),
           S.index = as.vector(S), PhysicalQuality = clas, simCI = NULL)
    }
    class(out) <- "S.index"
    return(out)
}


# ----------------------------------------
# print method
print.S.index <- function(x, ...)
{
    cat("\n          The S Index \n")
    cat("\nh_i :", round(x$h_i, 4),
        "\ntheta_i :", round(x$theta_i, 4),
        "\n|S| :", round(x$S.index, 4),
        "\nSoil physical quality :", x$PhysicalQuality, "\n")
    if (!is.null(x$simCI)) {
        cat("Simulated CI", " (", 100*x$conf.level, "%) : ",
           "[", x$simCI[1], ", ", x$simCI[2], "]\n", sep = "")
    }
    invisible(x)
}