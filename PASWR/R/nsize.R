nsize <-
function(b, sigma = NULL, p = 0.5, conf.level = 0.95, type = "mu")
{
    choices <- c("mu", "pi")
    alt <- pmatch(type, choices)
    type <- choices[alt]
    if(length(type) > 1 || is.na(type))
        stop("type must be one \"mu\", \"pi\"")
    if(type == "pi" && b > 1)
        stop("b must be less than 1")
    if(!missing(b))
        if(length(b) != 1 || is.na(b))
            stop("b must be a single number")
    if(type == "mu") {
        z <- qnorm(1 - (1 - conf.level)/2)
        n <- ((z * sigma)/b)^2
        n <- ceiling(n)
        cat("\n")
        cat("The required sample size (n) to estimate the population",
            "\n")
        cat("mean with a", conf.level,
            "confidence interval so that the margin", "\n")
        cat("of error is no more than", b, "is", n,".", "\n")
        cat("\n\n")
    }
    else if(type == "pi") {
        z <- qnorm(1 - (1 - conf.level)/2)
        n <- p * (1 - p) * (z/b)^2
        n <- ceiling(n)
        cat("\n")
        cat("The required sample size (n) to estimate the population",
            "\n")
        cat("proportion of successes with a", conf.level,
            "confidence interval", "\n")
        cat("so that the margin of error is no more than", b, "is",
            n,".", "\n")
        cat("\n\n")
    }
}

