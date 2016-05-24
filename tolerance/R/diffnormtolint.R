diffnormtol.int <- function (x1, x2, var.ratio = NULL, alpha = 0.05, P = 0.99, 
    method = c("HALL", "GK", "RG")) 
{
    method <- match.arg(method)
    n1 <- length(x1)
    n2 <- length(x2)
    x1.bar <- mean(x1)
    x2.bar <- mean(x2)
    s1.2 <- var(x1)
    s2.2 <- var(x2)
    z.p <- qnorm(P)
    if (is.null(var.ratio)){
        if (method == "HALL" | method == "GK") {
            q1 <- (s1.2 * (n2 - 3))/(s2.2 * (n2 - 1))
            q2 <- (s2.2 * (n1 - 3))/(s1.2 * (n1 - 1))
            f1 <- ((n1 - 1) * (q1 + 1)^2)/(q1^2 + (n1 - 1)/(n2 -1))
            f2 <- ((n2 - 1) * (q2 + 1)^2)/(q2^2 + (n2 - 1)/(n1 -1))
            nu1 <- (n1 * (1 + q1))/(q1 + (n1/n2))
            nu2 <- (n2 * (1 + q2))/(q2 + (n2/n1))
            lower <- x1.bar - x2.bar - (suppressWarnings(qt(1 - alpha, f1, ncp = (z.p * sqrt(nu1)))) * sqrt((s1.2 + s2.2)/nu1))
            upper <- x1.bar - x2.bar + (suppressWarnings(qt(1 - alpha, f1, ncp = (z.p * sqrt(nu1)))) * sqrt((s1.2 + s2.2)/nu1))
            if (method == "GK") {
                lower.alt <- x1.bar - x2.bar - (suppressWarnings(qt(1 - alpha, f2, ncp = (z.p * sqrt(nu2)))) * sqrt((s1.2 + s2.2)/nu2))
                upper.alt <- x1.bar - x2.bar + (suppressWarnings(qt(1 - alpha, f2, ncp = (z.p * sqrt(nu2)))) * sqrt((s1.2 + s2.2)/nu2))
                lower <- min(lower, lower.alt)
                upper <- max(upper, upper.alt)
            }
        } 
        else if (method == "RG"){
            q1 <- s1.2/s2.2
            f1 <- ((n1 - 1) * (q1 + 1)^2)/(q1^2 + (n1 - 1)/(n2 -1))
            nu1 <- (n1 * (1 + q1))/(q1 + (n1/n2))
            lower <- x1.bar - x2.bar - (suppressWarnings(qt(1 - alpha, f1, ncp = (z.p * sqrt(nu1)))) * sqrt((s1.2 + s2.2)/nu1))
            upper <- x1.bar - x2.bar + (suppressWarnings(qt(1 - alpha, f1, ncp = (z.p * sqrt(nu1)))) * sqrt((s1.2 + s2.2)/nu1))                    
        }
    } 
    else {
        q1 <- var.ratio
        nu1 <- (n1 * (1 + q1))/(q1 + (n1/n2))
        s.d <- sqrt(((1 + 1/q1) * ((n1 - 1) * s1.2 + (n2 - 1) * q1 * s2.2))/(n1 + n2 - 2))
        lower <- x1.bar - x2.bar - (suppressWarnings(qt(1 - alpha, n1 + n2 - 2, ncp = (z.p * sqrt(nu1)))) * s.d/sqrt(nu1))
        upper <- x1.bar - x2.bar + (suppressWarnings(qt(1 - alpha, n1 + n2 - 2, ncp = (z.p * sqrt(nu1)))) * s.d/sqrt(nu1))
    }
        temp <- data.frame(cbind(alpha, P, x1.bar - x2.bar, lower, upper))
        colnames(temp) <- c("alpha", "P", "diff.bar", "1-sided.lower", 
            "1-sided.upper")
        temp
} 










