poistol.int <- function (x, n, m = NULL, alpha = 0.05, P = 0.99, side = 1, method = c("TAB", 
    "LS", "SC", "CC", "VS", "RVS", "FT", "CSC")) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    method <- match.arg(method)
    if (length(x) > 1) 
        x <- sum(x)
    if (side == 2) {
        alpha <- alpha/2
        P <- (P + 1)/2
    }
    if(is.null(m)) m <- n
    if (method == "TAB") {
        lower.lambda <- 0.5 * qchisq(alpha, df = (2 * x))/n
        upper.lambda <- 0.5 * qchisq(1 - alpha, df = (2 * x + 
            2))/n
    }
    if (method == "LS") {
        lower.lambda <- (x/n) - (qnorm(1 - alpha) * sqrt(x))/n
        upper.lambda <- (x/n) + (qnorm(1 - alpha) * sqrt(x))/n
    }
    if (method == "SC") {
        k <- qnorm(1 - alpha)
        lower.lambda <- (x/n) + (k^2/(2 * n)) - (k/sqrt(n)) * 
            sqrt((x/n) + (k^2/(4 * n)))
        upper.lambda <- (x/n) + (k^2/(2 * n)) + (k/sqrt(n)) * 
            sqrt((x/n) + (k^2/(4 * n)))
    }
    if (method == "CC") {
        lower.lambda <- (x/n) - (qnorm(1 - alpha) * sqrt(x)/n + 0.5/n)
        upper.lambda <- (x/n) + (qnorm(1 - alpha) * sqrt(x)/n + 0.5/n)
    }
    if (method == "VS") {
        k <- qnorm(1 - alpha)
        lower.lambda <- (x/n) + (k^2/(4 * n)) - (k*sqrt(x)/n)
        upper.lambda <- (x/n) + (k^2/(4 * n)) + (k*sqrt(x)/n)
    }
    if (method == "RVS") {
        k <- qnorm(1 - alpha)
        lower.lambda <- (x/n) + (k^2/(4 * n)) - (k*sqrt((x/n+3/8)/n))
        upper.lambda <- (x/n) + (k^2/(4 * n)) + (k*sqrt((x/n+3/8)/n))
    }
    if (method == "FT") {
	  g <- function(z) ((z^2 - 1)/(2*z))^2
        k <- qnorm(1 - alpha)
	  TEMP.L <- sqrt(x/n) + sqrt((x/n) + 1) - k*(1/sqrt(n))
	  TEMP.U <- sqrt(x/n) + sqrt((x/n) + 1) + k*(1/sqrt(n))
	  if(TEMP.L >= 1){
		lower.lambda <- g(TEMP.L)
	  } else lower.lambda <- 0
        upper.lambda <- g(TEMP.U)
    }
    if (method == "CSC") {
        k <- qnorm(1 - alpha)
	  lam <- x/n
        lower.lambda <- lam - (1/(2*n)) + k^2/(2*n) - sqrt((lam-1/(2*n)+k^2/(2*n))^2-lam^2+lam/n-1/(4*n^2))
        upper.lambda <- lam + (1/(2*n)) + k^2/(2*n) + sqrt((lam+1/(2*n)+k^2/(2*n))^2-lam^2-lam/n-1/(4*n^2))
    }
    lower.lambda <- max(0,lower.lambda)
    lower <- qpois(1-P, lambda = (m*lower.lambda))
    upper <- qpois(P, lambda = (m*upper.lambda))
    if (side == 2) {
        alpha <- 2 * alpha
        P <- (2 * P) - 1
    }
    temp <- data.frame(cbind(alpha, P, x/n, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "lambda.hat", "2-sided.lower", 
            "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "lambda.hat", "1-sided.lower", 
            "1-sided.upper")
    }
    temp
}


