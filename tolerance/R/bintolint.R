bintol.int <- function (x, n, m = NULL, alpha = 0.05, P = 0.99, side = 1, method = c("LS", 
    "WS", "AC", "JF", "CP", "AS", "LO", "PR", "PO", "CL", "CC", "CWS"), a1 = 0.5, a2 = 0.5) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (side == 2) {
        alpha <- alpha/2
        P <- (P + 1)/2
    }
    method <- match.arg(method)
    if (length(x) > 1) 
        x <- sum(x)
    p.hat <- x/n
    k <- qnorm(1 - alpha)
    x.tilde <- (x + k^2/2)
    n.tilde <- n + k^2
    p.tilde <- x.tilde/n.tilde
    if(is.null(m)) m <- n
    if (method == "LS") {
        lower.p <- p.hat - k * sqrt(p.hat * (1 - p.hat)/n)
        upper.p <- p.hat + k * sqrt(p.hat * (1 - p.hat)/n)
    }
    else if (method == "WS") {
        lower.p <- p.tilde - (k * sqrt(n)/(n + k^2)) * sqrt(p.hat * 
            (1 - p.hat) + (k^2/(4 * n)))
        upper.p <- p.tilde + (k * sqrt(n)/(n + k^2)) * sqrt(p.hat * 
            (1 - p.hat) + (k^2/(4 * n)))
    }
    else if (method == "AC") {
        lower.p <- p.tilde - k * sqrt(p.tilde * (1 - p.tilde)/n.tilde)
        upper.p <- p.tilde + k * sqrt(p.tilde * (1 - p.tilde)/n.tilde)
    }
    else if (method == "JF") {
        if (is.null(a1) | is.null(a2)) 
            stop(paste("Must specify hyperparameters for this option!", 
                "\n"))
        lower.p <- suppressWarnings(qbeta(alpha, shape1 = x + 
            a1, shape2 = n - x + a2))
        upper.p <- suppressWarnings(qbeta(1 - alpha, shape1 = x + 
            a1, shape2 = n - x + a2))
    }
    else if (method == "CP") {
        lower.p <- (1 + ((n - x + 1) * qf(1 - alpha, df1 = 2 * 
            (n - x + 1), df2 = (2 * x)))/x)^(-1)
        upper.p <- (1 + (n - x)/((x + 1) * qf(1 - alpha, df1 = 2 * 
            (x + 1), df2 = 2 * (n - x))))^(-1)
    }
    else if (method == "AS") {
        p.sin <- (x + (3/8))/(n + (3/4))
        lower.p <- (sin(asin(sqrt(p.sin)) - 0.5 * k/sqrt(n)))^2
        upper.p <- (sin(asin(sqrt(p.sin)) + 0.5 * k/sqrt(n)))^2
    }
    else if (method == "LO") {
        l.hat <- log(x/(n - x))
        V.hat <- n/(x * (n - x))
        lower.lambda <- l.hat - k * sqrt(V.hat)
        upper.lambda <- l.hat + k * sqrt(V.hat)
        lower.p <- exp(lower.lambda)/(1 + exp(lower.lambda))
        upper.p <- exp(upper.lambda)/(1 + exp(upper.lambda))
    }
    else if (method == "PR") {
	  z.hat <- qnorm(p.hat)
	  TEMP <- pnorm(z.hat+c(-1,1)*k*sqrt((p.hat*(1-p.hat))/(n*dnorm(z.hat)^2)))
	  lower.p <- TEMP[1]
	  upper.p <- TEMP[2]
    }
    else if (method == "PO") {
	  mu.hat <- (-log(p.hat))
          TEMP <- exp(-(mu.hat + c(-1,1)*k*sqrt((1-p.hat)/(n*p.hat))))
	  lower.p <- TEMP[2]
	  upper.p <- TEMP[1]
    }
    else if (method == "CL") {
	  mu.hat <- -log(p.hat)
	  gamma.hat <- log(mu.hat)
          TEMP <- exp(-exp(gamma.hat + c(-1,1)*k*sqrt((1-p.hat)/(n*p.hat*mu.hat^2))))
	  lower.p <- TEMP[2]
	  upper.p <- TEMP[1]
    }
    if (method == "CC") {
        lower.p <- p.hat - k * sqrt(p.hat * (1 - p.hat)/n) - 1/(2*n)
        upper.p <- p.hat + k * sqrt(p.hat * (1 - p.hat)/n) + 1/(2*n)
    }
    else if (method == "CWS") {
        lower.p <- (2*n*p.hat + k^2 - 1 - k*sqrt(k^2 - 2 - (1/n) + 
		4*p.hat*(n*(1-p.hat) + 1)))/(2*(n + k^2))
        upper.p <- (2*n*p.hat + k^2 + 1 + k*sqrt(k^2 + 2 - (1/n) + 
		4*p.hat*(n*(1-p.hat) - 1)))/(2*(n + k^2))
    }
    lower.p <- max(0,lower.p)
    upper.p <- min(upper.p,1)
    lower <- qbinom(1-P, size = m, prob = lower.p)
    upper <- qbinom(P, size = m, prob = upper.p)    
	if (side == 2) {
        alpha <- 2 * alpha
        P <- (2 * P) - 1
    }
    temp <- data.frame(cbind(alpha, P, p.hat, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "p.hat", "2-sided.lower", 
            "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "p.hat", "1-sided.lower", 
            "1-sided.upper")
    }
    temp
}
