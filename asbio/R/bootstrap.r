bootstrap <- function (data, statistic, R = 1000, prob = NULL, matrix = FALSE) 
{
    
    theta.hat <- statistic(data)
    if (matrix == TRUE)  n <- nrow(data)
    if (matrix == FALSE) n <- nrow(as.matrix(data))
    
    
    
    if (matrix == TRUE) {
	boot.stat <- matrix(ncol = 1, nrow = R)
     samples <- matrix(ncol = n, nrow = R)
        for (i in 1:R) {
            samp <- data[sample(seq(1, n), replace = TRUE, prob = prob),]
            samples[i,] <- samp
            boot.stat[i] <- statistic(samples[i,])
        }
    }
    if (matrix == FALSE) {
            samp <- sample(data, R, replace = TRUE, prob = prob)
            samples <- matrix(nrow = R, ncol = n, data = samp, byrow = TRUE)
            boot.stat <- apply(samples, 1, statistic)
    }
    
    theta.hat.B <- mean(as.vector(boot.stat))
    bias <- theta.hat.B - theta.hat
    se <- sd(as.vector(boot.stat))
    res <- t(as.matrix(c(theta.hat, theta.hat.B, bias, se)))
    ends <- c("original", "theta.hat.B", "bias", "SE") 
    out <- list(dist = boot.stat, samples = samples , res = res, head = "Bootstrap summary :", ends = ends, data = data, statistic = statistic)
    class(out) <- "bootstrap"
    out
}


print.bootstrap <- function(x, digits = max(3, getOption("digits")), ...) 
{
    cat("\n")
    cat(x$head, "\n\n")
    rq <- structure(x$res, dimnames = list("",x$ends))
    print(rq, digits = digits, justify = "center")
    cat("\n")
    invisible(x)
}
 