umatol.int <- function (x, n = NULL, dist = c("Bin", "NegBin", "Pois"), N, 
    alpha = 0.05, P = 0.99) 
{
    dist <- match.arg(dist)
    if (length(x) == 1 & is.null(n) == TRUE) 
        stop("A value for n must be specified!")
    if (length(x) > 1) 
        n <- length(x)
    y = sum(x)
    if (dist == "Bin") {
	  if(y > 0){
        r.0 <- 1 - qbeta(alpha, N * n - y + 1, y)
        r.1 <- 1 - qbeta(alpha, N * n - y, y + 1)
        R <- max(r.0, r.1)
	   } 
	   else R <- 1-alpha^(1/(N * n))
        f.2 <- function(k, P, N) qbeta(1 - P, k + 1, N - k)
        k <- -1
        temp <- -1
        while (temp < R) {
            k <- k + 1
            temp <- ifelse(k < N, f.2(k = k, P = P, N = N), 1)
        }
        p.hat <- (y/n)/N
        temp.out <- data.frame(cbind(alpha, P, p.hat, k))
        colnames(temp.out) <- c("alpha", "P", "p.hat", "1-sided.upper")
    }
    if (dist == "NegBin") {
	   if(y > 0){
	       r.1 <- 1 - qbeta(alpha, n * N, y + 1)
		  r.0 <- 1 - qbeta(alpha, n * N, y)
		  R <- max(r.0, r.1)
	   } 
	   else R <- 1 - (alpha)^(1/(n*N))
        k <- -1
        temp <- 1.1
        while (temp > 1 - R) {
            k <- k + 1
            temp <- qbeta(P, N, k + 1)
        }
        nu.hat <- N/(N + (y/n))
        temp.out <- data.frame(cbind(alpha, P, nu.hat, k))
        colnames(temp.out) <- c("alpha", "P", "nu.hat", "1-sided.upper")
    }
    if (dist == "Pois") {
	   if(y > 0){
	       r.0 <- qchisq(1 - alpha, 2 * y + 2)/(2 * n)
		  r.1 <- qchisq(1 - alpha, 2 * y)/(2 * n)
		  R <- 2 * max(r.0, r.1)
		} 
	   else R <- -(1/n) * log(alpha)
        k <- -1
        temp <- -1
        while (temp < R) {
            k <- k + 1
            temp <- qchisq(1 - P, 2 * k + 2)
        }
        lambda.hat <- y/n
        temp.out <- data.frame(cbind(alpha, P, lambda.hat, k))
        colnames(temp.out) <- c("alpha", "P", "lambda.hat", "1-sided.upper")
    }
    temp.out
}
