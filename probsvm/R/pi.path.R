pi.path <- 
function (lambda, x, y, kernel = "linear", kparam = NULL, eps = 1e-05, Nmoves = 5 * n, ridge) 
{
    n <- length(y)
	
	
    part1 <- sub.pipath(lambda, x, y, kernel, kparam, eps, Nmoves, ridge)
	
	
    part2 <- sub.pipath(lambda, x, -y, kernel, kparam, eps, Nmoves, ridge)
	
	
    pi1 <- part1$Pi
    pi2 <- 1 - part2$Pi

    n1 <- length(pi1)
    n2 <- length(pi2)

    alpha01 <- part1$alpha0
    alpha02 <- -part2$alpha0

    alpha1 <- part1$alpha
    alpha2 <- part2$alpha
    alpha1 <- matrix(alpha1, length(alpha1)/length(pi1), length(pi1))
    alpha2 <- matrix(alpha2, length(alpha2)/length(pi2), length(pi2))

    result1 <- matrix(rbind(pi1, alpha01, alpha1), n + 2, n1)
    result2 <- matrix(rbind(pi2, alpha02, alpha2), n + 2, n2)
    result2 <- matrix(result2[, order(result2[1, ])], n + 2, n2)
    result <- cbind(c(0, 1, rep(0, n)), result2[, -n2], result1[, -1], c(1, 1, rep(0, n)))

    pi <- result[1, ]
    alpha0 <- result[2, ]
    alpha <- result[-(1:2), ]
    obj <- list(pi = pi, alpha0 = alpha0, alpha = t(alpha), x = x, 
        y = y, K = part1$K, lambda = lambda, kernel = kernel, kparam = kparam)
    obj
}
