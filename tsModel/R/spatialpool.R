distmat <- function(x, y, phi) {
        dmat <- dist(cbind(x, y))
        dmat <- as.matrix(dmat)
        exp(-phi * abs(dmat))
}


spatialgibbs <- function(b, v, x, y, phi = 0.1, scale = 1, maxiter = 1000,
                         burn = 500, a0 = 10, b0 = 100000) {
        d <- 1 * scale^2
        D <- diag(d, length(b))
	mu <- 0
        sigma2 <- mean(v)
        
	if(maxiter < burn)
		burn <- 0
        results <- matrix(nrow = maxiter, ncol = 2)
        H <- distmat(x, y, phi)
        R <- chol(H)
        Hinv <- solve(H)
        V <- diag(v)
        I <- diag(1, length(b))
	
	for(i in seq_len(maxiter)) {
                message(i)
                Ainv <- solve(V + sigma2 * H)
                sumAinv <- sum(Ainv)
                
		## Sample mu
                m <- (1/(sumAinv + 1/d)) * sum(Ainv %*% b)
                s2 <- 1/(sumAinv + 1/d)
                mu <- rnorm(1, m, sqrt(s2))
                
		## Sample sigma
                fullc <- makeSigma2FC(b, V, H, mu, a0, b0)
                prop <- 1/rgamma(1, shape = a0, scale = b0)
                
                lr <- fullc(prop) - fullc(sigma2)
                u <- log(runif(1))

                if(is.finite(lr) && !is.na(lr) && u < lr)
                        sigma2 <- prop

                results[i, ] <- c(mu, sigma2)
	}
	if(burn > 0)
		results <- results[-seq_len(burn), ]
        results
}

makeSigma2FC <- function(b, V, H, mu, a0, b0) {
        function(sigma2) {
                (ldmvnorm(b, rep(mu, length(b)), V + sigma2 * H))
                 ## -(a0 + 1) * log(sigma2) - b0 / sigma2)
        }
}
                         
ldmvnorm <- function(x, mean, sigma) {
        if(is.vector(x))
                x <- matrix(x, ncol = length(x))

        sigmaI <- solve(qr(sigma, LAPACK = TRUE))
        distval <- mahalanobis(x, center = mean, cov = sigmaI, inverted = TRUE)
        logdet <- as.numeric(determinant(sigma, logarithm = TRUE)$modulus)
        logretval <- -(ncol(x) * log(2 * pi) + logdet + distval) / 2
        logretval
}

