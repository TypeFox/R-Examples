NNclean <- 
function (data, k, distances = NULL, edge.correct = FALSE, wrap = 0.1, 
    convergence = 0.001, plot = FALSE, quiet=TRUE) 
{
#    require(mva)
    data <- as.matrix(data)
    d <- dim(data)[2]
    n <- dim(data)[1]
    if (n > 800) {
        options(object.size = 5e+07)
    }
    if ((d == 2) && (edge.correct == TRUE)) {
        r1 <- diff(range(data[, 1]))
        r2 <- diff(range(data[, 2]))
        tran2 <- matrix(c(rep(0, n), rep(r2, n)), byrow = FALSE, 
            nrow = n)
        tran1 <- matrix(c(rep(r1, n), rep(0, n)), byrow = FALSE, 
            nrow = n)
        aux.dat <- rbind(data + tran1, data - tran1, data + tran2, 
            data - tran2, data + tran1 + tran2, data - tran1 + 
                tran2, data - tran1 - tran2, data + tran1 - tran2)
        aux.dat <- aux.dat[aux.dat[, 1] < (max(data[, 1]) + wrap * 
            r1) & aux.dat[, 1] > (min(data[, 1]) - wrap * r1) & 
            aux.dat[, 2] < (max(data[, 2]) + wrap * r2) & aux.dat[, 
            2] > (min(data[, 2]) - wrap * r2), ]
        full.data <- rbind(data, aux.dat)
    }
    else {
        full.data <- data
    }
    dDk <- function(x, lambda, k, d, alpha.d) {
        (exp(-lambda * alpha.d * x^d) * 2 * (lambda * alpha.d)^k * 
            x^(d * k - 1))/gamma(k)
    }
    if (is.null(distances)) {
        distances <- dist(full.data)
    }
    kthNND <- rep(0, n)
    Labels <- 1:(n - 1)
    kthNND[1] <- sort(distances[Labels])[k]
    Labels[(2):(n - 1)] <- Labels[(2):(n - 1)] + (n - 1 - 1)
    for (i in 2:n) {
        kthNND[i] <- sort(distances[Labels])[k]
        Labels[1:(i - 1)] <- Labels[1:(i - 1)] + 1
        Labels[(i + 1):(n - 1)] <- Labels[(i + 1):(n - 1)] + 
            (n - i - 1)
    }
    kthNND <- kthNND[1:n]
    alpha.d <- (2 * pi^(d/2))/(d * gamma(d/2))
    delta <- rep(0, n)
    delta[kthNND > (min(kthNND) + diff(range(kthNND))/3)] <- 1
    delta[is.na(delta)] <- 0
    p <- 0.5
    lambda1 <- k/(alpha.d * mean((kthNND[delta == 0])^d))
    lambda2 <- k/(alpha.d * mean((kthNND[delta == 1])^d))
    loglik.old <- 0
    loglik.new <- 1
    while (abs((loglik.new - loglik.old)/loglik.new) > convergence) {
        delta <- (p * dDk(kthNND, lambda1, k = k, d = d, alpha.d = alpha.d))/(p * 
            dDk(kthNND, lambda1, k = k, d = d, alpha.d = alpha.d) + 
            (1 - p) * dDk(kthNND, lambda2, k = k, d = d, alpha.d = alpha.d))
        delta[is.na(delta)] <- 0
        p <- sum(delta)/n
        if (p>0)
          lambda1 <- (k * sum(delta))/(alpha.d * sum((kthNND^d) * 
             delta))
        if (p<1)
          lambda2 <- (k * sum((1 - delta)))/(alpha.d * sum((kthNND^d) * 
              (1 - delta)))
        if (p<=0)
          lambda1 <- lambda2/2
        if (p>=1)
          lambda2 <- 2*lambda1
        loglik.old <- loglik.new
        loglik.new <- sum(-p * lambda1 * alpha.d * ((kthNND^d) * 
            delta) - (1 - p) * lambda2 * alpha.d * ((kthNND^d) * 
            (1 - delta)) + delta * k * log(lambda1 * alpha.d) + 
            (1 - delta) * k * log(lambda2 * alpha.d))
        if (!quiet){
          cat("l.new=",loglik.new," p=",p," l1=",lambda1,
              " l2=",lambda2," ad=",alpha.d," knnd=", kthNND," delta=",delta,
              "\n")
        }
    }
    if (plot) {
        hist(kthNND, nclass = 20, axes = TRUE, ylab = "Estimate of Mixture", 
            xlim = c(0, max(kthNND)), probability = TRUE, xlab = paste("Distance to", 
                eval(k), "th nearest neighbour"))
        box()
        support <- seq(0, max(kthNND), length = 200)
        lines(support, (p * dDk(support, lambda1, k = k, d = d, 
            alpha.d = alpha.d) + (1 - p) * dDk(support, lambda2, 
            k = k, d = d, alpha.d = alpha.d)))
    }
    probs <- dDk(kthNND, lambda1, k = k, d = d, alpha.d = alpha.d)/(dDk(kthNND, 
        lambda1, k = k, d = d, alpha.d = alpha.d) + dDk(kthNND, 
        lambda2, k = k, d = d, alpha.d = alpha.d))
    probs[is.na(probs)] <- 1
    out <- list(z = round(probs), probs = probs, k = k, lambda1 = lambda1, 
        lambda2 = lambda2, p = p, kthNND = kthNND)
    class(out) <- "nnclean"
    return(out)
}

print.nnclean <- function(x, ...){
  cat("Nearest neighbor noise detection\n")
  cat("by Byers, S. and Raftery, A. E. (1998) Nearest-Neighbor Clutter\n")
  cat("Removal for Estimating Features in Spatial Point Processes,\n")
  cat("Journal of the American Statistical Association, 93, 577-584.\n")
  cat("Classification: ( 0 means noise)\n", x$z, "\n")
  cat("The Poisson process mixture of the distance to kth nearest neighbor\n")
  cat("is characterized by k=",x$k,", p=",x$p,"\n")
  cat("lambda1=",x$lambda1,", lambda2=",x$lambda2,"\n")
}
  
