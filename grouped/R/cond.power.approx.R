"cond.power.approx" <-
function(X, theta, sigma, m, a, grouping.mech = c("rounding", "equispaced")){
    n <- nrow(X)
    cc <- ncol(X)
    df. <- n - (cc - 1)
    mu <- X %*% theta
    if(!missing(grouping.mech)) grouping.mech <- match.arg(grouping.mech)        
    gmech <- switch(grouping.mech, rounding = 1, equispaced = 2)
    int. <- if(gmech == 1) rounding(0:m, m) else if(gmech == 2) equispaced(0:m, m)
    log.int <- qlogis(int.)
    w <- ww <- vector(length = n, mode = "list")
    mat1 <- matrix(0, ncol = nrow(int.), nrow = c(cc))
    for(i in 1:n){
        mu.i <- mu[i]
        ksi.m <- (qlogis(rowSums(int.) / 2) - mu.i)
        p.a <- pnorm(log.int[, 1], mean = mu.i, sd = sigma)
        p.b <- pnorm(log.int[, 2], mean = mu.i, sd = sigma)
        s.int1 <- (2 / sigma^3) * ksi.m * (p.b - p.a)
        s.int2 <- (3 / sigma^4) * ksi.m^2 * (p.b - p.a)
        for(j in 1:nrow(int.)){
            mat1[, j] <- c(s.int1[j] * c(X[i, ]))
        }
        w[[i]] <- rowSums(mat1)
        ww[[i]] <- -(1 / (sigma^2)) + sum(s.int2)
    }
    w. <- matrix(colSums(matrix(unlist(w), byrow = TRUE, nrow = n)), byrow = FALSE, ncol = cc)
    ww. <- sum(unlist(ww))
    I. <- rbind((1 / sigma^2) * crossprod(X), w.)
    I. <- cbind(I., c(w., ww.))
    se.delta <- sqrt(solve(I.)[2, 2])
    ncp. <- theta[2] / se.delta
    power. <- 1 - pt(qt(1 - a / 2, df = df.), df = df., ncp = ncp.) + pt(qt(a / 2, df = df.), df = df., ncp = ncp.)
    power.
}

