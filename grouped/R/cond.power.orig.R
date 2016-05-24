"cond.power.orig" <-
function(X, theta, sigma, m, a, grouping.mech = c("rounding", "equispaced")){
    n <- nrow(X)
    cc <- ncol(X)
    df. <- n - (cc - 1)
    mu <- X %*% theta
    g.mech <- match.arg(grouping.mech)       
    gmech <- switch(g.mech, rounding = 1, equispaced = 2)
    int. <- if(gmech == 1) rounding(0:m, m) else if(gmech == 2) equispaced(0:m, m)
    int. <- eval(int.)                     
    log.int <- qlogis(int.)
    qa <- ifelse(!is.finite(log.int[, 1]), 0, log.int[, 1])
    qb <- ifelse(!is.finite(log.int[, 2]), 0, log.int[, 2]) 
    w <- ww <- www <- vector(length = n, mode = "list")
    mat1 <- matrix(0, ncol = nrow(int.), nrow = c(cc * cc))
    mat2 <- matrix(0, ncol = nrow(int.), nrow = c(cc))
    mat3 <- matrix(0, ncol = nrow(int.), nrow = 1)
    for(i in 1:n){
        mu.i <- mu[i]
        qa.m <- (qa - mu.i)
        qb.m <- (qb - mu.i)
        f.a <- dnorm(log.int[, 1], mean = mu.i, sd = sigma)
        f.b <- dnorm(log.int[, 2], mean = mu.i, sd = sigma)
        p.a <- pnorm(log.int[, 1], mean = mu.i, sd = sigma)
        p.b <- pnorm(log.int[, 2], mean = mu.i, sd = sigma)
        s.int1 <-  ((1 / (sigma^2)) * (qb.m * f.b - qa.m * f.a) + (f.b - f.a)^2 / (p.b - p.a))
        s.int2 <- -(1 / sigma) * (f.b - f.a) + (1 / (sigma^3)) * ((qb.m^2) *f.b - (qa.m^2) * f.a) +
                        (1 / sigma) * (qb.m * f.b - qa.m * f.a) * (f.b - f.a) / (p.b - p.a)
        s.int3 <- -(2 / sigma^2) * (qb.m * f.b - qa.m * f.a) + (1 / (sigma^4)) * ((qb.m^3) * f.b -
                        (qa.m^3) * f.a) + (1 / sigma^2) * (qb.m * f.b - qa.m * f.a)^2 / (p.b - p.a)
        for(j in 1:nrow(int.)){
            mat1[, j] <- c(s.int1[j] * (c(X[i, ]) %*% t(X[i, ])))
            mat2[, j] <- c(s.int2[j] * (c(X[i, ])))
            mat3[, j] <- s.int3[j]
        }
        w[[i]] <- matrix(rowSums(mat1), byrow = FALSE, ncol = cc)                 
        ww[[i]] <- rowSums(mat2)
        www[[i]] <- rowSums(mat3)
    }
    w. <- matrix(colSums(matrix(unlist(w), byrow = TRUE, nrow = n)), byrow = FALSE, ncol = cc)
    ww. <- colSums(matrix(unlist(ww), byrow = TRUE, nrow = n))
    www. <- sum(unlist(www))
    I. <- rbind(w., t(ww.))
    I. <- cbind(I., c(ww., www.))
    se.delta <- sqrt(solve(I.)[2, 2])
    ncp. <- theta[2] / se.delta
    power. <- 1 - pt(qt( 1 - a / 2, df = df.), df = df., ncp = ncp.) + pt(qt(a / 2, df = df.), df = df., ncp = ncp.)
    power.
}

