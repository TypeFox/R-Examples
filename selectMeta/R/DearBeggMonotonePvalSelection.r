DearBeggMonotonePvalSelection <- function(y, u, theta0, sigma0, lam = 2, M = 1000, maxiter = 1000, test.stat = function(x){return(min(x))}){

cat("Now computing p-value to assess H_0 of no selection.\nMay take a moment (up to a few hours!),\ndepending on resources and number of studies involved.\n")

n <- length(y)
k <- 1 + floor(n / 2)

res.mono <- matrix(NA, ncol = M, nrow = k + 2)
colnames(res.mono) <- paste("loop", 1:M, sep = "")
rownames(res.mono) <- c(paste("w", 1:k, sep = ""), "theta", "sigma")
ran.num <- matrix(NA, ncol = M, nrow = n)

for (i in 1:M){

    y_i <- rep(NA, n)
    
    ## generate a sample of p-values according to the p-value density
    seed0 <- (i + 10) ^ 2    
    for (r in 1:n){
        rand <- rPval(n = 1, u[r], theta0, sigma0 ^ 2, seed = seed0 + r)
        y_i[r] <- rand$yn
    }

    ran.num[, i] <- y_i    
    set.seed(seed0)
    res.est <- DearBeggMonotone(y = y_i, u = u, lam = lam, maxiter = maxiter, CR = 0.9, NP = NA, trace = FALSE)
    res.mono[, i] <- c(res.est$w, res.est$theta, res.est$sigma)

    out <- paste("Run ", i, " of ", M, " runs done", sep = "")
    print(out)
}

## original estimate
est0 <- DearBeggMonotone(y = y, u = u, lam = lam, maxiter = maxiter, CR = 1, trace = FALSE)

## compute p-value
T0 <- test.stat(est0$w)
Ti <- apply(res.mono[1:k, ], 2, test.stat)
Ti <- Ti[is.na(Ti) == FALSE]
pval <- (1 + sum(Ti <= T0)) / (1 + length(Ti))

res <- list("pval" = pval, "res.mono" = res.mono, "mono0" = est0$w, "Ti" = Ti, "T0" = T0, "ran.num" = ran.num)
return(res)
}







