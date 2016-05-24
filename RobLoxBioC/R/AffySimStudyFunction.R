###############################################################################
## Function to perform simulation study comparing Tukey's biweight with
## rmx estimators
###############################################################################
AffySimStudy <- function(n, M, eps, seed = 123, eps.lower = 0, eps.upper = 0.05, 
                         steps = 3L, fsCor = TRUE, contD, 
                         plot1 = FALSE, plot2 = FALSE, plot3 = FALSE){
    stopifnot(n >= 3)
    stopifnot(eps >= 0, eps <= 0.5)
    if(plot1){
        from <- min(-6, q(contD)(1e-15))
        to <- max(6, q(contD)(1-1e-15))
        curve(pnorm, from = from, to = to, lwd = 2, n = 201, 
              main = "Comparison: ideal vs. real", ylab = "cdf")
        fun <- function(x) (1-eps)*pnorm(x) + eps*p(contD)(x)
        curve(fun, from = from, to = to, add = TRUE, col = "orange", 
              lwd = 2, n = 201, ylab = "cdf")
        legend("topleft", legend = c("ideal", "real"), 
              fill = c("black", "orange"))
    }

    set.seed(seed)
    r <- rbinom(n*M, prob = eps, size = 1)
    Mid <- rnorm(n*M)
    Mcont <- r(contD)(n*M)
    Mre <- matrix((1-r)*Mid + r*Mcont, ncol = n)
    ind <- rowSums(matrix(r, ncol = n)) >= n/2
    while(any(ind)){
        M1 <- sum(ind)
        cat("Samples to re-simulate:\t", M1, "\n")
        r <- rbinom(n*M1, prob = eps, size = 1)
        Mid <- rnorm(n*M1)
        Mcont <- r(contD)(n*M1)
        Mre[ind,] <- (1-r)*Mid + r*Mcont
        ind[ind] <- rowSums(matrix(r, ncol = n)) >= n/2
    }
    rm(Mid, Mcont, r, ind)


    if(plot2){
        ind <- if(M <= 20) 1:M else sample(1:M, 20)
        if(plot1) dev.new()
        M1 <- min(M, 20)
        print(
          stripplot(rep(1:M1, each = n) ~ as.vector(Mre[ind,]), 
                    ylab = "samples", xlab = "x", pch = 20,
                    main = ifelse(M <= 20, "Samples", "20 randomly chosen samples"))
        )
    }

    ## ML-estimator: mean and sd
    Mean <- rowMeans(Mre)
    Sd <- sqrt(rowMeans((Mre-Mean)^2))
    ## Median and MAD
    Median <- rowMedians(Mre)
    Mad <- rowMedians(abs(Mre - Median))/qnorm(0.75)
    ## Tukey 1-step + MAD
    Tukey <- apply(Mre, 1, function(x) tukey.biweight(x))
    Tukey <- cbind(Tukey, Mad)

    ## Radius-minimax estimator
    RadMinmax <- estimate(rowRoblox(Mre, eps.lower = eps.lower, 
                                    eps.upper = eps.upper, k = steps,
                                    fsCor = fsCor))

    if(plot3){
        Ergebnis1 <- list(Mean, Median, Tukey[,1], RadMinmax[,1])
        Ergebnis2 <- list(Sd, Mad, RadMinmax[,2])
        myCol <- brewer.pal(4, "Dark2")
        if(plot1 || plot2) dev.new()
        layout(matrix(c(1, 1, 1, 1, 3, 2, 2, 2, 2, 3), ncol = 2))
        boxplot(Ergebnis1, col = myCol, pch = 20, main = "Location")
        abline(h = 0)
        boxplot(Ergebnis2, col = myCol[c(1,2,4)], pch = 20, main = "Scale")
        abline(h = 1)
        op <- par(mar = rep(2, 4))
        plot(c(0,1), c(1, 0), type = "n", axes = FALSE)
        legend("center", c("ML", "Med/MAD", "biweight", "rmx"),
               fill = myCol, ncol = 4, cex = 1.5)
        on.exit(par(op))
    }

    ## ML-estimator
    MSE1.1 <- n*mean(Mean^2)
    ## Median + MAD
    MSE2.1 <- n*mean(Median^2)
    ## Tukey
    MSE3.1 <- n*mean(Tukey[,1]^2)
    ## Radius-minimax
    MSE4.1 <- n*mean(RadMinmax[,1]^2)
    empMSE <- data.frame(ML = MSE1.1, Med = MSE2.1, Tukey = MSE3.1, "rmx" = MSE4.1)
    rownames(empMSE) <- "n x empMSE (loc)"
    relMSE <- empMSE[1,]/empMSE[1,4]
    empMSE <- rbind(empMSE, relMSE)
    rownames(empMSE)[2] <- "relMSE (loc)"

    ## ML-estimator
    MSE1.2 <- n*mean((Sd-1)^2)
    ## Median + MAD
    MSE2.2 <- n*mean((Mad-1)^2)
    ## Tukey
    MSE3.2 <- MSE2.2
    ## Radius-minimax
    MSE4.2 <- n*mean((RadMinmax[,2]-1)^2)
    empMSE <- rbind(empMSE, c(MSE1.2, MSE2.2, MSE3.2, MSE4.2))
    rownames(empMSE)[3] <- "n x empMSE (scale)"
    relMSE <- empMSE[3,]/empMSE[3,4]
    empMSE <- rbind(empMSE, relMSE)
    rownames(empMSE)[4] <- "relMSE (scale)"
    empMSE <- rbind(empMSE, c(MSE1.1 + MSE1.2, MSE2.1 + MSE2.2, MSE3.1 + MSE3.2, 
                              MSE4.1 + MSE4.2))
    rownames(empMSE)[5] <- "n x empMSE (loc + scale)"
    relMSE <- empMSE[5,]/empMSE[5,4]
    empMSE <- rbind(empMSE, relMSE)
    rownames(empMSE)[6] <- "relMSE (loc + scale)"

    empMSE
}


