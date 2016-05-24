cor.comb.rep<-function (x1, x2, x3, x4, set1, set2, set3, set4, sims = 100, CI = 0.95)
{
    i=0
    
    
    LL <- (1 - CI)/2
    UL <- 1 - LL
    split.r <- rep(NA, sims)
    N <- nrow(data.frame(x1))
    half.n <- N/2
    for (j in 1:sims) {
        s1 <- sample(1:N, N, FALSE)
        r1 <- foreach(i = 1:ncol(set1), .combine = "rbind") %do%
            cor.comb(x1[s1 <= half.n], x2[s1 <= half.n], x3[s1 <=
                half.n], x4[s1 <= half.n], set1[s1 <= half.n, i],
                set2[s1 <= half.n, i], set3[s1 <= half.n, i], set4[s1 <=
                  half.n, i], simple = TRUE, seed = FALSE)
        r2 <- foreach(i = 1:ncol(set1), .combine = "rbind") %do%
            cor.comb(x1[s1 > half.n], x2[s1 > half.n], x3[s1 >
                half.n], x4[s1 > half.n], set1[s1 > half.n, i],
                set2[s1 > half.n, i], set3[s1 > half.n, i], set4[s1 >
                  half.n, i], simple = TRUE, seed = FALSE)
        split.r[j] <- cor(as.vector(r1), as.vector(r2))
    }
    reps <- split.r * 2/(split.r + 1)
    out <- cbind(N, mean(reps), sd(reps), quantile(reps, LL),
        quantile(reps, UL))
    colnames(out) <- c("N", "Rep", "SE", "LL", "UL")
    rownames(out) <- "Results"
    return(out)
}
