epi.herdtest <- function(se, sp, P, N, n, k){
   
   # Probability of testing positive:
   APpos <- P * se + (1 - P) * (1 - sp)
   APneg <- (1 - sp)
 
   # if(n/N < 0.2){
    # Binomial distribution:
    # HSe <- 1 - pbinom(k - 1, n, P)
    # HSp <- phyper(k - 1, N * APneg, N - N * APneg, n)
    # rval <- list(APpos = APpos, APneg = APneg, HSe = HSe, HSp = HSp)
    # }

   # else if(n/N >= 0.2){
    # Hypergeometric distribution:
    HSe <- 1 - phyper(k - 1, N * APpos, N - N * APpos, n)
    HSp <- phyper(k - 1, N * APneg, N - N * APneg, n)
    rval <- list(APpos = APpos, APneg = APneg, HSe = HSe, HSp = HSp)
    # }
rval
}

 