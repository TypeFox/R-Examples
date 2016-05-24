x <- c(1,1,0,4,2,1,3,0,0,2);  table(x)
m <- mean(x); n <- length(x); m
lrtStat <- function(x) {
    n <- length(x); m <- mean(x)  
    2 * ( -n * m + n * m * log(m) + n)  
    }
lrtStat(x)
pval <- 1 - pchisq(lrtStat(x), df=1); pval
