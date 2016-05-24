## last modified June 2002

groupstats <- function(mixdat) 
{
    m <- nrow(mixdat)
    n <- sum(mixdat[, 2])
    amid <- (c(3 * mixdat[1, 1] - 2 * mixdat[2, 1], mixdat[-m, 
        1]) + c(mixdat[-m, 1], 3 * mixdat[m - 1, 1] - 2 * mixdat[m - 
        2, 1]))/2
    amid[1] <- ifelse(amid[1] <= 0 & mixdat[1, 1] > 0, mixdat[1, 
        1]/2, amid[1])
    data.frame(pi = 1, mu = sum(amid * mixdat[, 2])/n, sigma = sqrt(var(rep(amid, 
        mixdat[, 2]))))
}
