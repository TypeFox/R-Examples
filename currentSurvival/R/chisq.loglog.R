
chisq.loglog <- function(st,ot.sq) 
{

A <- matrix(log(-log(st[1]))-log(-log(st[2:length(st)])),1,(length(st)-1)) # a vector with differences of complementary log-log transformed survival estimates of the 2nd, 3rd, ... groups from the complementary log-log transformed survival estimate of the first group
SIGMA <- array(ot.sq[1]/((log(st[1]))^2),c(length(st)-1,length(st)-1)) # a matrix containing estimated variances and their sums
diag(SIGMA) <- (ot.sq[1]/((log(st[1]))^2)) + (ot.sq[2:length(st)]/((log(st[2:length(st)]))^2))
chisq <- A %*% solve(SIGMA) %*% t(A) # chi-square test statistic

chisq.loglog <- chisq

}