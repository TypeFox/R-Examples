
chisq.log <- function(st,ot.sq) 
{

A <- matrix(log(st[1])-log(st[2:length(st)]),1,(length(st)-1)) # a vector with differences of log transformated survival estimates of the 2nd, 3rd, ... groups from the log transformed survival estimate of the first group 
SIGMA <- array(ot.sq[1],c(length(st)-1,length(st)-1)) # a matrix containing estimated variances and their sums
diag(SIGMA) <- ot.sq[1] + ot.sq[2:length(st)]
chisq <- A %*% solve(SIGMA) %*% t(A) # chi-square test statistic

chisq.log <- chisq

}