
chisq.naive <- function(st,ot.sq) 
{

A <- matrix(st[1]-st[2:length(st)],1,(length(st)-1)) # a vector with differences of survival estimates of the 2nd, 3rd, ... groups from survival estimates of the first group 
SIGMA <- array((st[1]^2)*ot.sq[1],c(length(st)-1,length(st)-1)) # a matrix containing estimated variances and their sums
diag(SIGMA) <- (st[1]^2)*ot.sq[1] + (st[2:length(st)]^2)*ot.sq[2:length(st)]
chisq <- A %*% solve(SIGMA) %*% t(A) # chi-square test statistic

chisq.naive <- chisq

}