simData <- function(n=50,m=2,rr=1.5,rseed=123) {

# returns data frame with n rows and m+2 columns x,y,z1,...,zm 
# x,y are 0/1 responses, z's are numeric N(0,1) distributed covariates
# n = number of observations
# m = the number of covariates 
# rr = the (common) OR calculated from P(X,Y|Z)
# rseed : seed for random number generator

p22.vf <- function(p1,p2,rr) {
# p1,p2 are marginals in 2x2 table with OR equal to rr 
# returns p22 in this table
  if(rr==1) p22 <- p1*p2 else {
    d <- ((rr-1)*(p1+p2)+1)^2 - 4*rr*(rr-1)*p1*p2
    p22 <- ((p1+p2)*(rr-1)+1-sqrt(d))/2/(rr-1)
  }
  p22
} # end of p22.vf

  set.seed(rseed)
  ddd <- matrix(NA, nrow=n, ncol=m+2)
  ddd[,3:(2+m)] <- round(rnorm(n*m,0,1),3)
  b1 <- runif(m,-1,1)
  b2 <- runif(m,-1,1)
  # calculate marginals
  p1 <- 1/(1+exp( ddd[,3:(2+m)] %*% b1 ))  # P(X=1)
  p2 <- 1/(1+exp( ddd[,3:(2+m)] %*% b2 ))  # P(Y=1)
  # calculate the 2x2 probs for each subjects and sample it
  for (i in 1:n) {
    p22 <- p22.vf(p1[i],p2[i],rr)  
    p12 <- p1[i]-p22
    p21 <- p2[i]-p22
    p11 <- 1-p22-p12-p21
    xy <- sample(1:4,1,prob=c(p11,p12,p21,p22) )
    x <- ifelse((xy==1 | xy==2), 0,1)
    y <- ifelse((xy==1 | xy==3), 0,1)
    ddd[i,1:2] <- c(x,y)
  }
  ddd <- data.frame(ddd)
  names(ddd) <- c("x","y",paste("z", 1:m, sep=""))
  return(ddd)
} # end simData.vf

# ddd <- simData()

