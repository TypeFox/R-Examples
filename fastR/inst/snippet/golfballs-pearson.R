# manual calculation
O <- golfballs; O
E <- rep(486/4,4); E
X <- sum ( (O - E)^2 / E); X
1 - pchisq(X,df=3)

# repeated using built-in method
chisq.test(O)
