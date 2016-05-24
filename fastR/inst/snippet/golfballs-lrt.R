# LRT calculation
O <- golfballs; O
E <- rep(486/4,4); E
G <- 2 * sum ( O * log(O/E)); G       # lrt Goodness of fit statistic
1 - pchisq(G,df=3)
