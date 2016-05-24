######################################################################
# Section 3.2 Normal Distribution with Known Mean but Unknown Variance
######################################################################
 
 library(LearnBayes)

 data(footballscores)
 attach(footballscores)
 d = favorite - underdog - spread
 n = length(d)
 v = sum(d^2)

 P = rchisq(1000, n)/v
 s = sqrt(1/P)
 hist(s)

 quantile(s, probs = c(0.025, 0.5, 0.975))
