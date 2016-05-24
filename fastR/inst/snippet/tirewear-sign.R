x <- sum(tirewear$weight > tirewear$groove)
n <- length(tirewear$weight)
binom.test(x,n)
prop.test(x,n)
