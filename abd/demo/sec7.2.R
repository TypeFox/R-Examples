# Table 7.2-1
dbinom(10:25, 25, 0.061)

# Pr[ number of successes >= 10 ]
sum( dbinom(10:25, 25, 0.061) )
pbinom(9, 25, 0.061, lower.tail=FALSE) 
1 - pbinom(9, 25, 0.061)

binom.test(10, 25, p=0.061)
