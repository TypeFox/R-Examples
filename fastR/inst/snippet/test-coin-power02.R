prob <- c( seq(0,0.4, by=0.10), 0.45, 0.5, 0.55, seq(0.6,1, by=0.10) )
power <- pbinom(39,100,prob) + 1- pbinom(60,100,prob)
print(cbind(prob,power))
