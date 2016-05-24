prob <- c( seq(0, 1, by=0.01) )
power <- pbinom(39,100,prob) + 1- pbinom(60,100,prob)
myplot <- xyplot(power~prob, type="l", 
            main="Power to detect a biased coin with 100 flips",
            xlab="true probability of heads")
