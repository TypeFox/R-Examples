fstats <- replicate(5000,
    {
    concrete.lm <- lm(sample(strength) ~ limestone * water, concrete)
    summary(concrete.lm)$fstat[1]
    }
)
mean( fstats > summary(concrete.lm1)$fstat[1] ) 
