tstats <- replicate(1000,
    {
    concrete.lm <- lm(sample(strength) ~ limestone * water, concrete)
    summary(concrete.lm)$coef[3,3]
    }
)
mean( abs(tstats) > abs(summary(concrete.lm1)$coef[3,3]) ) 
