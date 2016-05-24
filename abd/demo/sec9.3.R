## Table 9.3-1
xtabs(~ infection.status + eaten, Trematodes)

## Chi-squared Contingency Test
chisq.test( xtabs(~ infection.status + eaten, Trematodes) )
summary(chisq.test( xtabs(~ infection.status + eaten, Trematodes) ) )
