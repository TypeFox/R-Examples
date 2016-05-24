t.test(endurance$Vitamin,endurance$Placebo,paired=TRUE)
t.test(endurance$Vitamin-endurance$Placebo) # same as above
t.test(log(endurance$Vitamin),log(endurance$Placebo),paired=TRUE)
t.test(log(endurance$Vitamin)-log(endurance$Placebo)) # same as above
t.test(log(endurance$Vitamin/endurance$Placebo)) # same as above again
t.test(endurance$Vitamin/endurance$Placebo)
t.test(1/endurance$Vitamin,1/endurance$Placebo,paired=TRUE)
x <- sum(endurance$Vitamin > endurance$Placebo)
n <- nrow(endurance)
binom.test(x,n)
prop.test(x,n)
