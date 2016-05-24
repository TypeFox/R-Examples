xtabs(~Student+Parents,data=familySmoking) -> smokeTab; 
smokeTab
chisq.test(smokeTab)
observedStat <- chisq.test(smokeTab)$stat
stats <- replicate(2000,
    {
    chisq.test( xtabs(~sample(Student)+Parents,data=familySmoking))$stat
    }
)
sum( stats > observedStat ) -> x; x / length(stats)   # p-value
binom.test(x,length(stats),alternative="less")$conf.int
