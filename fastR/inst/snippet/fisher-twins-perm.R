
numSims <- 20000
ft <- data.frame(
        twin=rep(c("Di","Mono"),times=c(17,13)),
        conviction=rep(c("No","Yes","No","Yes"),times=c(15,2,3,10))
)
# check to see that table matches
xtabs(~twin+conviction,ft);
# test statistic is value in top left cell
xtabs(~twin+conviction,ft)[1,1];

# simulated data sets
testStats <- replicate(numSims, {
    sft <- ft;
    sft$conviction <- sample(ft$conviction); 
    xtabs(~twin+conviction,sft)[1,1];
    });

# for p-value 
table(testStats);
# tail probabilities
sum(testStats >= 15)/ numSims;
sum(testStats <= 5)/ numSims;
# 2-sided p-value
sum(testStats >= 15 | testStats <= 5)/ numSims ;
