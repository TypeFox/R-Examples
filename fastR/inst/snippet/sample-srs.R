sample(1:1000,25)               # 25 random numbers in 1-1000
sample(1:10,8,replace=TRUE)     # iid random sample
require(vcd)
sample(VonBort$deaths,10)              # show a sample of size 10
mean(sample(VonBort$deaths,10))        # mean of a (different) sample
replicate(10,mean(sample(VonBort$deaths,10))) # do it 10 times
mean(VonBort$deaths)                   # mean of entire data set
