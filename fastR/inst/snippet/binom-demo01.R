randomData <- rbinom(n=20,size=4,prob=0.5)
randomData
table(randomData)
dbinom(0:4,size=4,prob=0.5)          # matches earlier example 
dbinom(0:4,size=4,prob=0.5) * 20     # pretty close to our table above
pbinom(0:4,size=4,prob=0.5)          # same as cumsum(dbinom(...))
