######################################################################
## File: test2.R 
## Author: Ray Griner
## Date: November 11, 2011 (111111!)
## Desc: Short examples of predicted interval plots (PIPs)
######################################################################

set.seed(12345)
library(PIPS)

#######################################################
# EXAMPLE 2: ONE SAMPLE, BINARY OUTCOME
#######################################################
myY<-c(rep(1,times=20),rep(0,times=80),rep(1,times=25),rep(0,times=25))
myGroup<-c(rep('A',100),rep('B',50))

pips <- pred.int(y=myY, N=400, data.type="binary", iters=100)

png("test2.png")
plot(pips)
print(pips, pi.count=100)
dev.off()            # Save and write graph when done

