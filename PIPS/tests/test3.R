######################################################################
## File: test3.R
## Author: Ray Griner
## Date: November 11, 2011 (111111!)
## Desc: Short examples of predicted interval plots (PIPs)
######################################################################

set.seed(12345)
library(PIPS)

#######################################################
# EXAMPLE 3: TWO SAMPLE, CONTINUOUS OUTCOME
#######################################################
myY<-c(rnorm(n=100,mean=.4),rnorm(n=50,mean=.1))
myGroup<-c(rep('A',100),rep('B',50))

pips <- pred.int(y=myY, group=myGroup, N=c(400,400), data.type="t.test", iters=250)

png("test3.png")
plot(pips)
print(pips, pi.count=250)
dev.off()            # Save and write graph when done
