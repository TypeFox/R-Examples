######################################################################
## File: test1.R
## Author: Ray Griner
## Date: November 11, 2011 (111111!)
## Desc: Short examples of predicted interval plots (PIPs)
######################################################################

set.seed(12345)
library(PIPS)

#######################################################
# EXAMPLE 1: TWO SAMPLE, BINARY OUTCOME
#######################################################
myY<-c(rep(1,times=20),rep(0,times=80),rep(1,times=25),rep(0,times=25))
myGroup<-c(rep('A',100),rep('B',50))

pips <- pred.int(y=myY, group=myGroup, N=c(400,400), data.type="binary", iters=200)
                 

png("test1.png")
print(pips, pi.count=200)
plot(pips, vline=-.25)
dev.off()            # Save and write graph when done

