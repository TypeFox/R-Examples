######################################################################
## File: test4.R
## Author: Ray Griner
## Date: November 11, 2011 (111111!)
## Desc: Short examples of predicted interval plots (PIPs)
######################################################################

set.seed(12345)
library(PIPS)

#################################################################################
# EXAMPLE 4: Overriding default graphic options.  These are only examples of the
#  customizations that can be made.  Users should be able to customize almost any
#  feature of the graph.  (Type help(par) in an R session for a list of R graphics
#  parameters.  Any of these parameters can be passed to the plot function.)
#################################################################################
myY<-c(rep(1,times=20),rep(0,times=80),rep(1,times=25),rep(0,times=25))
myGroup<-c(rep('A',100),rep('B',50))

pips <- pred.int(y=myY, group=myGroup, N=c(400,400), data.type="binary", iters=500)

print(pips, pi.count=100)
png("test4.png")
plot(pips,                      # The object containing the pips data needs to be first 
     xlab="My x axis label",    # Change default x label 
     main="My graph for #BY#",  # Remember, #BY# is replaced with "B vs A". Not really useful
                                # here with one group but would be useful if multiple groups 
     xlim=c(.12,.6),            # Change limits of x-axis
     ylab="My y axis label", 
     cex.sub=.8,                # Change size of footnote text
     sub="This example shows how to customize the appearance of graphics elements",   # Footnote
     col.main="blue",           # Change color of main title
     axes=FALSE)                # Do not create default axes.  Will use custom axes with axis statements

###########################################
# Custom axes
###########################################
axis(1, at=seq(from=0,to=1,by=.04),las=1)
axis(2, at=seq(from=0,to=100,by=10), las=1)

###########################################
# Create a legend
###########################################
legend("topright",c("Actual","Predicted"), col=c("red","black"), lwd=c(2,1)) 

###########################################
# Create a couple vertical lines
###########################################
abline(v=c(.3,.4), lty=2, lwd=2)

dev.off()            # Save and write graph when done
