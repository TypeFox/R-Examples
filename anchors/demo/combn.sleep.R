#######################################################################
##
## Author  : Jonathan Wand <wand(at)stanford.edu>
## Created : 2006-11-03
##  
##
## MODIFIED:
##   2008-04-04 : JW
##   - updated to anchors 3.0 syntax
##
################# making a bar chart...
cat("Repl from King and Wand (2007) Fig3\n")

library(anchors)
data(sleep)
dim(sleep)

sleep$single  <- (sleep$q1008==1)*1
sleep$married <- (sleep$q1008==2)*1
sleep$height  <- sleep$q1006/100

fo <- list(self = q2080 ~ 1,
           vign = cbind(q2119c,q2103c,q2107c,q2115c,q2109c) ~ 1,
           cpolr= ~ q1001+q1002+q1004+q1010+height+married)


z0 <- anchors( fo, data=sleep, method="order")
summary(z0)

z0 <- anchors( fo, data=sleep, method="C", combn=TRUE)
summary(z0)

# fname <- "fig_entropy_sleep.ps"
# postscript( fname , height=8, width=8)
plot.anchors.rank(z0)
plot(z0, xy = c("minimum","estimated"))
plot(z0, xy = c("minimum","interval"))

## look at pairwise relationships
plot(z0$combn, cex=.8 )


                                        # dev.off()
# system(paste("epstopdf",fname))


