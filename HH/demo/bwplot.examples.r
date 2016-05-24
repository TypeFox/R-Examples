## This set of examples shows how to use control the spacing and
## colors of boxplots using the panel.bwplot.intermediate.hh and
## position functions in the HH package.


require(HH)  ## needed for both panel.bwplot.intermediate.hh and position

## 1. Basic example, twelve months and three groups within each month

## Critical items to notice:
##
## a. It is necessary to declare the months to be ordered; by default they are sorted alphabetically.
##
## b. The position() function in the HH package gives control of the spacing of factor levels
##    in an axis.
##
## c. The default axis labels are overprinting.  By displaying only the 2, 5, ... label, and by
##    grouping the three group levels within each month, you get a more readable graph.


mycolors <- c("lightgreen", "darkgreen", "red")

tmp <- data.frame(y=rnorm(612),
                  month=rep(ordered(month.abb, levels=month.abb), 51),
                  group=rep(ordered(mycolors, levels=mycolors), 3, each=17))
gm <- with(tmp, interaction(group, month))

bwplot(y ~ month | group, data=tmp, main="first: three panels with 12 boxes each", layout=c(3,1))
bwplot(y ~ gm, data=tmp, col=mycolors, main="second: one panel with 3*12 boxes")
bwplot(y ~ gm, data=tmp, col=mycolors, main="third: color control of boxes within months",
       panel=panel.bwplot.intermediate.hh)

bwplot(y ~ gm, data=tmp, col=mycolors, main="fourth: key and x-labels by month",
       panel=panel.bwplot.intermediate.hh,
       key=list(
         text=list(mycolors, col=mycolors),
         space="right", border=TRUE, title="Group Name", cex.title=.9),
       scales=list(x=list(at=(3*(1:12)-1), labels=month.abb)))

position(gm) <- as.vector(outer(c(-.8,0,.8), (3*(1:12)-1), "+"))  ## help("position", package="HH")
bwplot(y ~ gm, data=tmp, col=mycolors, panel=panel.bwplot.intermediate.hh,
       key=list(
         text=list(mycolors, col=mycolors),
         space="right", border=TRUE, title="Group Name",cex.title=.9),
       scales=list(x=list(at=(3*(1:12)-1), labels=month.abb)),
       main="fifth: The intended plot: boxes spaced closer within each month")



## 2. boxplots coded by week.  One treatment.

## simulated data
tmp2 <- data.frame(Y=rnorm(40, rep(c(20,25,15,22), 10), 5),
                  week=ordered(rep(1:4, 10)))
position(tmp2$week) <- c(1, 2, 4, 8)

bwplot(Y ~ week, horizontal=FALSE,
       scales=list(x=list(limits=c(0,9),
                     at=position(tmp2$week),
                     labels=position(tmp2$week))),
       data=tmp2, panel=panel.bwplot.intermediate.hh,
       xlab="Week",
       main="Only one treatment")



## 3. boxplots coded by week.  Two treatments.

data(lft.asat) ## Liver Function Tests: ASAT.  Real data.
lft.asat$week <- factor(lft.asat$week)

bwplot(asat ~ week, data=lft.asat,
       panel=panel.bwplot.intermediate.hh,
       xlab="week",
       main=list("Distribution of ASAT by Time: 0", cex=1.4))

position(lft.asat$week) <- as.numeric(levels(lft.asat$week))

bwplot(asat ~ week, data=lft.asat,
       panel=panel.bwplot.intermediate.hh,
       scales=list(
         x=list(
           at=position(lft.asat$week),  ## placement of tick labels and marks
           limits=c(-2,25),             ## x limits
           tck=1)),                     ## draw tick marks
       xlab="Week",
       main=list("Distribution of ASAT by Time: 1", cex=1.4))

lft.asat$wt <- interaction.positioned(lft.asat$week, lft.asat$trt,
                                      b.scale=.5) ## not a method, ".positioned" is necessary

bwplot(asat ~ wt, data=lft.asat,
       panel=panel.bwplot.intermediate.hh,
       scales=list(
          x=list(
            at=position(lft.asat$week), ## placement of tick labels and marks
           limits=c(-2,25),             ## x limits
           tck=1)),                     ## draw tick marks
       xlab="Week",
       main=list("Distribution of ASAT by Time and Treatment: 2", cex=1.4),
       pch=rep(c(17, 16), 7),
       col=rep(c("blue","red"), 7),
       key=list(
         text=list(c("A","B"), col=c("blue","red")),
         points=list(pch=c(17, 16), col=c("blue","red")),
       space="top", columns=2, border=TRUE,
       title="Treatment", cex.title=.9
         )
       )
