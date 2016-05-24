"pbb" <-
function(map="bss.background", add.plot=FALSE, ...)
{
# Plots the background map of the BSS project area
#
# add.plot ... if FALSE (default) create new plot
#
# compute range of plotting area:
all=get(map)
xrange=c(min(all[,1],na.rm=TRUE),max(all[,1],na.rm=TRUE))
yrange=c(min(all[,2],na.rm=TRUE),max(all[,2],na.rm=TRUE))
# setup plot
if (!add.plot){plot(1,1,xlim=xrange,ylim=yrange,xlab="",ylab="", ...)}
lines(all,...)
}

