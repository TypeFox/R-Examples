"pkb" <-
function(map="kola.background",which.map=c(1,2,3,4),map.col=c(5,1,3,4),map.lwd=c(2,1,2,1),
  add.plot=FALSE, ...)
{
# Plots the background map of the KOLA project area
#
# which==1 ... plot project boundary
# which==2 ... plot coast line
# which==3 ... plot country borders
# which==4 ... plot lakes and rivers
# add.plot ... if FALSE (default) create new plot
#
# compute range of plotting area:
all=get(map)
xrange=c(min(all[[1]][,1],na.rm=TRUE),max(all[[1]][,1],na.rm=TRUE))
yrange=c(min(all[[1]][,2],na.rm=TRUE),max(all[[1]][,2],na.rm=TRUE))
# setup plot
if (!add.plot){plot(1,1,xlim=xrange,ylim=yrange,xlab="",ylab="", ...)}
# plot desired lines:
for (i in 1:length(which.map)){
  lines(all[[which.map[i]]],col=map.col[i],lwd=map.lwd[i], ...)
}
}

