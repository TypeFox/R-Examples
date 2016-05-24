axisIntervals <- function(side=1, atTick1=fda::monthBegin.5, atTick2=fda::monthEnd.5,
              atLabels=fda::monthMid, labels=month.abb, cex.axis=0.9, ...)
{
#  Here's something trivial ...
# 1.  Interval start
  axis(side, at=atTick1, labels=FALSE, ...)
# 2.  Interval end
  if(any(!is.na(atTick2)))axis(side, at=atTick2, labels=FALSE, ...)
# 3.  Interval labels
  axis(side, at=atLabels, labels=labels, tick=FALSE,
       cex.axis=cex.axis, ...)
}

axesIntervals <- function(side=1:2, atTick1=fda::monthBegin.5,
                          atTick2=fda::monthEnd.5, atLabels=fda::monthMid,
                          labels=month.abb, cex.axis=0.9, las=1, ...)
{
  axisIntervals(side[1], atTick1=atTick1, atTick2=atTick2,
                labels=labels, cex.axis=cex.axis, las=las, ...)
  axis(side[2], las=las, ...)
}
