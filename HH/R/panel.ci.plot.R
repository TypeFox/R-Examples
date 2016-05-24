"panel.ci.plot" <-
function(x, y, newdata, newfit=newfit, ...) {
  tpgsl <- trellis.par.get("superpose.line")
  panel.xyplot(x, y, ...)
  panel.xyplot(x=newdata[[1]],
               y=newfit$fit,
               type="l", lty=tpgsl$lty[2], col=tpgsl$col[2], lwd=tpgsl$lwd[2])
  panel.xyplot(x=newdata[[1]],
               y=newfit$ci.fit[,"lower"],
               type="l", lty=tpgsl$lty[3], col=tpgsl$col[3], lwd=tpgsl$lwd[3])
  panel.xyplot(x=newdata[[1]],
               y=newfit$ci.fit[,"upper"],
               type="l", lty=tpgsl$lty[3], col=tpgsl$col[3], lwd=tpgsl$lwd[3])
  panel.xyplot(x=newdata[[1]],
               y=newfit$pi.fit[,"lower"],
               type="l", lty=tpgsl$lty[4], col=tpgsl$col[4], lwd=tpgsl$lwd[4])
  panel.xyplot(x=newdata[[1]],
               y=newfit$pi.fit[,"upper"],
               type="l", lty=tpgsl$lty[4], col=tpgsl$col[4], lwd=tpgsl$lwd[4])
  if.R(s={
    axis(1, at=mean(x), tck=-.03, labels=FALSE)
    axis(1, at=mean(x), ticks=FALSE,  labels="xbar", line=.9)
    axis(3, at=mean(x), tck=-.03, labels=FALSE)
    axis(3, at=mean(x), ticks=FALSE,  labels="xbar")
  },r={
    cpl <- current.panel.limits()
    if (mean(x) >= cpl$xlim[1] && mean(x) <= cpl$xlim[2]) {
      pushViewport(viewport(xscale = cpl$xlim,
                            yscale = cpl$ylim,
                            clip = "off"))
      ## put anything you want unclipped inside this:
      ## panel.abline(v=mean(x), lty=2)
      panel.axis("bottom", at=mean(x), tck=1.5, labels=FALSE)
      panel.axis("bottom", at=mean(x), ticks=FALSE,  labels=expression(bar(x)), rot=0, outside=TRUE)
      panel.axis("top", at=mean(x), tck=1.5, labels=FALSE)
      panel.axis("top", at=mean(x), ticks=FALSE,  labels=expression(bar(x)), rot=0, outside=TRUE)
      ## end of unclipped part
      upViewport()
    }
  })
}

## source("c:/HOME/rmh/HH-R.package/HH/R/panel.ci.plot.R")
