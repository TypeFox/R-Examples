mmc.explain <- function(group, n, ybar, ms.5, crit.point=1.96,
                        ylabel="ylabel",
                        factor.label="factor.label",
                        xlim,
                        ylim=xlim,
                        exit=FALSE,
                        col=c(
                          "gray30",      ## h-v axis
                          "navyblue",        ## m-d axis for data
                          "SlateBlue",   ## perpediculars to m-d
                          "darkgray", ## isomeans grid (xbar-ybar)
                          "royalblue",       ## m-d axis through c(0,0)
                          "red"       ## B-D contrast
                          )) {
  ybar.mat <- matrix(ybar, length(ybar), length(ybar))
  ry <- range(ybar)
  if (missing(xlim)) xlim <- ry + c(-1,1)*.7*diff(ry)
  xyplot(as.vector(ybar.mat) ~ as.vector(t(ybar.mat)),
         aspect=1, xlim=xlim, ylim=ylim,
         xlab=list("h=ybar", col=col[1]),
         ylab=list("v=ybar", col=col[1]),
         main=list(paste("Construction of MMC plot:", ylabel, "~", factor.label),
           col=col[1]),
         par.settings=list(
           clip=list(panel="off"),
           layout.widths=list(left.padding=10),
           layout.heights=list(bottom.padding=10),
           axis.line=list(col=col[1])),
         crit.point=crit.point,
         uy=ybar,
         uy.labels=group,
         ms.5=ms.5,
         exit=exit,         col=col,
         scales=list(col=col[1]),
         panel=function(x, y, ...,
           crit.point, uy, uy.labels, ms.5, exit, col) {
           ## (ybar, ybar) points
           panel.xyplot(x, y, ...,
                        col=col[4])
           ## square in (h,v) and (d,m) coordinates
           panel.segments(min(uy), uy, max(uy), uy, lty=2, col=col[4])
           panel.segments(uy, min(uy), uy, max(uy), lty=2, col=col[4])

           ## labels for constant v and h lines
           panel.text(x=uy, y=rep(min(y)-.5, length(uy)), uy.labels, col=col[4])
           panel.text(y=uy, x=rep(min(y)-.5, length(uy)), uy.labels, col=col[4])

           ## means on m axis
           panel.abline(a=0, b=1, col=col[5], lty=2)
           panel.text(x=max(y)+c(1.1,.6), y=max(y)+c(.5,1.5),
                      c("d = h-v = 0","(m-axis)"),
                      srt=45, adj=0,
                      col=col[5])
           ## ticks on m axis
           panel.segments(x1=uy-.05*diff(range(y)),
                          y1=uy+.05*diff(range(y)),
                          x2=uy+.05*diff(range(y)),
                          y2=uy-.05*diff(range(y)),
                          col=col[5])
           ## dotted line continuation of m axis
           panel.segments(x1=43.5, y1=43.5, x2=45.5, y2=45.5,
                          lty=3, col=col[5])

           ## (0,0) point and axes
           panel.segments(x1=44.5, y1=44.25, x2=44.5, y2=47, col=col[1])
           panel.segments(x1=44.25, y1=44.5, x2=47, y2=44.5, col=col[1])
           panel.text(x=c(44,44.5),y=c(44.5,44), c("0","0"), col=col[1], cex=.8)

           ## d axis
           panel.segments(x1=40.5, y1=48.5, x2=46.5, y2=42.5, col=col[5], lty=2)
           panel.text(x=c(42.1,41.9), y=c(47.5,46.5),
                      c("(d-axis)","m = (h+v)/2 = 0"), srt=-45, col=col[5])
           if (exit==1) return()
           ## export.eps(hh("mcomp/figure/mmc1-a.eps"))


           ## index for means on m axis
           panel.text(x=uy-.3,
                      y=uy+.3,
                      uy.labels, col=col[5])

           ## ticks in m coordinates
           panel.abline(a=8, b=1,
                        col=col[2])
           panel.segments(x1=uy-8/2,
                          y1=uy+8/2,
                          x2=uy-8/2+.05*diff(range(y)),
                          y2=uy+8/2-.05*diff(range(y)),
                          col=col[2])
           ## m labels, ybar
           panel.text(x=uy-8/2-.4,
                      y=uy+8/2+.4,
                      round(uy,1),
                      col=col[2])
           panel.text(x=uy[1]-8/2-.4+.8,
                      y=uy[1]+8/2+.4+.8,
                      "m",
                      col=col[2])
           panel.text(x=max(y)-9.3, y=max(y)+1,
                      "d = h-v = constant",
                      srt=45, adj=0,
                      col=col[2])
           ## perpendiculars from A to m-axis
           panel.segments(uy[1],uy[1], uy[1]-8/2, uy[1]+8/2, lty=2, lwd=3,
                          col=col[3]) ## to m-axis

           if (exit==2) return()
           ## export.eps(hh("mcomp/figure/mmc1-b0.eps"))


           ## differences on d axis
           panel.abline(a=2*min(y)-.2*diff(range(y))-.8, b=-1,
                        col=col[2])
           panel.text(x=min(y)-2.4, y=min(y)-1.5, srt=-45,
                      "m = (h+v)/2 = constant",
                      col=col[2])
           panel.segments(x1=uy-(uy-min(uy))/2-.1*diff(range(y))-.4,
                          y1=min(uy)-(uy-min(uy))/2-.1*diff(range(y))-.4,
                          x2=uy-(uy-min(uy))/2-.05*diff(range(y))-.4,
                          y2=min(uy)-(uy-min(uy))/2-.05*diff(range(y))-.4,
                          col=col[2])
           panel.text(x=uy-(uy-min(uy))/2+.0*diff(range(y))-.6,
                      y=min(uy)-(uy-min(uy))/2+.0*diff(range(y))-.4,
                      c(paste(uy.labels[1], uy[1]-min(uy), sep=": "),
                        uy.labels[2:4]),
                      adj=0,
                      col=col[2])
           ## ticks
           panel.segments(x1= (0:6)/2+min(uy)-.10*diff(range(y))-.4,
                          y1=-(0:6)/2+min(uy)-.10*diff(range(y))-.4,
                          x2= (0:6)/2+min(uy)-.15*diff(range(y))-.4,
                          y2=-(0:6)/2+min(uy)-.15*diff(range(y))-.4,
                          col=col[2])
           ## tick labels
           panel.text(x= (0:7)/2+min(uy)-.18*diff(range(y))-.4,
                      y=-(0:7)/2+min(uy)-.18*diff(range(y))-.4,
                      c(0:6, "d"),
                      col=col[2])
           ## perpendicular from A to d-axis
           panel.segments(uy[1],uy[4], uy[1]-3.8, uy[4]-3.8, lty=2, lwd=3,
                          col=col[3]) ## to d-axis
           if (exit==3) return()
           ## export.eps(hh("mcomp/figure/mmc1-b.eps"))

           ## CI for ybar2-ybar4
           panel.segments(x1=uy[2]-ms.5*sqrt(1/n[2]+1/n[4])*crit.point/2,
                          y1=uy[4]+ms.5*sqrt(1/n[2]+1/n[4])*crit.point/2,
                          x2=uy[2]+ms.5*sqrt(1/n[2]+1/n[4])*crit.point/2,
                          y2=uy[4]-ms.5*sqrt(1/n[2]+1/n[4])*crit.point/2,
                          lwd=4,
                          col=col[6])
           ## export.eps(hh("mcomp/figure/mmc1.eps"))
         }
         )
}
