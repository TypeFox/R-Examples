require(lattice)
require(latticeExtra)


bWA <- read.table(text="
color   x       y
red	1	6
blue	2	2
green	3	10
orange	4	5
brown	5	18
purple	6	12
", header=TRUE, stringsAsFactors=FALSE)


betaWeightedAverageLattice <-
  function(x, y,
           xbar=mean(x), ybar=mean(y),
           beta1.hat, beta0.hat,
           col=seq(length=length(x))+1,
           color.summary=1,
           summary=TRUE,
           summary.text=TRUE,
           cex=1.5, pch=19,
           ...) {
    if (missing(beta1.hat) || missing(beta0.hat)) {
      beta.hat <- coef(lm(y ~ x))
      beta1.hat <- beta.hat[2]
      beta0.hat <- beta.hat[1]
    }
    xyplot(y ~ x, col=col,
           xbar=xbar,
           ybar=ybar,
           beta1.hat=beta1.hat,
           beta0.hat=beta0.hat,
           col=col,
           color.summary=color.summary,
           summary=summary,
           summary.text=summary.text,
           cex=cex, pch=pch,
           ...,
           panel=function(x, y, ...) {
             panel.points(x=xbar, y=ybar, pch=8, cex=2, lwd=2, col=color.summary)
             panel.axis("bottom", at=xbar, labels=expression(bar(x)), rot=0, outside=TRUE)
             panel.axis("bottom", at=xbar, labels="")
             panel.axis("top", at=xbar, labels="")
             panel.axis("top", at=xbar, labels="", outside=TRUE)
             panel.axis("left", at=ybar, labels=expression(bar(y)), outside=TRUE)
             panel.axis("left", at=ybar, labels="")
             panel.axis("right", at=ybar, labels="")
             panel.axis("right", at=ybar, labels="", outside=TRUE)
             if (summary)
               panel.abline(a=beta0.hat, b=beta1.hat, lwd=2, col=color.summary)
             if (summary.text) {
               ## panel.text(x=xbar-1, y=ybar+2,
               ##            expression((bar(x) * ',' * bar(y))),
               ##            col=color.summary, cex=2)
               panel.text(x=6.6, y=14.2,
                          expression(y = hat(beta)[0] + hat(beta)[1] * x),
                          col=color.summary, cex=1.4,
                          xpd=NA, adj=0)
               panel.text(x=8, y=14.2,
                          paste("  =", round(beta0.hat, 3), "+", round(beta1.hat, 3), "x"),
                          col=color.summary, cex=1.4,
                          xpd=NA, adj=0)
             }
             panel.segments(x, y,    x,    ybar, col=col, lty=2, lwd=1)
             panel.segments(x, ybar, xbar, ybar, col=col, lty=2, lwd=1)
             panel.segments(x, y,    xbar, ybar, col=col, lty=1, lwd=3)
           },
           par.settings=list(clip=list(panel=FALSE))
           )
  }

## betaWeightedAverageLattice(bWA$x, bWA$y, col=bWA$color, xlim=c(-.2, 12.2))


dx <- bWA$x-mean(bWA$x)
dy <- bWA$y-mean(bWA$y)
betas <- dy / dx
wts <- dx^2 / sum(dx^2)

Six <- do.call(c,
lapply(1:nrow(bWA),
       function(i)
       betaWeightedAverageLattice(bWA$x[i],  bWA$y[i],
                                  col=bWA$color[i],
                                  xbar=mean(bWA$x), ybar=mean(bWA$y),
                                  summary=FALSE,
                                  xlim=c(0,7), ylim=c(0,20), xlab=NULL, ylab=NULL,
                                  summary.text=FALSE
                                  )
       )
)
## Six

print(position=c(0, .6, 1, 1), more=TRUE,
update(Six,
       strip=strip.custom(factor.levels=paste(
                            "slope = ", round(betas, 3), "\n",
                            "weight = ", round(wts, 3),
                            sep="")
         ),
       par.strip.text=list(lines=2.5, cex=.7),
       scales=list(alternating=0),
       between=list(x=1, y=1),
       as.table=TRUE,
       layout=c(6,1)
       )
      )


beta.hat <- coef(lm(y ~ x, data=bWA))
print(position=c(.2, 0, .8, .6), more=FALSE,
update(betaWeightedAverageLattice(bWA$x, bWA$y, col=bWA$color, xlim=c(-2, 9), summary.text=FALSE, xlab=NULL, ylab=NULL),
       strip=strip.custom(factor.levels=paste("slope =", as.character(round(beta.hat[2], 3)))))
)


sum(wts)
sum(betas*wts)
beta.hat
