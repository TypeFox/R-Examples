### R code from vignette source 'pBrackets_vignette.rnw'

###################################################
### code chunk number 1: pBrackets_vignette.rnw:13-15
###################################################
options(width=80, continuation="   ")
library(pBrackets)


###################################################
### code chunk number 2: plot1
###################################################

par(mar=c(1,1,1,1))
plot(0,0, type='n', xlim=c(0,20), ylim=c(0,20), axes=FALSE, xlab='', ylab='')
abline(h=seq(0,20), v=seq(0, 7), col=rgb(0.8, 0.9, 0.95))

brackets(0, 18, 7, 20, lwd=2)
text(8, 20, labels=expression(paste(bold('Braces:'), ' default')), adj=c(0,0))

brackets(0, 16, 7, 18, lwd=2, curvatur=1, type=2)
text(8, 18, labels=expression(paste(bold('Braces 2:'), ' curvatur=1, type=2')), adj=c(0,0))

brackets(0, 14, 7, 16, lwd=2, ticks=NA, curvatur=1, type=5)
text(8, 16, labels=expression(paste(bold('Parentheses:'), ' ticks=NA, curvature=1, type=5')), adj=c(0,0))

brackets(0, 12, 7, 14, lwd=2, ticks=NA, type=4, h=0.5)
text(8, 14, labels=expression(paste(bold('Square brackets:'), ' ticks=NA, type=4')), adj=c(0,0))

brackets(0, 10, 7, 12, lwd=2, ticks=NA, curvature=1, type=3)
text(8, 12, labels=expression(paste(bold('Chevrons:'), ' ticks=NA, curvature=1, type=3')), adj=c(0,0))

brackets(0, 8, 7, 10, lwd=2, ticks=NA, type=3, curvature=0.2, h=0.75)
text(8, 10, labels=expression(paste(bold('Stump brackets:'), ' ticks=NA, curvature=0.2, type=3')), adj=c(0,0))

brackets(0, 6, 7, 8, lwd=2, type=4)
text(8, 8, labels=expression(paste(bold('Square brackets with tick:'), ' type=4')), adj=c(0,0))

brackets(0, 4, 7, 6, lwd=2, ticks=c(0.25, 0.75))
text(8, 6, labels=expression(paste(bold('Double tick braces:'), ' ticks=c(0.25, 0.75)')), adj=c(0,0))

brackets(0, 2, 7, 4, lwd=2, ticks=-0.5, h=0.5)
text(8, 4, labels=expression(paste(bold('Negative tick braces:'), ' ticks=-0.5')), adj=c(0,0))

brackets(0, 0, 7, 2, lwd=2, ticks=c(-0.2, -0.4, -0.6, -0.8, 1), type=4)
text(8, 2, labels=expression(paste(bold('Multiples ticks:'), ' ticks=c(-0.2,-0.4,-0.6,-0.8, 1), type=4')), adj=c(0,0))


###################################################
### code chunk number 3: plot2
###################################################


par(mar=c(1,1,1,1))
plot(0,0, type='n', xlim=c(0,20), ylim=c(0,20), axes=F, xlab='', ylab='')
abline(h=seq(0,20), v=seq(0, 7), col=rgb(0.8, 0.9, 0.95))
#axis(side=2, at=seq(0, 20, 1))
#axis(side=1, at=seq(0, 20, 1))
brackets(0, 18, 6, 18, h=1, lwd=2)
text(8, 18, labels=expression(paste('x1=',bold('0'),',    y1=18,    x2=',bold('6'),',    y2=18,    h=1,    ...')), adj=c(0,-0.5))

brackets(6, 16, 0, 16, h=1, lwd=2)
text(8, 15, labels=expression(paste('x1=',bold('6'),',    y1=16,    x2=',bold('0'),',    y2=16,    h=1,    ...')), adj=c(0,-0.5))

brackets(0, 12, 6, 12, h=1,  lwd=2)
text(8, 12, labels=expression(paste('x1=0,    y1=12,    x2=6,    y2=12,    h=',bold('1'),',    ...')), adj=c(0,-0.5))

brackets(0, 10, 6, 10, h=-1, lwd=2)
text(8, 9, labels=expression(paste('x1=0,    y1=10,    x2=6,    y2=10,    h=',bold('-1'),',    ...')), adj=c(0,-0.5))

brackets(0, 2, 6, 8, h=sqrt(2),  lwd=2)
text(8, 6, labels=expression(paste('x1=0,    y1=4,     x2=6,    y2=8,    h=',bold('sqrt(2)'),',    ...')), adj=c(0,-0.5))

brackets(0, 0, 6, 6, h=-sqrt(2), lwd=2)
text(8, 3, labels=expression(paste('x1=0,    y1=2,    x2=6,    y2=6,    h=',bold('-sqrt(2)'),',    ...')), adj=c(0,-0.5))




###################################################
### code chunk number 4: plot3
###################################################

par(mar=c(1,1,1,1))
plot(0,0, type='n', xlim=c(0,20), ylim=c(0,20), axes=F, xlab='', ylab='')
abline(h=seq(0,20), v=seq(0, 7), col=rgb(0.8, 0.9, 0.95))
#axis(side=2, at=seq(0, 20, 1))
#axis(side=1, at=seq(0, 20, 1))
brackets(0, 19, 6, 19, h=1, lwd=2)
text(8, 19, labels=expression(paste('ticks = 0.5 (default)')), adj=c(0,-0.5))

brackets(0, 17, 6, 17, h=1, lwd=2, ticks=0.75)
text(8, 17, labels=expression(paste('ticks = 0.75')), adj=c(0,-0.75))

brackets(0, 15, 6, 15, h=1, lwd=2, ticks=0.9)
text(8, 15, labels=expression(paste('ticks = 0.9')), adj=c(0,-0.75))

brackets(0, 13, 6, 13, h=1, lwd=2, ticks=seq(0.2, 0.8, 0.2))
text(8, 13, labels=expression(paste('ticks = seq(0.2, 0.8, 0.2)')), adj=c(0,-0.5))

brackets(0, 11, 6, 11, h=1, lwd=2, ticks=-0.5)
text(8, 11, labels=expression(paste('ticks = -0.5')), adj=c(0,-0.75))

brackets(0, 9, 6, 9, h=0.5, lwd=2, ticks=-0.5)
text(8, 9, labels=expression(paste('ticks = -0.5,    h = 0.5')), adj=c(0,-0.5))

brackets(0, 7, 6, 7, h=1, lwd=2, ticks=c(-0.25, 0.5, -0.75))
text(8, 7, labels=expression(paste('ticks = c(-0.25,  0.5,  -0.75)')), adj=c(0,-0.5))

brackets(0, 5, 6, 5, h=1, lwd=2, ticks=NA)
text(8, 5, labels=expression(paste('ticks = NA')), adj=c(0,-0.75))

brackets(0, 3, 6, 3, h=0.5, lwd=2, ticks=NULL)
text(8, 3, labels=expression(paste('ticks = NULL,   h=0.5')), adj=c(0,-0.5))

brackets(0, 1, 6, 1,  h=1, lwd=2, ticks=0)
text(8, 1, labels=expression(paste('ticks = 0')), adj=c(0,-0.75))



###################################################
### code chunk number 5: plot4
###################################################

par(mar=c(1,1,1,1))
plot(0,0, type='n', xlim=c(0,20), ylim=c(3,20), axes=F, xlab='', ylab='')
abline(h=seq(0,20), v=seq(0, 7), col=rgb(0.8, 0.9, 0.95))
#axis(side=2, at=seq(0, 20, 1))
#axis(side=1, at=seq(0, 20, 1))
brackets(0, 19, 6, 19, h=1, lwd=2, curvature = 1, type=1)
text(8, 19, labels=expression(paste('curvature = ',bold('1'),',     type = 1')), adj=c(0,-0.5))

brackets(0, 17, 6, 17, h=1, lwd=2, curvature = 0.5, type=1)
text(8, 17, labels=expression(paste('curvature = ',bold('0.5'),',     type = 1')), adj=c(0,-0.75))

brackets(0, 15, 6, 15, h=1, lwd=2, curvature = 0.1, type=1)
text(8, 15, labels=expression(paste('curvature = ',bold('0.1'),',     type = 1')), adj=c(0,-0.75))

brackets(0, 13, 6, 13, h=1, lwd=2, curvature = 1, type=2)
text(8, 13, labels=expression(paste('curvature = ',bold('1'),',     type = 2')), adj=c(0,-0.5))

brackets(0, 11, 6, 11, h=1, lwd=2, curvature = 0.5, type=2)
text(8, 11, labels=expression(paste('curvature = ',bold('0.5'),',     type = 2')), adj=c(0,-0.75))

brackets(0, 9, 6, 9, h=1, lwd=2, curvature = 0.1, type=2)
text(8, 9, labels=expression(paste('curvature = ',bold('0.1'),',     type = 2')), adj=c(0,-0.5))

brackets(0, 7, 6, 7, h=1, lwd=2, curvature = 1, type=5, ticks=NA)
text(8, 7, labels=expression(paste('curvature = ',bold('1'),',     type = 5,    ticks=NA')), adj=c(0,-0.5))

brackets(0, 5, 6, 5, h=1, lwd=2, curvature = 0.5, type=5, ticks=NA)
text(8, 5, labels=expression(paste('curvature = ',bold('0.5'),',     type = 5,    ticks=NA')), adj=c(0,-0.75))

brackets(0, 3, 6, 3, h=1, lwd=2, curvature = 0.1, type=5, ticks=NA)
text(8, 3, labels=expression(paste('curvature = ',bold('0.1'),',     type = 5,    ticks=NA')), adj=c(0,-0.5))



