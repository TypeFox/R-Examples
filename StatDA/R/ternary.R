ternary <- function(x, nam=NULL, grid=FALSE, ...)
{
# Ternary plot
#
# x ... matrix with 3 columns
# nam ... names of the variables
# grid ... TRUE if grid should be plotted
# "..." ... further graphical parameters, see par

val=0.6 # value for grey 

if (is.null(nam)) { nam <- dimnames(x)[[2]]}
s <- rowSums(x)
if (any(s <= 0))
        stop("each row of the input `object' must have a positive sum")
dat <- x/s

xp <- dat[,2] + dat[,3]/2
yp <- dat[,3] * sqrt(3)/2

par(pty="s")
plot(xp,yp,xlim=c(0,1),ylim=c(0,0.9), 
   frame.plot=FALSE, xaxt="n", yaxt="n", xlab="", ylab="", ...)

segments(0,0,1,0)
segments(0,0,1/2,sqrt(3)/2)
segments(1/2,sqrt(3)/2,1,0)

mtext(nam[1],side=1, line=-1, at=-0.05,cex=1.2)
mtext(nam[2],side=1, line=-1, at=1.05,cex=1.2)
text(0.5, 0.9, nam[3],cex=1.2)

if(grid==TRUE)
{
segments(0.2,0,0.1,sqrt(0.03), col=grey(val), lty="dashed")
segments(0.4,0,0.2,sqrt(0.12), col=grey(val), lty="dashed")
segments(0.6,0,0.3,sqrt(0.27), col=grey(val), lty="dashed")
segments(0.8,0,0.4,sqrt(0.48), col=grey(val), lty="dashed")
segments(0.2,0,0.6,sqrt(0.48), col=grey(val), lty="dashed")
segments(0.4,0,0.7,sqrt(0.27), col=grey(val), lty="dashed")
segments(0.6,0,0.8,sqrt(0.12), col=grey(val), lty="dashed")
segments(0.8,0,0.9,sqrt(0.03), col=grey(val), lty="dashed")
segments(0.1,sqrt(0.03),0.9,sqrt(0.03), col=grey(val), lty="dashed")
segments(0.2,sqrt(0.12),0.8,sqrt(0.12), col=grey(val), lty="dashed")
segments(0.3,sqrt(0.27),0.7,sqrt(0.27), col=grey(val), lty="dashed")
segments(0.4,sqrt(0.48),0.6,sqrt(0.48), col=grey(val), lty="dashed")

text(0.5,0.66,"0.8", col=grey(val), cex = 0.6)
text(0.5,0.49,"0.6", col=grey(val), cex = 0.6)
text(0.5,0.32,"0.4", col=grey(val), cex = 0.6)
text(0.5,0.14,"0.2", col=grey(val), cex = 0.6)
text(0.95,0.21,"0.8", col=grey(val), cex = 0.6, srt = 60)
text(0.86,0.35,"0.6", col=grey(val), cex = 0.6, srt = 60)
text(0.75,0.54,"0.4", col=grey(val), cex = 0.6, srt = 60)
text(0.64,0.72,"0.2", col=grey(val), cex = 0.6, srt = 60)
text(0.05,0.21,"0.8", col=grey(val), cex = 0.6,srt = 300)
text(0.14,0.35,"0.6", col=grey(val), cex = 0.6,srt = 300)
text(0.25,0.54,"0.4", col=grey(val), cex = 0.6,srt = 300)
text(0.36,0.72,"0.2", col=grey(val), cex = 0.6,srt = 300)
}

}

