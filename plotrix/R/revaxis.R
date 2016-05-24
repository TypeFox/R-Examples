# Function revaxis --- rev for ``reverse'' (to reverse the direction of
# axes).
#
# Written by T. Rolf Turner, University of New Brunswick (now
# with the Starpath Project, Universtiy of Auckland).
#
# Some bugs fixed 7/2/02, with the assistance of Herberto Ghezzo
# of McGill University.
#
revaxis<-function(x,y,xrev=FALSE,yrev=TRUE,xside=if(yrev) 3 else 1,
 yside=if(xrev) 4 else 2,xlab=NULL,ylab=NULL,bty=NULL,...) {
 
 xname <- if(is.null(xlab)) deparse(substitute(x)) else xlab
 yname <- if(is.null(ylab)) deparse(substitute(y)) else ylab
 xlab <- if(yrev) "" else xname
 ylab <- if(xrev) "" else yname
 y1 <- if(yrev) -y else y
 x1 <- if(xrev) -x else x
 old.mar <- par()$mar
 on.exit(par(mar = old.mar))
 par(mar = old.mar[c(xside, yside, 4 - xside, 6 - yside)])
 plot(x1, y1, axes = FALSE, xlab = xlab, ylab = ylab, ...)
 if(xrev) {
  axis(xside, at = pretty(-x), labels = rev(pretty(x)))
  mtext(side = yside, line = 2, text = yname)
 }
 else axis(xside)
 if(yrev) {
  axis(yside,at=pretty(-y),labels=rev(pretty(y)),srt=90)
  mtext(side=xside,line=3,text=xname)
 }
 else axis(yside)
 if(!is.null(bty)) box(bty = bty)
 invisible()
}
