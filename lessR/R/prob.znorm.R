prob.znorm <- 
function(mu=0, sigma=1, color.border="gray10",
         r=.10, g=.34, b=.94, a=.20,
         xlab="", ylab="", main="", 
         y.axis=FALSE, z=TRUE, mag=.9, ...) {

# plot normal curve with integer SD lines
    
if ( (r<0 || r>1)  ||  (g<0 || g>1)  || (b<0 || b>1)  ||  (a<0) || a>1) { 
  cat("\n"); stop(call.=FALSE, "\n","------\n",
    "Values of r, g, b and a must all be between 0 and 1, inclusive.\n\n")
}

  dots <- list(...)  # check for deprecated parameters
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      if (substr(names(dots)[i], 1, 4) == "col.") {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "options that began with the abbreviation  col  now begin with  ",
          "color \n\n")
      }
    }
  }
 
if (mu==0  && sigma==1) z=FALSE

xmin <- mu - 4*sigma
xmax <- mu + 4*sigma

cuts <- seq(xmin,xmax,sigma)

color.sig <- rgb(r, g, b, a)
#par(mar=c(3, 3, 1, 2), mgp=c(2,.6,0))

# Normal Curve
x <- seq(xmin, xmax, length=200)
y <- dnorm(x ,mean=mu ,sd=sigma)

if (sys.nframe() == 1) .graphwin(1)  # do not open up new window if called from sim.CLT
plot(x,y, type="l", lwd=2, col=color.border, axes=FALSE, xlab="", ylab="", main=main)
if (z) title(xlab=xlab, line=3.5) else title(xlab=xlab)

abline(h=0)
axis(side=1, at=cuts, cex.axis=mag)
if (z) axis(side=1, at=cuts, cex.axis=mag, line=1.5, labels=-4:4, lwd=0, lwd.ticks=0)
if (y.axis) {
  axis(side=2, cex.axis=mag)
  if (ylab == "") ylab="Normal Density"
  title(ylab=ylab)
}

segments(mu, 0, mu, dnorm(mu, mean=mu, sd=sigma), col=color.border, lty="dotted")

xsub <- x>(mu-3*sigma) & x<(mu+3*sigma)
polygon(c(mu-3*sigma,x[xsub],mu+3*sigma), c(0,y[which(xsub)],0),
        col=color.sig, border=color.border, lty="dotted")

xsub <- x>(mu-2*sigma) & x<(mu+2*sigma)
polygon(c(mu-2*sigma,x[xsub],mu+2*sigma), c(0,y[which(xsub)],0),
        col=color.sig, border=color.border, lty="dotted")

xsub <- x>(mu-sigma) & x<(mu+sigma)
polygon(c(mu-sigma,x[xsub],mu+sigma), c(0,y[which(xsub)],0),
        col=color.sig, border=color.border, lty="dotted")

}
