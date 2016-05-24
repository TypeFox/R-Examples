rich.colors <- function(n,
                        palette="temperature",
                        alpha=1,
                        rgb=FALSE,
                        plot=FALSE)
{
  if(n <= 0)
    return(character(0))

  palette <- match.arg(palette, c("temperature","blues"))
  x <- seq(0, 1, length=n)

  if(palette == "temperature")
  {
    r <- 1 / (1+exp(20-35*x))
    g <- pmin(pmax(0,-0.8+6*x-5*x^2), 1)
    b <- dnorm(x,0.25,0.15) / max(dnorm(x,0.25,0.15))
  }
  else
  {
    r <-        0.6*x + 0.4*x^2
    g <-        1.5*x - 0.5*x^2
    b <- 0.36 + 2.4*x - 2.0*x^2
    b[x>0.4] <- 1
  }

  rgb.m <- matrix(c(r,g,b), ncol=3,
                  dimnames=list(NULL,c("red","green","blue")))
  col <- mapply(rgb, r, g, b, alpha)

  if(rgb) 
    attr(col, "rgb") <- cbind(rgb.m, alpha)
  
  if(plot)
  {
    opar <- par("fig", "plt")
    par(fig=c(0,1,0,0.7), plt=c(0.15,0.9,0.2,0.95))
    plot(NA, xlim=c(-0.01,1.01), ylim=c(-0.01,1.01), xlab="Spectrum", ylab="",
         xaxs="i", yaxs="i", axes=FALSE)
    title(ylab="Value", mgp=c(3.5,0,0))
    matlines(x, rgb.m, col=colnames(rgb.m), lty=1, lwd=3)
    matpoints(x, rgb.m, col=colnames(rgb.m), pch=16)
    axis(1, at=0:1)
    axis(2, at=0:1, las=1)
    par(fig=c(0,1,0.75,0.9), plt=c(0.08,0.97,0,1), new=TRUE)
    midpoints <- barplot(rep(1,n), col=col, border=FALSE, space=FALSE,
                         axes=FALSE)
    axis(1, at=midpoints, labels=1:n, lty=0, cex.axis=0.6)
    par(opar)
  }

  return(col)
}
