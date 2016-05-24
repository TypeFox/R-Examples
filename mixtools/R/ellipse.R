ellipse <- function(mu, sigma, alpha=.05, npoints=250,
                    newplot=FALSE, draw=TRUE, ...) {
  es <- eigen(sigma)
  e1 <- es$vec%*%diag(sqrt(es$val))
  r1 <- sqrt(qchisq(1-alpha,2))
  theta <- seq(0,2*pi,len=npoints)
  v1 <- cbind(r1*cos(theta),r1*sin(theta))
  pts=t(mu-(e1%*%t(v1)))
  if (newplot && draw) {
    plot(pts, ...)
  } else if (!newplot && draw) {
    lines(pts, ...)
  }
  invisible(pts)
}

  
