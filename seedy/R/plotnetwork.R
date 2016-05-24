plotnetwork <-
function(probmat, labels, arrlen=0.15, scale=1, ...) {
  n <- nrow(probmat)
  if (missing(labels)) {
    labels <- 1:n
  }
  plot(NULL, xlim=c(-1.2,1.2), ylim=c(-1.2,1.2), xaxt="n", yaxt="n", xlab="", ylab="", ...)
  points(cos(2*(1:n)*pi/n), sin(2*pi*(1:n)/n), pch=16, col="red")
  text(1.2*cos(2*(1:n)*pi/n), 1.2*sin(2*pi*(1:n)/n), labels)
  for (i in 1:n) { # TO i
    for (j in (1:n)[-i]) { # FROM j
      fct <- probmat[j,i]^scale
      
      arrows(cos(2*j*pi/n), sin(2*pi*j/n), cos(2*i*pi/n), sin(2*pi*i/n), col=rgb(0,0.5,0.75*fct,fct), 
             lwd=0.5+1.2*fct, angle=10, length=arrlen)
    }
  }
}
