density.draw <- function(panel) 
{
  print(panel$h[1])
  print(panel$h[2])
  plot(density(as.numeric(panel$x), bw=panel$h[1]*panel$h[2]))
  panel
}

panel <- rp.control(x=rnorm(50))
rp.slider(panel, h, c(0.1, 0.01), c(5, 1), title=c("Log", "Linear"), initval=c(2, 1), log=c(TRUE, FALSE), showvalue=TRUE, action=density.draw)
