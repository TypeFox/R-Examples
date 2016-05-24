density.draw <- function(panel, h) 
{
  plot(density(panel$x, bw = panel$h), col="blue")
  panel
}

panel <- rp.control(x=rnorm(50))
rp.doublebutton(panel, var=h, step=1.1, title="Density estimate", action=density.draw, range=c(0.1, 5), initval=1, showvalue=TRUE, log=TRUE, foreground="Blue", background="White")
