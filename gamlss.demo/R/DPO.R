demo.DPO <- function()
{
  if (interactive()) 
  {
    mu <- sigma <- NULL
    density.draw <- function(panel) 
    {
      
      op<-par(mfrow = c(2, 1))
      x <- seq(0, 50, 1)
      plot(x, dDPO(x, mu=panel$mu, sigma=panel$sigma), type="h", ylab="f(x)", main="DPO probability function")
      points(x, dDPO(x, mu=panel$mu, sigma=panel$sigma), col="blue")
      plot(x, pDPO(x, mu=panel$mu, sigma=panel$sigma), type="s", ylab="F(x)", main="NBI cumulative distribution function")
      par(op) 
      panel
    }
    ppoispanel <- rp.control('DPO',  mu= 5, sigma=1,  ymax=1 )
    rp.slider(ppoispanel,  variable=mu, from=0.1, to=30, resolution=0.01,  action = density.draw, title="mu",  showvalue = TRUE) 
    rp.slider(ppoispanel,  variable=sigma, from=0.1, to=5, resolution=0.01,  action = density.draw, title="sigma",  showvalue = TRUE)  
  }
}
#------------------------------------------------------------------------------