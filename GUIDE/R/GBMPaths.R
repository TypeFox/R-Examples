GBMPaths <-
function(){
  
  initialprice <- 100 # lets keep this fixed
  time <- 1
  steps <- 250
  dt <- time/steps
  #   paths <- 6
  
  
  
  my.draw <- function(panel) {
    mu <- panel$mu
    sigma <- panel$sigma
    paths <- panel$paths
    
    
    #   ds = S + S(mudt +sigmasqrtdt*epsilon)  
    GBMStockprices <- function(initialprice,steps,paths,sigma,mu){
      mudt <- mu*dt
      sigmasqrdt <- sigma^2*dt
      sigmasqrtdt <- sigma*sqrt(dt)
      
      
      prices <- matrix(data=NA, nrow=steps+1,ncol=paths) 
      prices[1,] <- initialprice
      for (i in 2:(steps+1)){
        #   prices[i,] <- prices[i-1,]+prices[i-1,]* (mudt + sigmasqrtdt*rnorm(paths))
        prices[i,] <- prices[i-1,]* exp(mudt - sigmasqrdt/2 + sigmasqrtdt*rnorm(paths))
      }
      return(prices)
    }
    
    prices <-GBMStockprices(initialprice,steps,paths,sigma,mu)
    
    x.axis <- seq(0,1,length=steps+1)  
    
    if (length(dev.list()) == 0) 
      dev.new()
    my.title <- paste(paths, " Geometric Brownian Motions", "(mu=", mu, ", sigma=", sigma,")")
    matplot(prices,main= my.title,xlab="time", ylab="price",type='l',lwd=2)
    panel
  }
  
  my.redraw <- function(panel) #not needed bcos we are not using tkr plot
  {
    rp.tkrreplot(panel, my.tkrplot)
    panel                                                                       
  }
  
  my.panel <- rp.control(title = "Geometric Brownian Motion", mu = 0.1, sigma = 0.25, paths=1, size=c(500,400))
  
  
  rp.doublebutton(panel = my.panel, variable= mu, step = 0.04, range = c(0.00, 0.20),
                  title = "Drift", action = my.draw)
  rp.doublebutton(panel = my.panel, variable= sigma, step = 0.05, range = c(0.00, 0.50),
                  title = "Volatility",  action = my.draw)
  #   rp.tkrplot(my.panel, my.tkrplot, my.draw) # doesnt appear very good
  rp.doublebutton(panel = my.panel, variable= paths, step = 1, range = c(1, 10),
                  title = "Paths", action = my.draw)
  rp.do(my.panel, my.draw)
}
