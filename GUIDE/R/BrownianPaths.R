BrownianPaths <-
function(){
  
  initialprice <- 0 # lets keep this fixed
  time <- 1
  steps <- 250
  dt <- time/steps
  #   paths <- 6
  sigma <- 1
  sigmasqrtdt <- sigma*sqrt(dt)
  
  
  my.draw <- function(panel) {
    paths <- panel$paths
    
    #   ds = S + S(mudt +sigmasqrtdt*epsilon)  
    Brownianprices <- function(initialprice,steps,paths,sigmasqrtdt){
      prices <- matrix(data=NA, nrow=steps+1,ncol=paths) 
      prices[1,] <- initialprice
      for (i in 2:(steps+1)){
        prices[i,] <- prices[i-1,]+ sigmasqrtdt*rnorm(paths)
      }
      return(prices)
    }
    
    prices <-Brownianprices(initialprice,steps,paths,sigmasqrtdt)
    
    x.axis <- seq(0,1,length=steps+1)  
    
    if (length(dev.list()) == 0) 
      dev.new()
    my.title <- paste(paths, " Brownian Motions")
    matplot(prices,main= my.title,xlab="time", ylab="price",type='l',lwd=2)
    panel
  }
  
  my.redraw <- function(panel) #not needed bcos we are not using tkr plot
  {
    rp.tkrreplot(panel, my.tkrplot)
    panel                                                                       
  }
  
  my.panel <- rp.control(title = "Brownian Motion", paths=1, size=c(500,400))
  
  
  #rp.doublebutton(panel = my.panel, variable= sigma, step = 0.05, range = c(0.00, 0.50),
  #                   title = "volatility",  action = my.draw)
  #rp.doublebutton(panel = my.panel, variable= mu, step = 0.04, range = c(0.00, 0.20),
  #                   title = "drift", action = my.draw)
  #   rp.tkrplot(my.panel, my.tkrplot, my.draw) # doesnt appear very good
  rp.doublebutton(panel = my.panel, variable= paths, step = 1, range = c(1, 10),
                  title = "Paths", action = my.draw)
  rp.do(my.panel, my.draw)
}
