if (getRversion() >= "2.15.1") utils::globalVariables(c("nu","delta","lambda"))


JDPaths <-
function(){
  #   http://www.thetaris.com/wiki/Jump_Diffusion
  
  
  initialprice <- 100 # lets keep this fixed
  time <- 1
  steps <- 250
  dt <- time/steps
  #   paths <- 6
  
  
  
  my.draw <- function(panel) {
    mu <- panel$mu
    sigma <- panel$sigma
    paths <- panel$paths
    
    nu <- panel$nu # mean of jumps
    delta <- panel$delta # std of jumps
    lambda <- panel$lambda # jump intensity or # jumps/year
    
    
    
    #   ds = S + S(mudt +sigmasqrtdt*epsilon)  
    JDStockprices <- function(initialprice,steps,paths,sigma,mu,nu,delta,lambda){
      
      
      prices <- matrix(data=NA, nrow=steps+1,ncol=paths) 
      prices[1,] <- initialprice
      for (i in 2:(steps+1)){
        P = rpois(n=paths,lambda=lambda*dt)
        U = exp(P*nu + sqrt(P)*delta*rnorm(1)) # the poisson process
        prices[i,] = prices[i-1,] * exp((mu-lambda*(exp(nu+0.5*delta^2)-1)-0.5*sigma^2)*dt + sigma*sqrt(dt)*rnorm(1))*U
        
      }
      return(prices)
    }
    
    prices <-JDStockprices(initialprice,steps,paths,sigma,mu,nu,delta,lambda)
    
    x.axis <- seq(0,1,length=steps+1)  
    
    if (length(dev.list()) == 0) 
      dev.new()
    my.title <- paste(paths, " Jump Diffusion Motions", "(mu=", mu, ", sigma=", sigma,")")
    matplot(prices,main= my.title,xlab="time", ylab="price",type='l',lwd=2)
    panel
  }
  
  my.redraw <- function(panel) #not needed bcos we are not using tkr plot
  {
    rp.tkrreplot(panel, my.tkrplot)
    panel                                                                       
  }
  
  my.panel <- rp.control(title = "Jump Diffusion", mu = 0.1, sigma = 0.25, nu=-0.1, delta=0.05,lambda=1,paths=1, size=c(500,400))
  
  
  rp.doublebutton(panel = my.panel, variable= mu, step = 0.04, range = c(0.00, 0.20),
                  title = "Drift", action = my.draw)
  rp.doublebutton(panel = my.panel, variable= sigma, step = 0.05, range = c(0.00, 0.50),
                  title = "Volatility",  action = my.draw)
  #   rp.tkrplot(my.panel, my.tkrplot, my.draw) # doesnt appear very good
  rp.doublebutton(panel = my.panel, variable= nu, step = 0.1, range = c(-0.2, 0.2),
                  title = "Mean of Jumps", action = my.draw)
  rp.doublebutton(panel = my.panel, variable= delta, step = 0.1, range = c(0, 0.3),
                  title = "Std Dev of jumps", action = my.draw)
  rp.doublebutton(panel = my.panel, variable= lambda, step = 2, range = c(1, 10),
                  title = "Jump Intensity", action = my.draw)
  rp.doublebutton(panel = my.panel, variable= paths, step = 1, range = c(1, 10),
                  title = "Paths", action = my.draw)
  rp.do(my.panel, my.draw)
}
