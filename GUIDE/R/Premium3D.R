if (getRversion() >= "2.15.1") utils::globalVariables(c("opttype","plottype"))

Premium3D <-
function(){
  
  # stockprice = 100 # lets keep this fixed
  
  
  
  my.draw <- function(panel) {
    
    rf <- panel$rf
    sigma <- panel$sigma
    opttype <- panel$opttype
    plottype <- panel$plottype
    
    BSMPrice <- function(x.axis.variable, time){
      d1 <- (log(stockprice/strike) + (rf + (sigma^2/2)*time))/(sigma * sqrt(time))
      d2 <- d1 - (sigma * sqrt(time))
      Nd1 <- pnorm(q=d1)
      Nd2 <- pnorm(q=d2)
      if (opttype=="Call"){
        price <- (stockprice * Nd1) - (strike * exp(-rf*time) * Nd2)
      }
      else{
        price <- (strike * exp(-rf*time) * (1-Nd2)) - (stockprice * (1-Nd1)) 
      }
      
      return(price)
    }
    
    if (panel$plottype=="Strike-Time"){
      stockprice <- 100
      strike <- seq(stockprice/2, stockprice*2, length=31) # lets keep this relative to stockprice
      x.axis.variable <- strike
      x.axis.variable.name <- "strike"
      }
    else{
      strike <- 100
      stockprice <- seq(strike*.3, strike*1.7, length=31)
      x.axis.variable <- stockprice
      x.axis.variable.name <- "stock price"
      }
    time <- seq(0, 1, length=31) # lets keep this fixed
    
    premium <- outer(x.axis.variable, time, FUN = BSMPrice)
    
    if (length(dev.list()) == 0) 
      dev.new()
    colors <- c("cyan", "steelblue", "green","greenyellow" , "lightgreen","deepskyblue" ,"darksalmon","gold", "skyblue", "orange", "violet")
    my.title <- paste(opttype, " Premium (Rf= ", rf, ", sigma= ", sigma,")")
    persp(x.axis.variable, time, premium, xlab = x.axis.variable.name, ticktype="detailed",main= my.title,theta= -40,phi= 10,col="cyan")
    panel
  }
  
  my.redraw <- function(panel) #not needed bcos we are not using tkr plot
  {
    rp.tkrreplot(panel, my.tkrplot)
    panel                                                                       
  }
  
  my.panel <- rp.control(title = "Premium 3D Plots", rf = 0.04, sigma = 0.20,size=c(500,400))
  
  rp.radiogroup(panel = my.panel, variable= opttype,
                vals = c("Call", "Put"), 
                action = my.draw, title = "Type of Option")
  rp.radiogroup(panel = my.panel, variable= plottype,
                vals = c("Stockprice-Time", "Strike-Time"), 
                action = my.draw, title = "X-Y axis: ")
  rp.doublebutton(panel = my.panel, variable= sigma, step = 0.05, range = c(0.00, 0.50),
                  title = "sigma",  action = my.draw)
  rp.doublebutton(panel = my.panel, variable= rf, step = 0.01, range = c(0.00, 0.10),
                  title = "risk free", action = my.draw)
  #rp.tkrplot(my.panel, my.tkrplot, my.draw) # doesnt appear very good
  rp.do(my.panel, my.draw)
}
