if (getRversion() >= "2.15.1") utils::globalVariables(c( 'opttype',
                                                         'greek'))
stockTimeGreeks <-
function(){
  strike = 100 # lets keep this fixed
  colors=c("greenyellow" ,  "deepskyblue" ,"darksalmon", "lightgreen", "gold")
  
  
  my.draw <- function(panel) {
    
    rf <- panel$rf
    sigma <- panel$sigma
    opttype <- panel$opttype
    greek <- panel$greek
    
    BSMGreek <- function(stockprice, time){
      d1 <- (log(stockprice/strike) + (rf + (sigma^2/2)*time))/(sigma * sqrt(time))
      d2 <- d1 - (sigma * sqrt(time))
      
      div <- 0
      if (opttype=="Call"){
        if (greek=="Delta"){
          greekval <- exp( -div*time ) * pnorm(d1)

        }
        else if (greek=="Gamma"){
         greekval <- exp(-div*time) * ( dnorm(d1)/( stockprice*sigma*sqrt(time)) )

         }
        else if (greek=="Vega"){
          greekval  <- stockprice * exp(-div*time) * dnorm(d1) * sqrt(time)

          }
        else if (greek=="Theta"){
          greekval <- exp(-div*time) * (stockprice*dnorm(d1)*sigma)/(2*sqrt(time)) -
          rf * strike * exp(-rf*time) * pnorm(d2)

          }
        else {
          greekval  <- strike * time * exp(-rf*time) * pnorm(d2)

          }
        greek.name <- greek
        
      }
#       in case of put
      else{
        if (greek=="Delta"){ 
       greekval <- -exp(-div*time ) * pnorm(-d1) 

       }
        else if (greek=="Gamma"){
          greekval <- exp(-div*time) * ( dnorm(d1)/( stockprice*sigma*sqrt(time)) )

          }
        else if (greek=="Vega"){
          greekval  <- stockprice * exp(-div*time) * dnorm(d1) * sqrt(time)

          }
        else if (greek=="Theta"){
          greekval <- exp(-div*time) * (stockprice*dnorm(d1)*sigma)/(2*sqrt(time)) +
          rf * strike * exp(-rf*time) * pnorm(-d2)

          }
        else {
          greekval  <- strike * time * exp(-rf*time) * pnorm(-d2)

          }
        greek.name <- greek
      }
      
      return(greekval)
    }
    
    
    
    stockprice <- seq(strike*0.3, strike*1.7, length=31) # lets keep this relative to stockprice
    time <- seq(0, 1, length=31) # lets keep this fixed
    
    greekval <- outer(stockprice, time, FUN = BSMGreek)

#     colors <- c("cyan", "steelblue", "green","greenyellow" , "lightgreen","deepskyblue" ,"darksalmon","gold", "skyblue", "orange", "violet")
#     colorNumber <- 1:length(colors)
#     thiscolor <- colors[sample(colorNumber)[1]]
    
#     choose a different color for each greek
    whichcolor <- function(greek){
      if (greek=="Delta"){
      thiscolor <- colors[1]
      }
      else if (greek=="Gamma"){
        thiscolor <- colors[2]
      }
      else if (greek=="Vega"){
        thiscolor <- colors[3]
      }
      else if (greek=="Theta"){
        thiscolor <- colors[4]
      }
      else {
        thiscolor <- colors[5]
      }
      return(thiscolor)
    }
    
    
    if (length(dev.list()) == 0) 
      dev.new()
    my.title <- paste(opttype, greek, "(Rf= ", rf, ", sigma= ", sigma,")")
    persp(stockprice, time, greekval, zlab= greek, ticktype = "detailed", main= my.title,theta= -40,phi= 10,col=whichcolor(greek))
    panel
  }
  
  my.redraw <- function(panel) #not needed bcos we are not using tkr plot
  {
    rp.tkrreplot(panel, my.tkrplot)
    panel                                                                       
  }
  
  my.panel <- rp.control(title = "stockprice-time-Delta", rf = 0.04, sigma = 0.20,size=c(500,400))
  
  rp.radiogroup(panel = my.panel, variable= opttype,
                vals = c("Call", "Put"), 
                action = my.draw, title = "Type of Option")
  rp.radiogroup(panel = my.panel, variable= greek,
                vals = c("Delta", "Gamma","Vega","Theta","Rho"), 
                action = my.draw, title = "Greek")
  rp.doublebutton(panel = my.panel, variable= sigma, step = 0.05, range = c(0.00, 0.50),
                  title = "sigma",  action = my.draw)
  rp.doublebutton(panel = my.panel, variable= rf, step = 0.01, range = c(0.00, 0.10),
                  title = "risk free", action = my.draw)
  #rp.tkrplot(my.panel, my.tkrplot, my.draw) # doesnt appear very good
  rp.do(my.panel, my.draw)
}
