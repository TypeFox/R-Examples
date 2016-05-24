if (getRversion() >= "2.15.1") utils::globalVariables(c('tkrp',
                                                        'longput',
                                                        'longcall',
                                                        'profit'))

strangle <-
function(){
  
  
  my.draw <- function(panel) 
  {
    
    with(panel, {
      ymax = k1 - minprice
      ymin =  -(c1 + p1)
      
      plot(c(minprice,k1,k2,maxprice),c(0,0,0,0),
           xlab = "Stock price",
           ylab = "Profit",type='l',ylim=c(ymin,ymax))
      drawprofits(panel)
      title(paste("Strangle with call and put"))
    })
    panel
  }
  
  
  drawprofits<-function(object)  {
    with(object,{
      
      if (longput){
        S = c(minprice,k1,k2,maxprice)
        profits1<- pmax(k1-S,0)-p1
        lines(S,profits1,type='l',col="blue", lwd=2)
      }
      
      if (longcall){
        S = c(minprice,k1,k2,maxprice)
        profits2<- pmax(S-k2,0)-c1
        lines(S,profits2,type='l',col="red", lwd=2)  
      }
      
      if (profit){
        if((!longput) | (!longcall)){
          rp.messagebox("Check the both Long Put and Long Call checkboxes to see profit graph.", title = "Insufficient Information")
        }
        else{
          S = c(minprice,k1,k2,maxprice)
          profits<- profits1+profits2
          lines(S,profits,type='l',col="black",lwd=2)  
        }
      }
    })
    object
    
  }
  
  
  my.redraw <- function(panel)
  {
    rp.tkrreplot(panel, tkrp)
    panel                                                                       
  }
  
  
  my.panel <- rp.control(title="Strangle with call and put", k1= 80, k2 = 120, maxprice=200, minprice=0,p1=10,c1=12,longput=FALSE,longcall=FALSE,profit=FALSE)
  rp.tkrplot(panel=my.panel , name=tkrp, plotfun=my.draw, pos="right",hscale=2,vscale=2)
  rp.checkbox(my.panel,variable=longput, action = my.redraw, title = "Long Put")
  rp.checkbox(my.panel,variable=longcall, action = my.redraw, title = "Long Call")
  rp.checkbox(my.panel,variable=profit, action = my.redraw, title = "Profit")
}
