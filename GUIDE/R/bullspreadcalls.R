if (getRversion() >= "2.15.1") utils::globalVariables(c("tkrp","longcall","shortcall","profit"))

bullspreadcalls <-
function(){
  
  
  my.draw <- function(panel) 
  {
    
    with(panel, {
      ymax = maxprice - k1
      ymin=  k2 - maxprice
      
      plot(c(minprice,k1,k2,maxprice),c(0,0,0,0),
           xlab = "Stock price",
           ylab = "Profit",type='l',ylim=c(ymin,ymax))
      drawprofits(panel)
      title(paste("Bull Spread with calls"))
    })
    panel
  }
  
  
  drawprofits<-function(object)  {
    with(object,{
      
      if (longcall){
        S = c(minprice,k1,k2,maxprice)
        profits1<- pmax(S-k1,0)-c1
        lines(S,profits1,type='l',col="blue", lwd=2)
      }
      
      if (shortcall){
        S = c(minprice,k1,k2,maxprice)
        profits2<- -pmax(S-k2,0)+c2
        lines(S,profits2,type='l',col="red", lwd=2)  
      }
      
      if (profit){
        if((!longcall) | (!shortcall)){
          rp.messagebox("Check the both Long Call and Short Call checkboxes to see profit graph.", title = "Insufficient Information")
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
  
  
  my.panel <- rp.control(title="Bull spread with calls", k1= 50, k2 = 120, maxprice=200, minprice=0,c1=10,c2=7,longcall=FALSE,shortcall=FALSE,profit=FALSE)
  rp.tkrplot(panel=my.panel , name=tkrp, plotfun=my.draw, pos="right",hscale=2,vscale=2)
  rp.checkbox(my.panel,variable=longcall, action = my.redraw, title = "Long Call")
  rp.checkbox(my.panel,variable=shortcall, action = my.redraw, title = "Short Call")
  rp.checkbox(my.panel,variable=profit, action = my.redraw, title = "Profit")
}
