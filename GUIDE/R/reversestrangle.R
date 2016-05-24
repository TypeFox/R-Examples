if (getRversion() >= "2.15.1") utils::globalVariables(c('tkrp',
                                                        'shortput',
                                                        'shortcall',
                                                        'profit'))

reversestrangle <-
function(){
  
  
  my.draw <- function(panel) 
  {
    
    with(panel, {
      ymin = k1 - maxprice
      ymax =  (c1 + p1)
      
      plot(c(minprice,k1,k2,maxprice),c(0,0,0,0),
           xlab = "Stock price",
           ylab = "Profit",type='l',ylim=c(ymin,ymax))
      drawprofits(panel)
      title(paste("Reverse Strangle with call and put"))
    })
    panel
  }
  
  
  drawprofits<-function(object)  {
    with(object,{
      
      if (shortput){
        S = c(minprice,k1,k2,maxprice)
        profits1<- -pmax(k1-S,0)+p1
        lines(S,profits1,type='l',col="blue", lwd=2)
      }
      
      if (shortcall){
        S = c(minprice,k1,k2,maxprice)
        profits2<- -pmax(S-k2,0)+c1
        lines(S,profits2,type='l',col="red", lwd=2)  
      }
      
      if (profit){
        if((!shortput) | (!shortcall)){
          rp.messagebox("Check the both Short Put and Short Call checkboxes to see profit graph.", title = "Insufficient Information")
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
  
  
  my.panel <- rp.control(title="Reverse Strangle with call and put", k1= 80, k2 = 120, maxprice=200, minprice=0,p1=10,c1=12,shortput=FALSE,shortcall=FALSE,profit=FALSE)
  rp.tkrplot(panel=my.panel , name=tkrp, plotfun=my.draw, pos="right",hscale=2,vscale=2)
  rp.checkbox(my.panel,variable=shortput, action = my.redraw, title = "Short Put")
  rp.checkbox(my.panel,variable=shortcall, action = my.redraw, title = "Short Call")
  rp.checkbox(my.panel,variable=profit, action = my.redraw, title = "Profit")
}
