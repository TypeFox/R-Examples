if (getRversion() >= "2.15.1") utils::globalVariables(c('tkrp',
                                                        'shortcall1',
                                                        'shortcall2',
                                                        'longcalls',
                                                        'profit'))

reversebutterfly <-
function(){
  
  
  my.draw <- function(panel) 
  {
    
    with(panel, {
      ymin = minprice - k1
      ymax =  maxprice - k2
      
      plot(c(minprice,k1,kmid,k2,maxprice),c(0,0,0,0,0),
           xlab = "Stock price",
           ylab = "Profit",type='l',ylim=c(ymin,ymax))
      drawprofits(panel)
      title(paste("Reverse Butterfly"))
    })
    panel
  }
  
  
  drawprofits<-function(object)  {
    with(object,{
      
      if (shortcall1){
        S = c(minprice,k1,kmid,k2,maxprice)
        profits1<- -pmax(S-k1,0)+c1
        lines(S,profits1,type='l',col="blue", lwd=2)
      }
      
      if (shortcall2){
        S = c(minprice,k1,kmid,k2,maxprice)
        profits2<- -pmax(S-k2,0)+c2
        lines(S,profits2,type='l',col="red", lwd=2)  
      }
      
      if (longcalls){
        S = c(minprice,k1,kmid,k2,maxprice)
        profitsmid<- 2*(pmax(S-kmid,0))-2*cmid
        lines(S,profitsmid,type='l',col="green",lwd=2)  
      }
      
      if (profit){
        if((!shortcall1) | (!shortcall2) | (!longcalls)){
          rp.messagebox("Check the three checkboxes Short Call 1,Short Call 2 and Two long Calls to see profit graph.", title = "Insufficient Information")
        }
        else{
          S = c(minprice,k1,kmid,k2,maxprice)
          profits<- profits1+profits2+profitsmid
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
  
  
  my.panel <- rp.control(title="Reverse Butterfly", k1= 80, k2 = 120, kmid=100, 
                         maxprice=200, minprice=0,c1=20,c2=10,cmid=13, 
                         shortcall1=FALSE,shortcall2=FALSE,longcalls=FALSE, profit=FALSE)
  rp.tkrplot(panel=my.panel , name=tkrp, plotfun=my.draw, pos="right",hscale=2,vscale=2)
  rp.checkbox(my.panel,variable=shortcall1, action = my.redraw, title = "Short Call 1")
  rp.checkbox(my.panel,variable=shortcall2, action = my.redraw, title = "Short Call 2")
  rp.checkbox(my.panel,variable=longcalls, action = my.redraw, title = "Long two Calls")
  rp.checkbox(my.panel,variable=profit, action = my.redraw, title = "Profit")
}
