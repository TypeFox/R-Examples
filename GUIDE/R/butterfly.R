if (getRversion() >= "2.15.1") utils::globalVariables(c("tkrp","longcall1","longcall2","shortcalls","profit"))

butterfly <-
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
      title(paste("Butterfly"))
    })
    panel
  }
  
  
  drawprofits<-function(object)  {
    with(object,{
      
      if (longcall1){
        S = c(minprice,k1,kmid,k2,maxprice)
        profits1<- pmax(S-k1,0)-c1
        lines(S,profits1,type='l',col="blue", lwd=2)
      }
      
      if (longcall2){
        S = c(minprice,k1,kmid,k2,maxprice)
        profits2<- pmax(S-k2,0)-c2
        lines(S,profits2,type='l',col="red", lwd=2)  
      }
      
      if (shortcalls){
        S = c(minprice,k1,kmid,k2,maxprice)
        profitsmid<- -2*(pmax(S-kmid,0))+2*cmid
        lines(S,profitsmid,type='l',col="green",lwd=2)  
      }
      
      if (profit){
        if((!longcall1) | (!longcall2) | (!shortcalls)){
          rp.messagebox("Check the three checkboxes Long Call 1,Long Call 2 and Two short Calls to see profit graph.", title = "Insufficient Information")
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
  
  
  my.panel <- rp.control(title="Butterfly", k1= 80, k2 = 120, kmid=100, 
                         maxprice=200, minprice=0,c1=20,c2=10,cmid=13, 
                         longcall1=FALSE,longcall2=FALSE,shortcalls=FALSE, profit=FALSE)
  rp.tkrplot(panel=my.panel , name=tkrp, plotfun=my.draw, pos="right",hscale=2,vscale=2)
  rp.checkbox(my.panel,variable=longcall1, action = my.redraw, title = "Long Call 1")
  rp.checkbox(my.panel,variable=longcall2, action = my.redraw, title = "Long Call 2")
  rp.checkbox(my.panel,variable=shortcalls, action = my.redraw, title = "Short two Calls")
  rp.checkbox(my.panel,variable=profit, action = my.redraw, title = "Profit")
}
