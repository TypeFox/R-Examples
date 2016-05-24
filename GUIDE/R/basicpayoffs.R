if (getRversion() >= "2.15.1") utils::globalVariables(c("tkrp","position","opttype","plottype"))

basicpayoffs <-
function(){
  
  
  my.draw <- function(panel) 
  {
    
    
    with(panel, {
      ymax = maxprice - k1
      ymin=  k2 - maxprice
      
      plot(c(minprice,k1,k2,maxprice),c(0,0,0,0),
           xlab = "Stock price",  
           ylab = "Profit / Payoff",
           type='l',ylim=c(ymin,ymax))
      
      
      
      titletail = ""
      # position <- as.character(panel$position)
      
      
      
      #         minprice <- panel$minprice
      #         maxprice <- panel$maxprice
      #         plottype <- panel$plottype
      #         k1 <- panel$k1
      #         k2 <- panel$k2
      #         c1 <- panel$c1
      #         p1 <- panel$p1
      #         opttype <- panel$opttype
      #         position <- panel$position
      
      opttype<- as.character(opttype)
      
      if (opttype=="Call"){ 
        
        
        if (plottype == "Payoff"){
          titlehead = paste("Payoffs","\n", "K= ", k1) 
          c1=0; p1=0
        }
        else {
          titlehead = paste("Profit & Loss", "\n", "K= ", k1) 
          titletail = paste(", C= ",c1)
        }
        
        
        if (panel$position[1]){
          
          #  position = "Long Call"
          S = c(minprice,k1,k2,maxprice)
          profits<- pmax(S-k1,0)-c1
          lines(S,profits,type='l',col="blue", lwd=2)
          
          title(paste(titlehead,titletail ))
        }
        
        if(panel$position[2]){
          
          #  position = "Short Call"
          S = c(minprice,k1,k2,maxprice)
          profits<- -pmax(S-k1,0)+c1
          lines(S,profits,type='l',col="red", lwd=2) 
          
          title(paste(titlehead,titletail ))
        }
        
      }
      
      
      else  { # it is a put
        
        
        if (plottype == "Payoff"){
          titlehead = paste("Payoff","\n", "K= ", k2) 
          c1=0; p1=0
        }
        else {
          titlehead = paste("Profit & Loss", "\n", "K= ", k2) 
          titletail = paste(", P= ",p1)
        }
        
        
        if (panel$position[1]){
          
          
          #   position = "Long Put"
          S = c(minprice,k1,k2,maxprice)
          profits<- pmax(k2-S,0)-p1
          lines(S,profits,type='l',col="blue", lwd=2)
          
          title(paste(titlehead,titletail ))
        }
        
        if(panel$position[2]) {
          
          
          #   position = "Short Put"
          S = c(minprice,k1,k2,maxprice)
          profits<- -pmax(k2-S,0)+p1
          lines(S,profits,type='l',col="red", lwd=2)  
          
          title(paste(titlehead,titletail ))
        }
        
      }
      
    })
    panel
  }
  
  
  
  
  my.redraw <- function(panel)
  {
    rp.tkrreplot(panel, tkrp)
    panel                                                                       
  }
  
  
  my.panel <- rp.control(title=paste("Payoff / Profit & Loss Graphs"), plottype="Payoff",opttype="Call",k1= 100, k2 = 100, maxprice=200, minprice=0,c1=15,p1=12)
  rp.checkbox(my.panel,variable=position, labels=c("Long","Short"), action = my.redraw, title = "Position")
  rp.tkrplot(panel=my.panel, name=tkrp, plotfun=my.draw, hscale=1.5, vscale=1.5, pos="right")
  rp.radiogroup(my.panel,variable=opttype,labels=c("Call","Put"),vals=c("Call","Put"),action=my.redraw,title="Option Type")
  rp.radiogroup(my.panel,variable=plottype,title="Plot Type", labels=c("Payoff","Profit & Loss"),vals=c("Payoff","Profit & Loss"),action=my.redraw)
  #rp.do(panel=my.panel,action=my.draw)
}
