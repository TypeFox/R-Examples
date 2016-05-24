if (getRversion() >= "2.15.1") utils::globalVariables(c("tkrp","facevalue","dur","conv","yldchange"))

bondchange <-
function(){
  
  
  
  my.draw <- function(panel) 
  {
    
    
    faceval<-as.numeric(panel$facevalue)
    dur= as.numeric(panel$dur)
    conv <- as.numeric(panel$conv)
    yldchange <- as.numeric(panel$yldchange)/10000
    
    
    
    #     effrate = discrate/(100*freq)
    #     effperiods = freq*maturity
    #     pv_coupons<-(coupon/effrate)*(1-(1+effrate)^(-effperiods)) # PV of coupons
    #     pv_face<-faceval*(1+effrate)^(-effperiods) # PV of face value
    #     price<-pv_coupons+pv_face # bond price is the sum of both
    #     price <- round(price,2)
    
    
    
    if(panel$approx =="Duration & Convexity"){
      change <- (-dur*yldchange + 0.5*conv*(yldchange)^2)*faceval
    }
    else{
      change <- -dur*yldchange*faceval  
    }
    
    change<- round(change,4)
    
    plot(1:10, 1:10, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(5, 5, paste("Change: ", change),cex=1.4)
    #cat(pv)
    panel
    
  }
  
  
  my.redraw <- function(panel)
  {
    rp.tkrreplot(panel, tkrp)
    panel                                                                       
  }
  
  
  
  
  my.panel <- rp.control("Bond Price Change")
  rp.textentry(panel = my.panel, variable= facevalue,
               labels = "Face Value:        ", action = my.redraw, initval="1000")
  rp.textentry(panel = my.panel, variable= dur,
               labels = "Mod. Duration: ", action = my.redraw, initval="2")
  rp.textentry(panel = my.panel, variable= conv,
               labels = "Convexity:         ", action = my.redraw, initval="8")
  rp.doublebutton(my.panel,variable=yldchange,step=1,title="Change in Yield (bp)",initval=1,range=c(1,30),showvalue=TRUE,action=my.redraw)
  rp.radiogroup(panel = my.panel, variable= approx,
                vals = c("Duration only", "Duration & Convexity"),
                action = my.redraw, title = "Formula/Approximation")
  rp.tkrplot(panel=my.panel , name=tkrp, plotfun=my.draw)
  
  #rp.do(my.panel, my.draw)
}
