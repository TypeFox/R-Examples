if (getRversion() >= "2.15.1") utils::globalVariables(c("tkrp","facevalue","couprate","discrate","maturity","ratefreq","durtype"))

bonddur <-
function(){
  
  
  
  my.draw <- function(panel) 
  {
    
    
    faceval<-as.numeric(panel$facevalue)
    discrate = as.numeric(panel$discrate)/100
    maturity <- panel$maturity
    
    if (panel$frequency == "quarterly"){
      freq<-4 # change thru radio button
      times<-seq(from=0.25,by=0.25,length.out=maturity*freq)
    }
    else if (panel$frequency == "semi-annual"){
      freq<-2 # change thru radio button
      times<-seq(from=0.5,by=0.5,length.out=maturity*freq)
    }
    else {
      freq<-1 # change thru radio button
      times<-seq(from=1,by=1,length.out=maturity*freq)
    }
    
    
    
    
    if (panel$ratefreq=="continuous comp"){
      pvfactors=exp(-discrate*times) 
    }
    else if(panel$ratefreq=="annual comp"){
      pvfactors=1/(1+discrate)^times
    }
    else{
      pvfactors=1/(1+discrate/freq)^(freq*times)
    }
    
    
    #     effrate = discrate/(100*freq)
    #     effperiods = freq*maturity
    #     pv_coupons<-(coupon/effrate)*(1-(1+effrate)^(-effperiods)) # PV of coupons
    #     pv_face<-faceval*(1+effrate)^(-effperiods) # PV of face value
    #     price<-pv_coupons+pv_face # bond price is the sum of both
    #     price <- round(price,2)
    
    coupon<-panel$couprate*faceval/(100*freq)
    cashflows <- rep(coupon,maturity*freq)
    cashflows[length(cashflows)] = cashflows[length(cashflows)]+faceval
    
    price<-sum(cashflows*pvfactors)
    dur = sum(cashflows*pvfactors*times)/price
    
    if (panel$durtype=="Modified"){
      dur<- dur/(1+discrate/freq)  
    }
    
    dur<-round(dur,2)
    
    plot(1:10, 1:10, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(5, 5, paste("Duration: ", dur),cex=1.4)
    #cat(pv)
    panel
    
  }
  
  
  my.redraw <- function(panel)
  {
    rp.tkrreplot(panel, tkrp)
    panel                                                                       
  }
  
  
  
  
  my.panel <- rp.control("Bond Duration", frequency="quarterly",couprate= 8,discrate=10, maturity = 10)
  rp.textentry(panel = my.panel, variable= facevalue,
               labels = "Face Value:   ", action = my.redraw, initval="1000")
  rp.doublebutton(my.panel,variable=couprate,step=0.25,title="Coupon (%  p.a.)",initval=10,range=c(1,15),showvalue=TRUE,action=my.redraw)
  rp.doublebutton(my.panel,variable=discrate,step=0.25,title="Discount Rate (%  p.a.)",initval=10,range=c(1,15),showvalue=TRUE,action=my.redraw)
  rp.doublebutton(my.panel,variable=maturity,step=0.25,title="Maturity (Yrs)",initval=10,range=c(1,25),showvalue=TRUE,action=my.redraw)
  rp.radiogroup(panel = my.panel, variable= frequency,
                vals = c("quarterly", "semi-annual", "annual"),
                action = my.redraw, title = "Coupon payments")
  rp.radiogroup(panel = my.panel, variable= ratefreq,
                vals = c("continuous comp", "same as coupon freq","annual comp"),
                action = my.redraw, title = "Frequency of discount rate")
  rp.radiogroup(panel = my.panel, variable= durtype,
                vals = c("Macaulay", "Modified"),
                action = my.redraw, title = "Duration formula")
  rp.tkrplot(panel=my.panel , name=tkrp, plotfun=my.draw)
  
  #rp.do(my.panel, my.draw)
}
