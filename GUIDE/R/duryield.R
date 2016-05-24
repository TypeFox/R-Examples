if (getRversion() >= "2.15.1") utils::globalVariables(c('tkrp',
                                                        'couprate',
                                                        'maturity'))

duryield <-
function(){
  
  
  #discrate<-0.10
  maxdiscrate<-0.20
  #discrate <- 10/ # lets keep this fixed
  faceval<-1000 # lets keep this fixed
  
  
  
  #faceval<-as.numeric(panel$facevalue)
  
  # maturity <- panel$maturity
  
  bonddur<- function(panel,discrate,faceval){
    
    
    couprate<-as.numeric(panel$couprate)/100
    maturity<-as.numeric(panel$maturity)
    #discrate = as.numeric(panel$discrate)/100
    freq<-1 # change thru radio button
    times<-seq(from=1,by=1,length.out=maturity*freq)
    
    pvfactors=1/(1+discrate/freq)^(freq*times)
    
    
    
    #     effrate = discrate/(100*freq)
    #     effperiods = freq*maturity
    #     pv_coupons<-(coupon/effrate)*(1-(1+effrate)^(-effperiods)) # PV of coupons
    #     pv_face<-faceval*(1+effrate)^(-effperiods) # PV of face value
    #     price<-pv_coupons+pv_face # bond price is the sum of both
    #     price <- round(price,2)
    
    coupon<-couprate*faceval/100
    cashflows <- rep(coupon,maturity*freq)
    cashflows[length(cashflows)] = cashflows[length(cashflows)]+faceval
    
    price<-sum(cashflows*pvfactors)
    dur = sum(cashflows*pvfactors*times)/price
    dur<-round(dur,2)
    return(dur)
  }
  
  
  
  
  my.draw <- function(panel) 
  {
    with(panel, {
      bonddurs<-rep(0,20)
      discrates<-seq(1,20,1)
      for (i in 1:20){
        bonddurs[i]=bonddur(panel,discrates[i],faceval)
      }
      
      plot(discrates, bonddurs, type='l',lwd=2, xlab="Yield", ylab="Duration",
           frame = TRUE)
      title(paste("Duration and Yield"))
      #cat(pv)
    })
    panel
    
  }
  
  
  
  
  
  my.redraw <- function(panel)
  {
    rp.tkrreplot(panel, tkrp)
    panel                                                                       
  }
  
  
  
  
  my.panel <- rp.control("Duration and Yield")
  rp.doublebutton(my.panel,variable=couprate,step=0.5,title="Coupon (%  p.a.)",initval=10,range=c(1,15),showvalue=TRUE,action=my.redraw)
  #rp.doublebutton(my.panel,variable=discrate,step=0.25,title="Discount Rate (%  p.a.)",initval=10,range=c(1,15),showvalue=TRUE,action=my.redraw)
  rp.doublebutton(my.panel,variable=maturity,step=1,title="Maturity (Yrs)",initval=10,range=c(1,25),showvalue=TRUE,action=my.redraw)
  rp.tkrplot(panel=my.panel , name=tkrp, plotfun=my.draw)
  
  #rp.do(my.panel, my.draw)
}
