if (getRversion() >= "2.15.1") utils::globalVariables(c( 'tkrp',
                                                         'tkrp',
                                                         'couprate',
                                                         'maturity'))

priceyield <-
function(){
  
  #   maxmaturity <- 30 # lets keep this fixed
  faceval<-1000 # lets keep this fixed
  discrates = seq(0,30,1)
  prices = seq(0,0.30,length.out=length(discrates))
  
  
  bondvalues<-function(panel)
  {
    
    maturity <- panel$maturity
    if (panel$frequency == "quarterly"){
      freq<-4 # change thru radio button
    }
    else if (panel$frequency == "semi-annual"){
      freq<-2 # change thru radio button
    }
    else {
      freq<-1 # change thru radio button
    }
    coupon<-panel$couprate*faceval/(100*freq)
    
    for (i in 1:length(discrates)){
      effrate = discrates[i]/(100*freq)
      effperiods = freq*maturity
      pv_coupons<-(coupon/effrate)*(1-(1+effrate)^(-effperiods)) # PV of coupons
      pv_face<-faceval*(1+effrate)^(-effperiods) # PV of face value
      prices[i]<-pv_coupons+pv_face # bond price is the sum of both
    }
    return (prices)
  }
  
  
  
  my.redraw <- function(panel)
  {
    rp.tkrreplot(panel, tkrp)
    panel                                                                       
  }
  
  
  my.draw <- function(panel) 
  {
    with(panel, 
{
  plot(discrates, bondvalues(panel),
       xlab = "Yield(%)",
       ylab = "Bond Price",,type='l',lwd=2)
  title(paste("Par Value = ", faceval, "\nCoupon = ",couprate,", Maturity = ", maturity,"\n"))
})
    panel
  }
  
  panel <- rp.control(" Price- yield Relationship ", frequency="quarterly",couprate= 8, maturity = 10)
#   rp.slider(panel, discrate, 1,30,my.redraw, title = "Discount Rate")
  rp.tkrplot(panel=panel , name=tkrp, plotfun=my.draw,pos="right",hscale=2,vscale=2)
  rp.slider(panel,couprate,1,20,my.redraw, title = "Coupon Rate")
  rp.doublebutton(panel,variable=maturity,step=3,title="Maturity",initval=10,range=c(1,22),action=my.redraw)
  rp.radiogroup(panel = panel, variable= frequency,
                vals = c("quarterly", "semi-annual", "annual"),
                action = my.redraw, title = "Coupon frequency")
  
  
  #rp.do(my.panel, my.draw)
}
