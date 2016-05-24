if (getRversion() >= "2.15.1") utils::globalVariables(c( 'tkrp',
                                                         'tkrp',
                                                         'couprate',
                                                         'discrate'))

pricematurity <-
function(){
  
  maxmaturity <- 30 # lets keep this fixed
  faceval<-1000 # lets keep this fixed
  
  
  bondval<-function(panel)
  {
    
    maturity <- 0:maxmaturity
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
    effrate = panel$discrate/(100*freq)
    effperiods = freq*maturity
    pv_coupons<-(coupon/effrate)*(1-(1+effrate)^(-effperiods)) # PV of coupons
    pv_face<-faceval*(1+effrate)^(-effperiods) # PV of face value
    price<-pv_coupons+pv_face # bond price is the sum of both
    return (price)
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
  plot(c(0:maxmaturity), bondval(panel),
       xlab = "Maturity(Years)",
       ylab = "Bond Price",,type='l',lwd=2)
  title(paste("Par Value = ", faceval, "\nCoupon = ",couprate,", Disc. Rate = ", discrate,"\n"))
})
    panel
  }
  
  panel <- rp.control(" Price- Maturity Relationship ", frequency="quarterly",couprate= 8, discrate = 10, ylim=2000)
  rp.tkrplot(panel=panel , name=tkrp, plotfun=my.draw, pos="right",hscale=2,vscale=2)
  rp.slider(panel, discrate, 1,30,my.redraw, title = "Discount Rate")
  rp.slider(panel,couprate,1,30,my.redraw, title = "Coupon Rate")
  rp.radiogroup(panel = panel, variable= frequency,
                vals = c("quarterly", "semi-annual", "annual"),
                action = my.redraw, title = "Coupon frequency")
  #rp.do(my.panel, my.draw)
}
