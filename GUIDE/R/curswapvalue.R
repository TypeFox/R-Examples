if (getRversion() >= "2.15.1") utils::globalVariables(c('Notional1',                                                    
                                                        'fixed1',
                                                        'homerate',
                                                        'Notional2',
                                                        'fixed2',
                                                        'foreignrate',
                                                        'SpotRate',
                                                        'nperiods',
                                                        'settlement'))

curswapvalue <-
function(){
  
  my.draw <- function(panel) {
    Notional1 = as.numeric(panel$Notional1)
    Notional2 = as.numeric(panel$Notional2)
    fixedrate1 = as.numeric(panel$fixed1)
    fixedrate2 = as.numeric(panel$fixed2)
    fixedpayment1 = Notional1*fixedrate1
    fixedpayment2 = Notional2*fixedrate2
    spotrate = as.numeric(panel$SpotRate)
    nperiods = as.numeric(panel$nperiods)
#   spotrates = spotrates[is.na(spotrates)== F]
#   spotrates = as.numeric(spotrates)
    
    start = as.numeric(panel$start)/12
   homerate = as.numeric(panel$homerate)
  foreignrate = as.numeric(panel$foreignrate)  
    
    
    fixedpayments1 = rep(fixedpayment1,times=nperiods)
    fixedpayments2 = rep(fixedpayment2,times=nperiods)
    #   count number of non NA values
    #   create a vector of fixedpayments of that length 
    
    if (panel$settlement=="quarterly"){
      Effperiods = seq(from=start,by=0.25,length.out=nperiods) 
      fixedpayments1 = fixedpayments1/4
      fixedpayments2 = fixedpayments2/4
    }
    else if (panel$settlement=="semi-annual"){
      Effperiods = seq(from=start,by=0.5,length.out=nperiods) 
      fixedpayments1 = fixedpayments1/2
      fixedpayments2 = fixedpayments2/2
    }
    else{
      Effperiods = seq(from=start,by=1,length.out=nperiods) 
    }
    
    
#   assuming spot rates are on continuous comp basis  
    Discrates1 = exp(-homerate*Effperiods)
    Discrates2 = exp(-foreignrate*Effperiods)
     
    #     val = npayments 
    fixval1=sum(Discrates1*fixedpayments1) + Notional1*Discrates1[length(Discrates1)]
    fixval2=sum(Discrates2*fixedpayments2) + Notional2*Discrates2[length(Discrates2)]
    
    val=abs(fixval1-fixval2*spotrate)
    val <- round(val, 2)
    
    plot(1:30, 1:30, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(15, 15, paste("Value of Swap = ", val, sep=""),cex=1.2)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "Value of Currency Swap")
  rp.textentry(panel=my.panel,variable=Notional1,title="Notional (Home):                ",action=my.redraw,initval=175)
  rp.textentry(panel=my.panel,variable=fixed1,title="Payment Rate (Home):       ",action=my.redraw,initval=0.05)
  rp.textentry(panel=my.panel,variable=homerate,title="Interest Rate (Home):         ",action=my.redraw,initval=0.02)
  rp.textentry(panel=my.panel,variable=Notional2,title="Notional (Foreign):             ",action=my.redraw,initval=100)
  rp.textentry(panel=my.panel,variable=fixed2,title="Payment Rate (Foreign):    ",action=my.redraw,initval=0.06)
  rp.textentry(panel=my.panel,variable=foreignrate,title="Interest Rate (Foreign):       ",action=my.redraw,initval=0.04)
  rp.textentry(panel=my.panel,variable=start,title="Months for 1st payment:   ",action=my.redraw,initval=12)
  rp.textentry(panel=my.panel,variable=SpotRate,title="Spot Exchange rate:           ",action=my.redraw,initval=1.5)
  rp.textentry(panel=my.panel,variable=nperiods,title="Number of Periods:            ",action=my.redraw,initval=3)
  rp.radiogroup(panel = my.panel, variable= settlement,
                vals = c("quarterly", "semi-annual", "annual"),
                action = my.redraw, title = "Settlement frequency")
  rp.tkrplot(panel = my.panel, pos="bottom",name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
}
