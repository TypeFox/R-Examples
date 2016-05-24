if (getRversion() >= "2.15.1") utils::globalVariables(c(  'Notional',
                                                          'fixed',
                                                          'float0',
                                                          'SpotRate',
                                                          'spotfreq',
                                                          'settlement'))

irswapvalue <-
function(){
  
  my.draw <- function(panel) {
    Notional = as.numeric(panel$Notional)
    fixedrate = as.numeric(panel$fixed)
    fixedpayment = Notional*fixedrate
    lastfloat = as.numeric(panel$float0)
    floatpayment = lastfloat*Notional
    floatingleg = Notional + Notional*lastfloat
    spotrate = as.numeric(strsplit(panel$SpotRate,",")[[1]])
    spotrates = spotrate[!is.na(spotrate)]
    #     spotrates = as.numeric(spotrates)
    npayments = length(spotrates)
    start = as.numeric(panel$start)/12
    
    fixedpayments = rep(fixedpayment,times=npayments)
    
    #   count number of non NA values
    #   create a vector of fixedpayments of that length 
    
    if (panel$settlement=="quarterly"){
      Effperiods = seq(from=start,by=0.25,length.out=npayments) 
      fixedpayments = fixedpayments/4
      floatpayment = floatpayment/4
    }
    else if (panel$settlement=="semi-annual"){
      Effperiods = seq(from=start,by=0.5,length.out=npayments) 
      fixedpayments = fixedpayments/2
      floatpayment = floatpayment/2
    }
    else{
      Effperiods = seq(from=start,by=1,length.out=npayments) 
    }
    
    
    
    if (panel$spotfreq=="continuous")
    {
      Discrates = exp(-spotrates*Effperiods)
    }
    else {
      Discrates = (1+spotrates)^(-Effperiods)
      
    }
    
    #     val = npayments 
    fixval=sum(Discrates*fixedpayments) + Notional*Discrates[length(Discrates)]
    floatval= (Notional + floatpayment)*Discrates[1]
    val=abs(floatval-fixval)
    val <- round(val, 0)
    
    plot(1:30, 1:30, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(15, 15, paste("Value of Swap = ", val, sep=""),cex=1.2)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "Value of Interest Rate Swap")
  rp.textentry(panel=my.panel,variable=Notional,title="Notional:                           ",action=my.redraw,initval=1000000)
  rp.textentry(panel=my.panel,variable=fixed,title="Fixed Rate:                         ",action=my.redraw,initval=0.06)
  rp.textentry(panel=my.panel,variable=float0,title="Last spot rate:                   ",action=my.redraw,initval=0.05)
  rp.textentry(panel=my.panel,variable=start,title="Months for 1st payment:",action=my.redraw,initval=3)
  rp.textentry(panel=my.panel,variable=SpotRate,title="Spot Rates:                        ",action=my.redraw,initval="0.054, 0.056, 0.058")
  rp.radiogroup(panel = my.panel, variable= spotfreq,
                vals = c("continuous", "quarterly", "semi-annual", "annual"),
                action = my.redraw, title = "Frequency of spot rates")
  rp.radiogroup(panel = my.panel, variable= settlement,
                vals = c("quarterly", "semi-annual", "annual"),
                action = my.redraw, title = "Settlement frequency")
  rp.tkrplot(panel = my.panel,name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
  
  # specify rates with commas
  # times =  seq based on frequency
  # sum(exp(rates * times))
}
