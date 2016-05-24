if (getRversion() >= "2.15.1") utils::globalVariables(c('t1',
                                                        'r1',
                                                        't2',
                                                        'r2'))

fra <-
function(){
  my.draw <- function(panel) {
    r1 <-as.numeric(panel$r1)
    r2 <-as.numeric(panel$r2)
    quote = paste(panel$t1,"x", panel$t2)
    t1 <-as.numeric(panel$t1)
    t2 <-as.numeric(panel$t2)
    
    
    if (panel$frequency=="continuous")
    {
      t2=t2/12
      t1=t1/12
      f12 = (r2*t2-r1*t1)/(t2-t1)
      f12 = round(12/(t2-t1)*(exp(f12*(t2-t1)/12)-1), 4)
      #    fwd rates are expressed in comp frequency of t2-t1   
    }
    else{
      r2=r2*(t2-t1)/12
      r1=r1*(t2-t1)/12
      t2new=t2/(t2-t1)
      t1new=t1/(t2-t1)
      t12=1
      f12 <- round(12/(t2-t1)*(((1+r2)^t2new)/((1+r1)^t1new) - 1), 4)
    }
    
    plot(1:30, 1:30, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(15, 15, paste(quote, " Fwd Rate = ", f12, sep=""),cex=1.2)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "Forward Rate")
  rp.textentry(panel=my.panel,variable=t1,title="Months1:",action=my.redraw,initval=3)
  rp.textentry(panel=my.panel,variable=r1,title="Rate1:      ",action=my.redraw,initval=0.09)
  rp.textentry(panel=my.panel,variable=t2,title="Months2:",action=my.redraw,initval=6)
  rp.textentry(panel=my.panel,variable=r2,title="Rate2:      ",action=my.redraw,initval=0.12)
  rp.radiogroup(panel = my.panel, variable= frequency,
                vals = c("Continuous", "Loan period"),
                action = my.redraw, title = "Frequency of spot rates")
  rp.tkrplot(panel = my.panel, pos="bottom",name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
  
}
