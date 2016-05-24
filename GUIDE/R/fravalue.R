if (getRversion() >= "2.15.1") utils::globalVariables(c('Notional',
                                                        'f12',
                                                        'fixed',
                                                        't1',
                                                        't2'))

fravalue <-
function(){
  my.draw <- function(panel) {
    f12 <-as.numeric(panel$f12)
    fixed <- as.numeric(panel$fixed)
    t1 <-as.numeric(panel$t1)/12
    t2 <-as.numeric(panel$t2)/12
    Notional<-as.numeric(panel$Notional)
    
    val <- abs(round(Notional*(fixed-f12)*(t2-t1), 2))
    
    
    plot(1:30, 1:30, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(15, 15, paste("Value = ", val, sep=""),cex=1.5)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "Value of FRA")
  rp.textentry(panel=my.panel,variable=Notional,title="Notional:  ",action=my.redraw,initval=100000)
  rp.textentry(panel=my.panel,variable=f12,title="Fwd Rate:  ",action=my.redraw,initval=0.09)
  rp.textentry(panel=my.panel,variable=fixed,title="Fixed Rate:",action=my.redraw,initval=0.10)
  rp.textentry(panel=my.panel,variable=t1,title="Months1:    ",action=my.redraw,initval=3)
  rp.textentry(panel=my.panel,variable=t2,title="Months2:    ",action=my.redraw,initval=6)
  rp.tkrplot(panel = my.panel, pos="bottom",name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
  
}
