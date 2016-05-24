if (getRversion() >= "2.15.1") utils::globalVariables(c('S',
                                                        'r',
                                                        'incometype',
                                                        'income',
                                                        'incometime'))

forwardstock <-
function(){
  my.draw <- function(panel) {
    
    S <-as.numeric(panel$S)
    r <-as.numeric(panel$r)
    t <-as.numeric(panel$t)
    
    incometype <- panel$incometype
    
    
      if (incometype == "Yield"){
      income <-as.numeric(panel$income)
      price =   exp((r-income)*t)*S
    }
    else{
      income <- as.numeric(strsplit(as.character(panel$income),",")[[1]])
      incometime <- as.numeric(strsplit(as.character(panel$incometime),",")[[1]])
      price = (S-sum(income*exp(-r*incometime)))*exp(r*t)
	 }

    
    plot(1:20, 1:20, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(10, 10, paste("Price = ", round(price,3), sep=""),cex=1.5)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "Stock Forward")
  rp.textentry(panel=my.panel,variable=S,title="Spot:                       ",action=my.redraw,initval=100)
  rp.textentry(panel=my.panel,variable=r,title="Risk free:                ",action=my.redraw,initval=0.05)
  rp.textentry(panel=my.panel,variable=t,title="Maturity:                ",action=my.redraw,initval=0.5)
  rp.textentry(panel=my.panel,variable=income,title="Dividend(s):            ",action=my.redraw,initval=0)
  rp.textentry(panel=my.panel,variable=incometime,title="Dividend times(s): ",action=my.redraw,initval=0)
  rp.radiogroup(panel = my.panel, variable= incometype, vals = c("Yield","Cash"), 
                action = my.redraw, title = "Type of Income")
  rp.tkrplot(panel = my.panel, name = my.tkrplot, plotfun = my.draw)
  
  #rp.do(my.panel, my.draw)
  
  
}
