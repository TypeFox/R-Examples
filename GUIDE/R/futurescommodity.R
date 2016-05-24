if (getRversion() >= "2.15.1") utils::globalVariables(c('S',
                                                        'r',
                                                        'storagetype',
                                                        'storage',
                                                        'storagetime',
                                                        'convenience'))

futurescommodity <-
function(){
  my.draw <- function(panel) {
    
    S <-as.numeric(panel$S)
    r <-as.numeric(panel$r)
    t <-as.numeric(panel$t)
    
    storagetype <- panel$storagetype
    convenience <- as.numeric(panel$convenience)
    
    if (storagetype == "Yield"){
      storage <- as.numeric(panel$storage)
      price <-   (exp(r+storage-convenience)*t)*S
    }
    else {
      storage <- as.numeric(strsplit(as.character(panel$storage),",")[[1]])
      storagetime <- as.numeric(strsplit(as.character(panel$storagetime),",")[[1]])
      price = (S+sum(storage*exp(-r*storagetime)))*(exp(r-convenience)*t)
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
  
  my.panel <- rp.control(title = "Commodity Futures")
  rp.textentry(panel=my.panel,variable=S,title="Spot:                           ",action=my.redraw,initval=100)
  rp.textentry(panel=my.panel,variable=r,title="Risk free:                    ",action=my.redraw,initval=0.05)
  rp.textentry(panel=my.panel,variable=t,title="Maturity:                    ",action=my.redraw,initval=0.5)
  rp.textentry(panel=my.panel,variable=convenience,title="Convenience yield:    ",action=my.redraw,initval=0)
  rp.textentry(panel=my.panel,variable=storage,title="Storage Cost(s):          ",action=my.redraw,initval=0)
  rp.textentry(panel=my.panel,variable=storagetime,title="Storage Cost Time(s):",action=my.redraw,initval=0)
  rp.radiogroup(panel = my.panel, variable= storagetype, vals = c("Yield","Cash"),
                initval="Yield",action = my.redraw, title = "Type of Storage Cost")
  rp.tkrplot(panel = my.panel, name = my.tkrplot, plotfun = my.draw)
  
  #rp.do(my.panel, my.draw)
  
  
}
