pval <-
function(){
  my.draw <- function(panel) {
    zval <-as.numeric(panel$zval)
    
    
    pval<-round(pnorm(zval),4)
    plot(1:20, 1:20, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(10, 10, paste("p-value = ", pval, sep=""),cex=1.5)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "P-value calculator")
  rp.textentry(panel=my.panel,variable=zval,title="z value: ",action=my.redraw,initval=2.33)
  rp.tkrplot(panel = my.panel, name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
  
}
