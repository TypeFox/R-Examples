zval <-
function(){
  my.draw <- function(panel) {
    pval <-as.numeric(panel$pval)
    
    
    zval<-round(qnorm(pval),4)
    plot(1:20, 1:20, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(10, 10, paste("z-value = ", zval, sep=""),cex=1.5)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "z-value calculator")
  rp.textentry(panel=my.panel,variable=pval,title="p value: ",action=my.redraw,initval=0.10)
  rp.tkrplot(panel = my.panel, name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
  
}
