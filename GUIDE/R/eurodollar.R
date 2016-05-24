eurodollar <-
function(){
  
  my.draw <- function(panel) {
    quote <-as.numeric(panel$quote)
    
    val <- round(10000*(100-0.25*(100-quote)), 2)
    
    
    plot(1:30, 1:30, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(15, 15, paste("Value = ", val, sep=""),cex=1.5)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "Eurodollar futures contract value")
  rp.textentry(panel=my.panel,variable=quote,title="CME Quote:",action=my.redraw,initval=97.8)
  rp.tkrplot(panel = my.panel, pos="bottom",name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
  
}
