if (getRversion() >= "2.15.1") utils::globalVariables(c("quoted","convfac","accint"))

cashprice <-
function(){
  
  my.draw <- function(panel) {
    quoted <-as.numeric(panel$quoted)
    convfac <-as.numeric(panel$convfac)
    accint <-as.numeric(panel$accint)
    
    val <- round(quoted*convfac+accint, 2)
    
    
    plot(1:30, 1:30, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(15, 15, paste("Value = ", val, sep=""),cex=1.5)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "Cash price of a T Bond Futures")
  rp.textentry(panel=my.panel,variable=quoted,title="Quoted Price:",action=my.redraw,initval=102)
  rp.textentry(panel=my.panel,variable=convfac,title="Conv. Factor:",action=my.redraw,initval=1.08)
  rp.textentry(panel=my.panel,variable=accint,title="Acc. Interest: ",action=my.redraw,initval=4)
  rp.tkrplot(panel = my.panel, pos="bottom",name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
  
}
