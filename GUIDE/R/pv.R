if (getRversion() >= "2.15.1") utils::globalVariables(c( 'FV',
                                                         'r'))
pv <-
function(){
  my.draw <- function(panel) {
    
    if(panel$frequency=="continuous"){
      pv<-as.numeric(panel$FV)*exp(-1*as.numeric(panel$r)*as.numeric(panel$t))
    }
    else if(panel$frequency=="quarterly"){
      effrate<-as.numeric(panel$r)/(4)
      effperiod<-(as.numeric(panel$t)*4)
      pv<-as.numeric(panel$FV)/((1+effrate)^(effperiod))
    }
    else if(panel$frequency=="semi-annual"){
      effrate<-as.numeric(panel$r)/(2)
      effperiod<-(as.numeric(panel$t)*2)
      pv<-as.numeric(panel$FV)/((1+effrate)^(effperiod))
    }
    else{
      effrate<-as.numeric(panel$r)/(1)
      effperiod<-(as.numeric(panel$t)*1)
      pv<-as.numeric(panel$FV)/((1+effrate)^(effperiod))
    }
    pv<-round(pv,2)
    plot(1:10, 1:10, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(5, 5, paste("PV: ", pv),cex=1.4)
    #cat(pv)
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "Present Value")
  rp.textentry(panel = my.panel, variable= FV,
               labels = "Future Value: ", action = my.redraw, initval="100")
  rp.textentry(panel = my.panel, variable= r,
               labels = "Rate:              ", action = my.redraw, initval="0.10")
  rp.textentry(panel = my.panel, variable= t,
               labels = "Time:             ", action = my.redraw, initval="1")
  rp.radiogroup(panel = my.panel, variable= frequency,
                vals = c("continuous", "quarterly", "semi-annual", "annual"),
                action = my.redraw, title = "Compounding frequency")
  rp.tkrplot(panel = my.panel, name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
}
