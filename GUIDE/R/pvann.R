if (getRversion() >= "2.15.1") utils::globalVariables(c( 'Instalment',
                                                         'r'))

pvann <-
function(){
  my.draw <- function(panel) {
    coupon <- as.numeric(panel$Instalment)
    
    if(panel$frequency=="monthly"){
      effrate<-as.numeric(panel$r)/(12)
      effperiods<-(as.numeric(panel$t)*12)
      pv<-(coupon/effrate)*(1-(1+effrate)^(-effperiods))
    }
    else if(panel$frequency=="quarterly"){
      effrate<-as.numeric(panel$r)/(4)
      effperiods<-(as.numeric(panel$t)*4)
      pv<-(coupon/effrate)*(1-(1+effrate)^(-effperiods))
    }
    else if(panel$frequency=="semi-annual"){
      effrate<-as.numeric(panel$r)/(2)
      effperiods<-(as.numeric(panel$t)*2)
      pv<-(coupon/effrate)*(1-(1+effrate)^(-effperiods))
    }
    else{
      effrate<-as.numeric(panel$r)/(1)
      effperiods<-(as.numeric(panel$t)*1)
      pv<-(coupon/effrate)*(1-(1+effrate)^(-effperiods))
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
  
  my.panel <- rp.control(title = "Present Value of Annuity")
  rp.textentry(panel = my.panel, variable= Instalment,
               labels = "Installment:   ", action = my.redraw, initval="1000")
  rp.textentry(panel = my.panel, variable= r,
               labels = "Rate:              ", action = my.redraw, initval="0.10")
  rp.textentry(panel = my.panel, variable= t,
               labels = "Time:             ", action = my.redraw, initval="1")
  rp.radiogroup(panel = my.panel, variable= frequency,
                vals = c("monthly", "quarterly", "semi-annual", "annual"),
                action = my.redraw, title = "Payment frequency")
  rp.tkrplot(panel = my.panel, name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
}
