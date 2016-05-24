if (getRversion() >= "2.15.1") utils::globalVariables(c('fromfrequency',
                                                        'tofrequency',
                                                        'fromrate'))
rate <-
function(){
  my.draw <- function(panel) {
    fromrate <-as.numeric(panel$fromrate)
    
    if(panel$fromfrequency=="continuous"){
      if (panel$tofrequency=="continuous"){
        torate<- fromrate
      }
      else if(panel$tofrequency=="quarterly"){
        torate<- 4*(exp(fromrate/4)-1)
      }
      else if(panel$tofrequency=="semi-annnual"){
        torate<- 2*(exp(fromrate/2)-1)
      }
      else {
        torate<- 1*(exp(fromrate/1)-1)
      }}
    else if (panel$fromfrequency=="quarterly"){      
      if (panel$tofrequency=="continuous"){
        torate<- 4*log(1+ fromrate/4) 
      }
      else if(panel$tofrequency=="quarterly"){
        torate<- fromrate
      }
      else if(panel$tofrequency=="semi-annnual"){
        torate<- 2*((1 + fromrate/4)^(4/2) - 1)
      }
      else {
        torate<- (1 + fromrate/4)^4 - 1
      }} 
    else if (panel$fromfrequency=="semi-annual"){      
      if (panel$tofrequency=="continuous"){
        torate<- 2*log(1+ fromrate/2) 
      }
      else if(panel$tofrequency=="quarterly"){
        torate<- 4*((1 + fromrate/2)^(2/4) - 1)
      }
      else if(panel$tofrequency=="semi-annnual"){
        torate<- fromrate
      }
      else {
        torate<- (1 + fromrate/2)^2 - 1
      }} 
    else {      
      if (panel$tofrequency=="continuous"){
        torate<- 1*log(1+ fromrate/1) 
      }
      else if(panel$tofrequency=="quarterly"){
        torate<- 4*((1 + fromrate)^(1/4) - 1)
      }
      else if(panel$tofrequency=="semi-annnual"){
        torate<- 2*((1 + fromrate)^(1/2) - 1)
      }
      else {
        torate<- fromrate
      }}
    
    torate<-round(torate,4)
    plot(1:20, 1:20, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(10, 10, paste("Rate = ", torate, sep=""),cex=1.5)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "Rate Converter")
  rp.radiogroup(panel = my.panel, variable= fromfrequency,
                vals = c("continuous", "quarterly", "semi-annual", "annual"),
                action = my.redraw, title = "Given frequency")
  rp.radiogroup(panel = my.panel, variable= tofrequency,
                vals = c("continuous", "quarterly", "semi-annual", "annual"),
                action = my.redraw, title = "Required frequency")
  rp.textentry(panel=my.panel,variable=fromrate,title="Given rate: ",action=my.redraw,initval=0.05)
  rp.tkrplot(panel = my.panel, name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
}
