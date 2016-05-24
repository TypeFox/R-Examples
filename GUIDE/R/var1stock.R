if (getRversion() >= "2.15.1") utils::globalVariables(c('price',
                                                        'conf',
                                                        'horizon',
                                                        'distribution'))

var1stock <-
function(){

  my.draw <- function(panel) {
    price <-as.numeric(panel$price)
    mu <- as.numeric(panel$mu)
    sigma <- as.numeric(panel$sigma)
    conf <- as.numeric(panel$conf)
    horizon <- as.numeric(panel$horizon)/12
    
   if (panel$distribution=="normal"){
     var= (qnorm(conf)*sigma*sqrt(horizon)-mu)*price
   }
    else{
     var= (1-(exp(mu-qnorm(conf)*sigma*sqrt(horizon)))) *price
    }
    var<-round(var,2)
    plot(1:20, 1:20, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(10, 10, paste("VaR= ", var, sep=""),cex=1.5)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "Value at Risk- single stock/Portfolio")
  rp.textentry(panel=my.panel,variable=price,title="Value:          ",action=my.redraw,initval=120)
  rp.textentry(panel=my.panel,variable=mu,title="mu:             ",action=my.redraw,initval=0.10)
  rp.textentry(panel=my.panel,variable=sigma,title="sigma:        ",action=my.redraw,initval=0.30)
  rp.textentry(panel=my.panel,variable=conf,title="Conf level: ",action=my.redraw,initval=0.95)
  rp.doublebutton(my.panel,variable=horizon,step=1,title="Horizon (Months)",initval=12,showvalue=T,range=c(1,12),action=my.redraw)
#   rp.slider(panel=my.panel,variable=horizon,from=0.25,to=1,resolution=0.1,title = "Horizon:",action=my.redraw,showvalue=T,initval=1)
  rp.radiogroup(panel = my.panel, variable= distribution,
                vals = c("normal", "log-normal"),
                action = my.redraw, title = "Distribution")
  
  rp.tkrplot(panel = my.panel, name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
}
