if (getRversion() >= "2.15.1") utils::globalVariables(c("S","K","r","div","opttype"))

blackscholes <-
function(){
  my.draw <- function(panel) {
    
    S <-as.numeric(panel$S)
    K <-as.numeric(panel$K)
    r <-as.numeric(panel$r)
    t <-as.numeric(panel$t)
    sigma <-as.numeric(panel$sigma)
    div <-as.numeric(panel$div)
    opttype <- panel$opttype
    
    
    d1 = (log(S/K) + (r -div +  sigma^2/2)*t)/(sigma*sqrt(t))
    #d2 = d1 - sigma*sqrt(t)
    d2 = (log(S/K) + (r -div -  sigma^2/2)*t)/(sigma*sqrt(t))
    Nd1 = pnorm(d1)
    Nd2 = pnorm(d2)
    
   
     if (opttype == "Call"){
      price = S*Nd1*exp(-div*t) - K*exp(-r*t)*Nd2  
    }
    else{
      price =   K*exp(-r*t)*(1-Nd2) -  S*(1-Nd1)*exp(-div*t) 
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
  
  my.panel <- rp.control(title = "Black Scholes price")
  rp.textentry(panel=my.panel,variable=S,title="Spot:            ",action=my.redraw,initval=100)
  rp.textentry(panel=my.panel,variable=K,title="Strike:          ",action=my.redraw,initval=110)
  rp.textentry(panel=my.panel,variable=r,title="Risk free:     ",action=my.redraw,initval=0.05)
  rp.textentry(panel=my.panel,variable=t,title="Maturity:     ",action=my.redraw,initval=0.5)
  rp.textentry(panel=my.panel,variable=sigma,title="Sigma:         ",action=my.redraw,initval=0.30)
  rp.textentry(panel=my.panel,variable=div,title="Div yield:     ",action=my.redraw,initval=0.0)
  rp.radiogroup(panel = my.panel, variable= opttype, vals = c("Call", "Put"), 
                action = my.redraw, title = "Type of Option",initval="Call")
  #   rp.button(panel=my.panel,title="calculate", action=my.redraw)
  rp.tkrplot(panel = my.panel, name = my.tkrplot, plotfun = my.draw)
  
  #rp.do(my.panel, my.draw)
  
   
}
