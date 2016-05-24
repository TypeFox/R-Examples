if (getRversion() >= "2.15.1") utils::globalVariables(c("S","K","r","div","opttype","greek"))

calcgreeks <-
function(){
  my.draw <- function(panel) {
    
    S <-as.numeric(panel$S)
    K <-as.numeric(panel$K)
    r <-as.numeric(panel$r)
    t <-as.numeric(panel$t)
    sigma <-as.numeric(panel$sigma)
    div <-as.numeric(panel$div)
    opttype <- panel$opttype
    greek <- panel$greek
    
    d1 = (log(S/K) + (r +  sigma^2/2)*t)/(sigma*sqrt(t))
    d2 = d1 - sigma*sqrt(t)
    Nd1 = pnorm(d1)
    Nd2 = pnorm(d2)
    
    
    bsgreek <- function(S,K,r,t,sigma,div,opttype,greek){
      
      if (opttype=="Call"){
        if (greek=="Delta"){
          greekval <- exp( -div*t ) * pnorm(d1)
          
        }
        else if (greek=="Gamma"){
          greekval <- exp(-div*t) * ( dnorm(d1)/(S*sigma*sqrt(t)) )
          
        }
        else if (greek=="Vega"){
          greekval  <- S * exp(-div*t) * dnorm(d1) * sqrt(t)
          
        }
        else if (greek=="Theta"){
          greekval <- exp(-div*t) * (S*dnorm(d1)*sigma)/(2*sqrt(t)) -
            r * K * exp(-r*t) * pnorm(d2)
          
        }
        else {
          greekval  <- K * t * exp(-r*t) * pnorm(d2)
          
        }
        
        
      }
      #       in case of put
      else{
        if (greek=="Delta"){ 
          greekval <- -exp(-div*t ) * pnorm(-d1) 
          
        }
        else if (greek=="Gamma"){
          greekval <- exp(-div*t) * ( dnorm(d1)/(S*sigma*sqrt(t)) )
          
        }
        else if (greek=="Vega"){
          greekval  <- S * exp(-div*t) * dnorm(d1) * sqrt(t)
          
        }
        else if (greek=="Theta"){
          greekval <- exp(-div*t) * (S*dnorm(d1)*sigma)/(2*sqrt(t)) +
            r * K * exp(-r*t) * pnorm(-d2)
          
        }
        else {
          greekval  <- K * t * exp(-r*t) * pnorm(-d2)
          
        }
        
      }
      
      return(greekval)
    }
    
    greekval <- bsgreek(S,K,r,t,sigma,div,opttype,greek)
    
    plot(1:20, 1:20, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(10, 10, paste(greek,"= ", round(greekval,3), sep=""),cex=1.5)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "Option Greeks")
  rp.textentry(panel=my.panel,variable=S,title="Spot:            ",action=my.redraw,initval=100)
  rp.textentry(panel=my.panel,variable=K,title="Strike:          ",action=my.redraw,initval=110)
  #rp.textentry(panel=my.panel,variable=r,title="Risk free:     ",action=my.redraw,initval=0.05)
  rp.textentry(panel=my.panel,variable=t,title="Maturity:     ",action=my.redraw,initval=0.5)
  #  rp.textentry(panel=my.panel,variable=sigma,title="Sigma:         ",action=my.redraw,initval=0.30)
  rp.textentry(panel=my.panel,variable=div,title="Div yield:     ",action=my.redraw,initval=0.0)
  
  rp.radiogroup(panel = my.panel, variable= opttype,
                vals = c("Call", "Put"), 
                action = my.redraw, title = "Type of Option")
  rp.radiogroup(panel = my.panel, variable= greek,
                vals = c("Delta", "Gamma","Vega","Theta","Rho"), 
                action = my.redraw, title = "Greek")
  rp.doublebutton(panel = my.panel, showvalue=TRUE,variable= sigma, 
                  initval=0.25,step = 0.05, range = c(0.00, 0.50),
                  title = "sigma",  action = my.redraw)
  rp.doublebutton(panel = my.panel, showvalue=TRUE,variable= r, 
                  initval=0.04,step = 0.01, range = c(0.00, 0.10),
                  title = "risk free", action = my.redraw)
  rp.tkrplot(my.panel, my.tkrplot, my.draw) # doesnt appear very good
  #rp.do(my.panel, my.draw)
}
