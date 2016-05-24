if (getRversion() >= "2.15.1") utils::globalVariables(c("r","pd","defaulttimes","recovery"))

cdswap <-
function(){
  my.draw <- function(panel) {
    
    #S <-as.numeric(panel$S)
    r <-as.numeric(panel$r)
    t <-as.numeric(panel$t)
    pd <-as.numeric(panel$pd)
    recovery<-as.numeric(panel$recovery)
    
    discrates=rep(0,t)
    allyears<-1:t
    discrates<- exp(-r*allyears)
    PD<- rep(0,t)
    PD[1]<- pd
    SP<- rep(0,t)
    SP[1]<-1-pd
    
    for (i in 2:t){
      PD[i]<-PD[i-1]*(1-pd)
      SP[i]<-SP[i-1]*(1-pd)
    }
    
    if (panel$defaulttimes=="End of Q1"){
      spreadfactor <- 0.25
      times<-seq(from=0.25,by=1,length.out=t)
      pvfactors<-exp(-r*times)  
    }
    else if(panel$defaulttimes=="End of half year"){
      spreadfactor <- 0.5
      times<-seq(from=0.5,by=1,length.out=t)
      pvfactors<-exp(-r*times)  
    }
    else if(panel$defaulttimes=="End of Q3"){
      spreadfactor <- 0.75
      times<-seq(from=0.75,by=1,length.out=t)
      pvfactors<-exp(-r*times)  
    }
    else{
      spreadfactor <- 1
      times<-seq(from=1,by=1,length.out=t)
      pvfactors<-exp(-r*times)  
    }
    
    PVpayoffs<- sum(PD*pvfactors)*(1-recovery)
    PVaccruals<-sum(PD*pvfactors)*spreadfactor
    PVpayments<-sum(SP*discrates)
    
    price=PVpayoffs/(PVaccruals+PVpayments)
    price= round(price*10000,0)
    
    plot(1:20, 1:20, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(10, 10, paste("Spread:", price, " bp"),cex=1.4)
    
    panel
  }
  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  my.panel <- rp.control(title = "CDS Spread")
  #rp.textentry(panel=my.panel,variable=S,title="Spot:               ",action=my.redraw,initval=100)
  rp.textentry(panel=my.panel,variable=r,title="Risk free:        ",action=my.redraw,initval=0.05)
  rp.textentry(panel=my.panel,variable=t,title="Maturity:        ",action=my.redraw,initval=5)
  rp.textentry(panel=my.panel,variable=pd,title="Prob. Default:",action=my.redraw,initval=0.02)
  rp.radiogroup(panel=my.panel,variable=defaulttimes, title="Default Assumption", 
                vals=c("End of Q1", "End of half year","End of Q3", "End of Year"),action=my.redraw)
  rp.doublebutton(panel=my.panel,variable=recovery,step=0.05,range=c(0,1),initval=0.4,action=my.redraw,showvalue=T)
  #   rp.slider(panel=my.panel,variable=recovery,from=0,to=1,initval=0.4,
  #             resolution=.05,showvalue=TRUE,action=my.redraw,title="Recovery rate")
  rp.tkrplot(panel = my.panel, pos="right", name = my.tkrplot, plotfun = my.draw)
  
  #rp.do(my.panel, my.draw)
  
  
}
