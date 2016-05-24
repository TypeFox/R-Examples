if (getRversion() >= "2.15.1") utils::globalVariables(c('val1',
                                                        'val2',
                                                        'mu1',
                                                        'mu2',
                                                        'sigma1',
                                                        'sigma2',
                                                        'conf',
                                                        'corr',
                                                        'horizon'))

var2stocks <-
function(){
  
  my.draw <- function(panel) {
    val1=as.numeric(panel$val1)
    val2=as.numeric(panel$val2)
    mu1=as.numeric(panel$mu1)
    sigma1=as.numeric(panel$sigma1)
    mu2=as.numeric(panel$mu2)
    sigma2=as.numeric(panel$sigma2)
    corr=as.numeric(panel$corr)
    conf = as.numeric(panel$conf)
    horizon <- as.numeric(panel$horizon)/12
    
    portsize = sum(val1+val2)
    weights = matrix(c(val1,val2)/portsize,nrow=2)
    portreturn = mu1*weights[1,1] + mu2*weights[2,1]
    cov12 = corr*sigma1*sigma2
    cov = matrix(c(sigma1^2,cov12,cov12,sigma2^2),nrow=2,byrow=F)
    portrisk = t(weights) %*% (cov %*% weights)
    portvar = (qnorm(conf)*portrisk*sqrt(horizon) - portreturn)*portsize
    
    portvar<-round(portvar,2)
    plot(1:20, 1:20, type="n", xlab="", ylab="",
         axes=FALSE, frame = TRUE)
    text(10, 10, paste("VaR= ", portvar, sep=""),cex=1.5)
    
    panel
  }

  
  my.redraw <- function(panel) {
    rp.tkrreplot(panel, my.tkrplot)
    panel
  }
  
  
      
  my.panel <- rp.control(title = "Value at Risk- Two stocks")
  rp.textentry(panel=my.panel,variable=val1,title="Value1:          ",action=my.redraw,initval=1000)
  rp.textentry(panel=my.panel,variable=val2,title="Value2:          ",action=my.redraw,initval=2000)
  rp.textentry(panel=my.panel,variable=mu1,title="mu1:             ",action=my.redraw,initval=0.10)
  rp.textentry(panel=my.panel,variable=mu2,title="mu2:             ",action=my.redraw,initval=0.08)
  rp.textentry(panel=my.panel,variable=sigma1,title="sigma1:        ",action=my.redraw,initval=0.30)
  rp.textentry(panel=my.panel,variable=sigma2,title="sigma2:        ",action=my.redraw,initval=0.25)
  rp.textentry(panel=my.panel,variable=conf,title="Conf level:   ",action=my.redraw,initval=0.95)
  rp.textentry(panel=my.panel,variable=corr,title = "Correlation:",action=my.redraw,initval=0.5)
  rp.doublebutton(my.panel,variable=horizon,step=1,title="Horizon (Months)",initval=12,showvalue=T,range=c(1,12),action=my.redraw)
#   rp.slider(panel=my.panel,variable=horizon,from=0.25,to=1,resolution=0.1,title = "Horizon:",action=my.redraw,showvalue=T,initval=1)
  rp.tkrplot(panel = my.panel, name = my.tkrplot, plotfun = my.draw)
  #rp.do(my.panel, my.draw)
#   rp.radiogroup(panel=my.panel,variable=plottype,title="Plot Type",vals=c("Confidence level-Correlation","Confidence level-Horizon"),action=my.draw1)
#  rp.do(my.panel, my.draw1)
  }
