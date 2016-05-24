if (getRversion() >= "2.15.1") utils::globalVariables(c('weight1',
                                                        'mu1',
                                                        'sigma1',
                                                        'mu2',
                                                        'sigma2',
                                                        'plottype'))


varbehavior <-
function(){
  
  
   my.draw1<- function(panel){
    
     with(panel, {
    
    weight1=as.numeric(panel$weight1)
    
    mu1=as.numeric(panel$mu1)
    sigma1=as.numeric(panel$sigma1)
    mu2=as.numeric(panel$mu2)
    sigma2=as.numeric(panel$sigma2)
    
    
    VaR1 <- function(conflvl, correl){
      
      horizon1=12
      
      weights = matrix(c(weight1,1-weight1),nrow=2)
      portreturn = mu1*weights[1,1] + mu2*weights[2,1]
     cov12 = correl*sigma1*sigma2
#       cov = matrix(c(sigma1^2,cov12,cov12,sigma2^2),nrow=2,byrow=F)
#     portrisk = t(weights) %*% (cov %*% weights)
      portrisk = weights[1,1]^2*sigma1^2 + weights[2,1]^2*sigma2^2 + 2*cov12*weights[1,1]*weights[2,1]
      portvar= (qnorm(conflvl)*portrisk*sqrt(horizon1/12) - portreturn)
      return(portvar)
    }
    
    
    VaR2 <- function(conflvl, horizon1){
      
      correl = 0.5
      
      weights = matrix(c(weight1,1-weight1),nrow=2)
      portreturn = mu1*weights[1,1] + mu2*weights[2,1]
      cov12 = correl*sigma1*sigma2
#       cov = matrix(c(sigma1^2,cov12,cov12,sigma2^2),nrow=2,byrow=F)
#       portrisk = t(weights) %*% (cov %*% weights)
      portrisk = weights[1,1]^2*sigma1^2 + weights[2,1]^2*sigma2^2 + 2*cov12*weights[1,1]*weights[2,1]
      portvar= (qnorm(conflvl)*portrisk*sqrt(horizon1/12) - portreturn)
      return(portvar)
    }
    
    
    
    if (panel$plottype=="Confidence level-Correlation"){
      
      conflvl <- seq(0.9, 0.999, length=31) 
      correl <- seq(-1, 1, length=31)
      x.axis.variable = conflvl
      x.axis.variable.name <- "Confidence level"
      y.axis.variable = correl
      y.axis.variable.name <- "Correlation"
      
      var <- outer(conflvl, correl, FUN = VaR1)
      
      if (length(dev.list()) == 0) 
        dev.new()
      colors <- c("cyan", "steelblue", "green","greenyellow" , "lightgreen","deepskyblue" ,"darksalmon","gold", "skyblue", "orange", "violet")
      my.title <- paste("Value at Risk")
      persp(x.axis.variable, y.axis.variable, var, xlab = x.axis.variable.name, ylab = y.axis.variable.name,ticktype="detailed",main= my.title,theta= -40,phi= 10,col="darksalmon")
      panel
    }
    
    else if (panel$plottype=="Confidence level-Horizon"){
      
      conflvl <- seq(0.9, 0.999, length=31) 
      horizon1 <- seq(1,12,length=31)
      x.axis.variable = conflvl
      x.axis.variable.name <- "Confidence level"
      y.axis.variable = horizon1
      y.axis.variable.name <- "Horizon"
      
      var <- outer(conflvl, horizon1, FUN = VaR2)
      
      if (length(dev.list()) == 0) 
        dev.new()
      colors <- c("cyan", "steelblue", "green","greenyellow" , "lightgreen","deepskyblue" ,"darksalmon","gold", "skyblue", "orange", "violet")
      my.title <- paste("Value at Risk")
      persp(x.axis.variable, y.axis.variable, var, xlab = x.axis.variable.name, ylab = y.axis.variable.name,ticktype="detailed",main= my.title,theta= -40,phi= 10,col="gold")
      panel
    }
    
    panel
     })
    
  }
  
  
  my.panel <- rp.control(title = "Value at Risk- Two stocks")

#    change below to doublebuttons/sliders


 rp.doublebutton(panel=my.panel,variable=weight1,title="Weight 1",step=0.05,showvalue=T,
                   action=my.draw1,range=c(0,1), initval=0.5)
 rp.doublebutton(panel=my.panel,variable=mu1,title="Mu 1",step=0.05,showvalue=T,
                 action=my.draw1,range=c(0,0.2), initval=0.1)
 rp.doublebutton(panel=my.panel,variable=sigma1,title="Sigma 1",step=0.05,showvalue=T,
                 action=my.draw1,range=c(0,0.5), initval=0.2)
 rp.doublebutton(panel=my.panel,variable=mu2,title="Mu 2",step=0.05,showvalue=T,
                 action=my.draw1,range=c(0,0.2), initval=0.1)
 rp.doublebutton(panel=my.panel,variable=sigma2,title="Sigma 2",step=0.05,showvalue=T,
        action=my.draw1,range=c(0,0.5), initval=0.2)  
   
 
  rp.radiogroup(panel=my.panel,variable=plottype,title="Plot Type",
      vals=c("Confidence level-Correlation","Confidence level-Horizon"),action=my.draw1)
  rp.do(my.panel, my.draw1)
}
