################################################################################
#  Plot of objects of class 'select.parfm'                                     #
################################################################################
#                                                                              #
#                                                                              #
#                                                                              #
#   Date: December 22, 2011                                                    #
#   Last modification on: October 17, 2012                                     #
################################################################################

plot.select.parfm <- function(x, 
                              mar=c(2.5, 2, 1.5, .5),
                              ...){
  par(mfrow=c(1, 3))
  
  ### --- AIC --- ###
  par(mar=mar)
  plot(0,0, ty="n", xlab="", ylab="", main="AIC", xaxt="n",
    ylim=c(min(x$AIC, na.rm=TRUE) *  .9975,
           max(x$AIC, na.rm=TRUE) * 1.0025),
    xlim=c(.5, ncol(x$AIC) + .5), cex.lab=1.5)
  abline(v=1:ncol(x$AIC), col="grey")    
  
  mtext(c(none="No",
          gamma="Ga",
          ingau="IG",
          possta="PS",
          lognor="LN")[colnames(x$AIC)],
        side=1, at=1:ncol(x$AIC), padj=1)
  
  for (i in 1:nrow(x$AIC)) points(
    (1:ncol(x$AIC)), x$AIC[i, ],
    pch=i, cex=1.5)
  
  
  ### --- names --- ###
  par(mar=mar)
  plot(0:2, 0:2, xaxt="n", yaxt="n", bty="n", xlab="", ylab="",
       ty="n")
  
  legend(c(.3, 1.7), c(1, 1.75),
         c(exponential="exponential", weibull="Weibull", 
           inweibull="inverse Weibull",
           gompertz="Gompertz", loglogistic="loglogistic", 
           lognormal="lognormal")[rownames(x$AIC)],
         pch=1:nrow(x$AIC), 
         bg="white", bty="n",
         ncol=1, cex=1.5, xjust=.5)
  
  legend(c(0, 2), c(.25, 1), yjust=1,
         mapply(paste, 
                c(none="No",
                  gamma="Ga",
                  ingau="IG",
                  possta="PS",
                  lognor="LN")[colnames(x$AIC)],
                c(none="no frailty",
                  gamma="gamma",
                  ingau="inverse Gaussian",
                  possta="positive stable",
                  lognor="lognormal")[colnames(x$AIC)],
                sep=" = "),
         bg="white", bty="n",
         ncol=1, cex=1.5, xjust=.5)
  ### --- end names --- ###
  
  
 
  ### --- BIC --- ###
  par(mar=mar)
  plot(0,0, ty="n", xlab="", ylab="", main="BIC", xaxt="n",
    ylim=c(min(x$BIC, na.rm=TRUE) *  .9975,
           max(x$BIC, na.rm=TRUE) * 1.0025),
    xlim=c(.5, ncol(x$BIC) + .5), cex.lab=1.5)
  abline(v=1:ncol(x$BIC), col="grey")    
  
  mtext(c(none="No",
          gamma="Ga",
          ingau="IG",
          possta="PS",
          lognor="LN")[colnames(x$BIC)],
        side=1, at=1:ncol(x$BIC), padj=1)
  
  for (i in 1:nrow(x$BIC)) points(
    (1:ncol(x$BIC)), x$BIC[i, ],
    pch=i, cex=1.5)
}

