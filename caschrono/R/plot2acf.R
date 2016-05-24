plot2acf= function(y1, y2, lag.max=40, main=c("",""))
{

acfa = stats::acf(y1, lag.max = lag.max, plot=FALSE)$acf[-1]
acfb = stats::acf(y2, lag.max = lag.max, plot=FALSE)$acf[-1]

titre1 = main[1]; titre2=main[2]

# une meme echelle pour les 2 graphiques
# 2 acf calcules au mm retard
maxu = min(1, 1.1*max(c(abs(acfa),abs(acfb))))
maxlag= length(acfa)
LAG = 1:maxlag
# en haut 
op <- par(no.readonly = TRUE)
par(fig= c(0,1,0.5,1), cex=.9)
par(mai=c(0,0.810,0.5,0.324))  # mai remplace mar

# par()$fin
# ACF   Haut   et titre ##############################################################
plot(LAG, acfa, type="h",ylim=c(-maxu,maxu),xlab="",xaxt="n",ylab="",
     cex=.5,cex.lab=.6,cex.axis=.8,las=1)
abline(h=c(-1.96/sqrt(length(y1)),1.96/sqrt(length(y1))),lty=2,col='blue')
abline(h=0)
# text(.8*maxlag,0.9*maxu, labels=titre1, pos=3, cex=.9,font=3)   
legend("topright", legend=titre1, cex=.9, bty="n")  
# ACF   bas    ##############################################################   
par(new=TRUE)
par(fig= c(0,1,0,0.5)) #, cex=.8, mex=.6,cex.lab=.85)
par(mai=c(0.5,0.810,0,0.324))                                                                                        
plot(LAG, acfb, type="h",ylim=c(-maxu,maxu),xlab="",xaxt="n",ylab="",
     cex=.5,cex.lab=.6,cex.axis=.8,las=1)
abline(h=c(-1.96/sqrt(length(y2)),1.96/sqrt(length(y2))),lty=2,col='blue')     
axis(1, at = c(1,12,24,36), labels = c("1","12","24","36"),  tick = TRUE, line = NA,
     pos = NA, outer = FALSE, font = NA, lty = "solid",
      lwd = 1, lwd.ticks = 0.5, col = NULL, col.ticks = NULL,
      hadj = NA, padj = NA)
abline(h=0)
legend("topright", legend=titre2, cex=.9, bty="n")  
par(op) 
}
