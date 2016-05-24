plotacfthemp= function(y, ar = numeric(0), ma = numeric(0), lag.max = 20, titre="")
{

# computes theoretical ACF and PACF
acf.th =ARMAacf(ar=ar, ma=ma, lag.max=lag.max, pacf=FALSE)                         
pacf.th =ARMAacf(ar=ar, ma=ma, lag.max=lag.max, pacf=TRUE)   
                 
# computes empirical ACF and PACF                                                                                   
pacf.emp=stats::acf(y,lag.max=lag.max,type="partial",plot=FALSE)                                                                                                                   
acf.emp= stats::acf(y,lag.max=lag.max,type="correlation",plot=FALSE)                      

n = lag.max
ly = length(y)

# Initialisation
LAG = 1:n
pacfemp= list(y= pacf.emp$acf[,1,1], bsup=(1.96/sqrt(ly))*rep(1,n))
acfemp=  list(y= acf.emp$acf[-1,1,1], bsup=(1.96/sqrt(ly))*rep(1,n))
pacfth=  list(y= pacf.th, bsup= rep(0,n))
acfth=   list(y= acf.th[-1], bsup= rep(0,n))



# une même échelle pour les 4 graphiques
maxu = min(1, 1.1*max(c(abs(pacfemp$y),abs(acfemp$y), abs(pacfth$y), abs(acfth$y))))
maxlag= max(LAG)
L = -acfemp$bsup[1]; U=-L
# en haut à gauche
op <- par(no.readonly = TRUE)
par(fig=c(0,0.5,0.5,1))
par(mai=c(0,0.5,0.5,0))
par()$fin
# ACF TH Haut gauche et titre ##############################################################
plot(LAG, acfth$y, type="h",ylim=c(-maxu,maxu),xlab="",xaxt="n",ylab="",
     cex=.5,cex.lab=.6,cex.axis=.8,las=1)
abline(h=0)
text(.8*maxlag,0.9*maxu, labels=titre, pos=3, cex=.9,font=3)
text(.8*maxlag,0.7*maxu, labels="ACF th.", pos=3, cex=.9)
par(new=TRUE)
par(fig=c(0.5,1,0.5,1))
par(mai=c(0,0,0.5,0))
par()$fin
# PACF th haut droit ##############################################################
 plot(LAG, pacfth$y, type="h",ylim=c(-maxu,maxu),xlab="",xaxt="n",yaxt="n",
    cex=.5,cex.lab=.6,cex.axis=.6,las=1)
  abline(h=0)
  text(.8*maxlag,0.7*maxu, labels="PACF th.", pos=3, cex=.9)
par(new=TRUE)
par(fig=c(0,0.5,0,0.5))
par(mai=c(0.5,0.5,0,0))
par()$fin
# ACF emp bas gauche   ##############################################################
plot(LAG, acfemp$y, type="h",ylim=c(-maxu,maxu),xlab="",xaxt="n",ylab="",
     cex=.5,cex.lab=.6,cex.axis=.8,las=1)
axis(1, at = c(1,5,10,15,20), labels = c("1","5","10","15","20"),  tick = TRUE, line = NA,
     pos = NA, outer = FALSE, font = NA, lty = "solid",
     lwd = 1, lwd.ticks = 0.5, col = NULL, col.ticks = NULL,
     hadj = NA, padj = NA)
abline(h=0)
abline(h=L,lty=2,col="blue")
abline(h=U,lty=2,col="blue")
text(.8*maxlag,0.7*maxu, labels="ACF emp.", pos=3, cex=.9)
# PACF emp bas droit ##############################################################
par(new=TRUE)
par(fig=c(0.5,1,0,0.5))
par(mai=c(0.5,0,0,0))
par()$fin
plot(LAG, pacfemp$y, type="h",ylim=c(-maxu,maxu),xlab="",xaxt="n",ylab="",yaxt="n",
     cex=.5,cex.lab=.6,cex.axis=.6,las=1)
axis(1, at = c(1,5,10,15,20), labels = c("1","5","10","15","20"),  tick = TRUE, line = NA,
     pos = NA, outer = FALSE, font = NA, lty = "solid",
     lwd = 1, lwd.ticks = 0.5, col = NULL, col.ticks = NULL,
     hadj = NA, padj = NA)
abline(h=0)
abline(h=L,lty=2,col="blue")
abline(h=U,lty=2,col="blue")
text(.8*maxlag,0.7*maxu, labels="PACF emp.", pos=3, cex=.9)
par(op)
}
