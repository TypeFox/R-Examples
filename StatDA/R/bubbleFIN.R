######################################################################
# Bubbleplot due to Finnish method:
######################################################################
circles <- function(x, y, radius, col=NA, border=par("fg")) {
#draw circles
    nmax <- max(length(x), length(y));
    if (length(x) < nmax) x <- rep(x, length=nmax);
    if (length(y) < nmax) y <- rep(y, length=nmax);
    if (length(col) < nmax) col <- rep(col, length=nmax);
    if (length(border) < nmax) border <- rep(border, length=nmax);
    if (length(radius) < nmax) radius <- rep(radius, length=nmax);
    theta <- 2* pi * seq(0, 355, by=5) / 360;
    ct <- cos(theta);
    st <- sin(theta);
    for(i in 1:nmax)
	polygon(x[i] + ct * radius[i], y[i] + st * radius[i], col=col[i], border=border[i]);
}

######################################################################
bubbleFIN <-function(x,y,z,radi=10000,S=9,s=0.9,wa=0,wb=0.95,wc=0.05, plottitle="BubblePlot",legendtitle="Legend", text.cex=1,
        legtitle.cex=1,backgr="kola.background",leg=TRUE,ndigits=1)
{
# ndigits ... number of digits for the legend
#
	data<-data.frame(x=x,y=y,z=z);

#Start bubbling
#maximum radius =~ 1000m
   if (min(z)<0) {
        zerofactor<-abs(min(z));
    }else {
        zerofactor<-0;
    }
    zz <- z+zerofactor
    q1 <- quantile(zz,0.1)
    q2 <- quantile(zz,0.99)
    C <- q2/((q1/q2)^(wc/wb))
    c <- q1/((q2/q1)^(wa/wb))
    xi <- pmax(pmin(zz,C),c)
    di <- s*(S/s)^(log10(xi/c)/log10(C/c))

    dataf <- radi/S
    ord<-order(z,decreasing=TRUE);
    circles(x[ord],y[ord],dataf*di[ord],col=1,border=9)
#Add legend, top right
if (leg==TRUE){
    dataq<-quantile(di,probs=seq(1,0,-0.1))
    zq<-quantile(z,probs=seq(1,0,-0.1))
    xc<-max(x)-3*radi
    yc<-max(y)-radi
#draw legend bubbles, largest diameter first
    for (i in 1:length(dataq)) {
	diameter<- dataq[i]*dataf
	yc<-yc-1.3*radi
	circles(xc,yc,diameter, col=1,border=9)
#	text(xc+dataq[1]*dataf+2*radi,yc,format(round(zq[i],digits=2),nsmall=1),cex=text.cex)
	text(xc+dataq[1]*dataf+2*radi,yc,roundpretty(zq[i],maxdig=ndigits),cex=text.cex)
    } 
    text(xc+radi,max(y),legendtitle,cex=legtitle.cex)
}
#end bubble function
}

