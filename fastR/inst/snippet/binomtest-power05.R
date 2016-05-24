p <- seq(0,1,by=0.02);
power <- 1 - ( pbinom(60,100,p) - pbinom(39,100,p) );
myplot <- xyplot(power~p,ylab="power",xlab=expression(pi[a]),
            type='l', lwd=2);

