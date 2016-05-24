p <- rep(c(0.52,0.55,0.60), each=2000);
plab <- paste("alt prob =", as.character(p));
n <- rep(1:2000,times=3);
critical <- qbinom(0.025,size=n,prob=p);
power <- 1 - ( pbinom(n-critical+1,n,p) - pbinom(critical-1,n,p) );
myplot <- xyplot(power~n|plab,ylab="power",xlab="number of coin tosses",
            ylim=c(0,1.1), type='l', lwd=2);
