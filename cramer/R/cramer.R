#require(boot)

# ".cramer.statistic" calculates the value of the test statistic T
# it is build so it can be called by "boot"
#   daten   is the original data set (is just given since boot wants it like that))
#   indexe  is the bootstrapped vector of indices (first m are observations of X, the next n of Y)
#   m       is the number of observations of X
#   n       is the number of observations of Y
#   lookup  is a matrix containing the distances of points i und j
.cramer.statistic<-function(daten,indexe,mm,nn,lookup) {
    xind<-indexe[1:mm]
    yind<-indexe[(mm+1):(mm+nn)]
    mm*nn/(mm+nn)*(2*sum(lookup[xind,yind])/(mm*nn)-sum(lookup[xind,xind])/(mm^2)-sum(lookup[yind,yind])/(nn^2))
}

.cramer.characteristicfunction<-function(lambdasquare,t) { z<--0.5*log(1-2i*lambdasquare%*%t(t));return(exp(complex(length(t),rowsum(Re(z),rep(1,length(lambdasquare))),rowsum(Im(z),rep(1,length(lambdasquare)))))) }

.cramer.kritwertfft<-function(lambdasquare,conf.level,maxM,K) {
	M<-2^11
	while (150*pi*M/K^2<(2*sum(lambdasquare)+lambdasquare[1])) M<-M*2
	M<-min(c(M,maxM))
	goodlimit<-150*pi*M/K^2
	a<-0                    
	t<-0:(M-1)*K/M                    
	x<-0:(M-1)*2*pi/K
	#Verteilungsfunktion
	t[1]<-1                    
	h<-.cramer.characteristicfunction(lambdasquare,t)/t*exp(-a*1i*t);
	h[1]<-complex(1,0,1)*sum(lambdasquare)
	Fx<-1/2-Im(K/(M*pi)*fft(h,inverse=FALSE))+K/(2*M*pi)*(sum(lambdasquare)+x+a)
	#quantile
        xindex<-which.min(abs(Fx[1:(M/2)]-conf.level))       
        #linear interpolation between two evaluated points of the distribution function (much better than just taking x[xindex]
        if (Fx[xindex]>conf.level) xindex<-xindex-1
        if (xindex<1) xindex<-1
        quantile<-x[xindex]+(conf.level-Fx[xindex])*(x[xindex+1]-x[xindex])/(Fx[xindex+1]-Fx[xindex])
	#warnings?
	if (Fx[M/2]<conf.level) warning("Quantile calculation discrepance. Try to increase K!")
	if (quantile>goodlimit) warning("Quantile beyond good approximation limit. Try to increase maxM or decrease K!")
	kritwert <- list(quantile=quantile,hypdist.x=x,hypdist.Fx=Fx)
	return(kritwert)
}	

cramer.test<-function(x,y,conf.level=0.95,replicates=1000,sim="ordinary",just.statistic=FALSE,kernel="phiCramer",maxM=2^14,K=160) {    
    RVAL <- list(method = paste("nonparametric Cramer-Test with kernel",kernel,"\n(on equality of two distributions)"),
                 d = 0,
                 m = 0,
                 n = 0,
                 statistic = 0,
                 conf.level = conf.level,
                 crit.value = 0,
                 p.value = 0,                 
                 result = 0,
                 sim = sim,
                 replicates = replicates,
		 hypdist = 0,
		 ev = 0) 
    if ((is.vector(x))&&(is.vector(y))) RVAL$d<-1
    if ((is.matrix(x))&&(is.matrix(y))) if (ncol(x)==ncol(y)) RVAL$d<-ncol(x)
    if (RVAL$d==0) stop("types of x and y incompatible or inappropriate.")    
    if (RVAL$d==1) {
        RVAL$m<-length(x)
        RVAL$n<-length(y)
        daten<-matrix(c(x,y),ncol=1,byrow=TRUE)
    } else {
        RVAL$m<-nrow(x)
        RVAL$n<-nrow(y)
        daten<-matrix(c(t(x),t(y)),ncol=ncol(x),byrow=TRUE)
    }
    lookup<-matrix(rep(0,(RVAL$m+RVAL$n)^2),ncol=(RVAL$m+RVAL$n))
    for (i in 2:(RVAL$m+RVAL$n))
         for (j in 1:(i-1)) {
             lookup[i,j]<-sum((daten[i,]-daten[j,])^2)
             lookup[j,i]<-lookup[i,j]
         }
    lookup<-eval(call(kernel,lookup))
#print(date()) #!!!
    if (just.statistic) {
        RVAL$statistic<-.cramer.statistic(daten,1:(RVAL$m+RVAL$n),RVAL$m,RVAL$n,lookup)
    } else if (sim!="eigenvalue") {
        b<-boot(data=daten,statistic=.cramer.statistic,mm=RVAL$m,nn=RVAL$n,lookup=lookup,sim=RVAL$sim,stype="i",R=RVAL$replicates)
        RVAL$statistic<-b$t0
        RVAL$p.value<-1-rank(c(b$t0,b$t))[1]/(replicates+1)
        RVAL$crit.value<-sort(b$t)[round(RVAL$conf.level*RVAL$replicates)]
        if (RVAL$statistic>RVAL$crit.value) RVAL$result<-1
	RVAL$hypdist.x<-sort(b$t)
	RVAL$hypdist.Fx<-(1:RVAL$replicates)/RVAL$replicates
    } else {
	RVAL$statistic<-.cramer.statistic(daten,1:(RVAL$m+RVAL$n),RVAL$m,RVAL$n,lookup)
	N<-RVAL$m+RVAL$n
	# alternative Berechnung von B (DEUTLICH schneller)
	C1<-rep(0,N)
	for (i in 1:N)     
	    for (j in 1:N) C1[i]<-C1[i]+lookup[i,j] 
	C1<-C1/N
	C2<-0
	for (i in 1:N)     
	    for (j in 1:N) C2<-C2+lookup[i,j] 
	C2<-C2/N^2
	
	B<-matrix(rep(0,N^2),ncol=N)
	for (i in 1:N) 
	    for (j in 1:N) B[i,j]<-C1[i]+C1[j]-C2-lookup[i,j] 
	B<-B/N

	RVAL$ev<-eigen(B,FALSE) # <--only.values   TRUE)
	kritwert<-.cramer.kritwertfft(Re((RVAL$ev)$values),conf.level,maxM,K)

	RVAL$p.value<-1-(kritwert$hypdist.Fx)[which.min(abs((kritwert$hypdist.x)[1:(3*length(kritwert$hypdist.x)/4)]-RVAL$statistic))];
	RVAL$crit.value<-kritwert$quantile
	RVAL$hypdist.x<-kritwert$hypdist.x
	RVAL$hypdist.Fx<-kritwert$hypdist.Fx
	if (RVAL$statistic>RVAL$crit.value) RVAL$result<-1
    }
    class(RVAL) <- "cramertest"
#print(date()) #!!!
    return(RVAL)
}

print.cramertest<-function(x,...) {
    cat("\n",x$d,"-dimensional ",x$method,"\n\n")
    cat("\tx-sample: ",x$m," values        ")
    cat("y-sample: ",x$n," values\n\n")
    if (x$crit.value>0) {
        cat("critical value for confidence level ",format(100 * x$conf.level),"% : ",x$crit.value,"\n")
        cat("observed statistic ",x$statistic,", so that\n\t hypothesis (\"x is distributed as y\") is ")
        cat(ifelse(x$result==0," ACCEPTED"," REJECTED"),".\n")
        cat("estimated p-value = ",x$p.value,"\n\n")
	if (x$sim!="eigenvalue") {
          cat("\t[result based on ",x$replicates," ",x$sim," bootstrap-replicates]\n\n");
	} else {
          cat("\t[result based on eigenvalue decomposition and inverse fft]\n\n");
	}
    } else {
        cat("observed statistic ",x$statistic,"\n\n")
    }
    invisible(x)
}


phiCramer<-function(x) return(sqrt(x)/2)
phiBahr<-function(x) return(1-exp(-x/2))
phiLog<-function(x) return(log(1+x))
phiFracA<-function(x) return(1-1/(1+x))
phiFracB<-function(x) return(1-1/(1+x)^2)
