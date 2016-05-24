`sensitivity` <-
function(x, numberf, order=4,make.plot=FALSE, show.legend=TRUE,
plot.max=max(ff[-1]), include.total.variance=FALSE, cukier=TRUE,
names=paste(sep="", "P", 1:numberf), main="", xlab="frequency",
ylab="Fourier Coef", pch=rep(0,numberf), col=(1:numberf)+1, reorder=1:numberf, ...){
if(cukier){
        t.runs= min_runs_cukier75[numberf]
        t.freq = freq_cukier(numberf)[reorder]
} else {
        t.runs <- min_runs_mcrae82[numberf]
        t.freq <- freq_mcrae82(numberf)[reorder]
}
if(NROW(x)<t.runs){
cat("x is too short. Expected number of values: ", t.runs,
    "\n but found ", NROW(x), " values")
    return
}
y<- na2mean(x)
ff <- abs(fft(double_serie(y)))/NROW(x)
#Debug-Message 
#cat("Frequencies: ", t.freq, "Using data up to ", (NROW(x)+1), "\n")
ff <- ff[1:(NROW(x)+1)]

#Frequenzen und deren vielfache 
# McRae garantiert unabh. bis Ordnung 4
freq <- t.freq %o% 1:order

if(make.plot){
    freqcol= rep(col,order)
    freqpch= rep(pch,order)
    plot(c(1,NROW(ff))-1,c(0,plot.max),
              t="n", xlab=xlab,
              ylab=ylab, main=main, ...)
    points(((1:NROW(ff))-1)[-as.vector(freq)],ff[-as.vector(1+freq)])
    points(as.vector(freq), ff[as.vector(1+freq)], col=freqcol, pch=freqpch)
    if(show.legend){
        legend("topright", inset=0.1, names, pch=pch, col=col)
    }
}

#Sensitivitaeten berechnen
sens<-function(frequency){
sum((ff[drop(frequency)])^2, na.rm=TRUE)
}
# Sensitivitaeten fuer Frequenzen berechnen
# Zu beachten, dass fft[1] dem A0 entspricht
# und nicht der Frequenz 1. Fuer Total auch A0 wegnehmen?

total<-sum(ff[-1]^2)
#browser()
#total<-sum(ff^2)
#total<-1
toReturn <- apply(1 +freq, 1, sens)/total
if(include.total.variance) toReturn <- c(toReturn ,total)
return(toReturn)

}

