boot_autocov <-
function(series,replic=5000,plot=TRUE,alpha=0.05){

if (is.ts(series)==TRUE){


library(boot)
kas = tsboot(series,statistic,R=replic,sim="scramble")


quantiles  = matrix(0,length(kas$t[1,]),3)

for (i in 2:length(kas$t[1,])){
cp            = kas$t[,i]
quantiles[i,1]= quantile(cp,alpha) 
quantiles[i,2]= quantile(cp,1-alpha/2)
quantiles[i,3]= mean(cp)
}

quantiles  = quantiles[-1,]

if (plot==TRUE){

par(mfrow=c(1,2))

x=seq(1,length(quantiles[,1]),1)/frequency(series)
plot(x,quantiles[,1],type="l",col="blue",main="Bootstraped Correlogram",ylab="value",lwd=1,xlab="lag")
polygon(c(x,rev(x)),c(quantiles[,1],rev(quantiles[,2])),col="skyblue")
lines(x,quantiles[,3],type="o",col="black",pch=20)

abline(a=0,b=0)

plot(acf(series,plot=FALSE),main="Asymptotic Correlogram",ylim=c(-1,1))

}

lista = list(average=quantiles[,1],upper=quantiles[,2],lower=quantiles[,3])
return(lista)

}else{

return ('Object is not a time-series')

}

}
