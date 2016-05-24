Auto.Q <-
function(y,lags=10)
{
data <- y - mean(y)
T <- length(data)
ac <- matrix(acf(data,lag.max=lags,plot=F)$acf[, , 1])
ac <- ac[2:(lags+1),1]


ac1 <- matrix(acf(data,lag.max=lags,type="covariance",plot=F)$acf[, , 1])
ac2 <- matrix(NA,nrow=lags)
    for( i in 1:lags){
    y <- data[(i+1):T]^2
    x <- data[1:(T-i)]^2
    t <- length(y)
    ac2[i] <- crossprod(x,y)/t
    }
    
ac3 <- (ac1[2:(lags+1),1]^2) /ac2
#Qps <- T*sum(ac3)
BP<- T * cumsum(ac3)

aux <- matrix(1:lags)
q <- 2.4
maxro <- sqrt(T)*max(sqrt(ac3))
pin <- 2
if(maxro <= sqrt(q*log(T)) ) pin <-log(T)

Lp <- BP-aux*pin;
phat <- which.max(Lp)
Tn <- BP[phat]
pvalue <- 1-pchisq(Tn,1)
return(list(Stat=Tn,Pvalue=pvalue))
}
