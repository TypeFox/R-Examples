Boot.PI <-
function(x,p,h,nboot,prob)
{
#set.seed(12345)
n <- nrow(x)
B <- OLS.AR(x,p,h,prob)
bb <- LSMB(x,p)$coef
eb <- sqrt( (n-p) / ( (n-p)-length(bb)))*LSMB(x,p)$resid
ef <- sqrt( (n-p) / ( (n-p)-length(bb)))*LSM(x,p)$resid

fore <- matrix(NA,nrow=nboot,ncol=h)
for(i in 1:nboot)
    {    
        index <- as.integer(runif(n-p, min=1, max=nrow(eb)))
        es <- eb[index,1]
        xs <- ysb(x, bb, es)
        bs <- LSM(xs,p)$coef 
        fore[i,] <- AR.ForeB(xs,bs,h,ef,length(bs)-1)
    }

Interval <- matrix(NA,nrow=h,ncol=length(prob),dimnames=list(1:h,prob))
for( i in 1:h)
Interval[i,] <- quantile(fore[,i],probs=prob)
return(list(PI=Interval,forecast=B$forecast))
}
