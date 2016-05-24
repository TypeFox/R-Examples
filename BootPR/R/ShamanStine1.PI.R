ShamanStine1.PI <-
function(x,p,h,nboot,prob)
{
#set.seed(12345)
n <- nrow(x)

BC <- Shaman.Stine(x,p,h)
BCB <- Shaman.StineB(x,p,h)

bb <- BCB$coef
eb <- sqrt( (n-p) / ( (n-p)-length(bb)))*BCB$resid
ef <- sqrt( (n-p) / ( (n-p)-length(bb)))*BC$resid

fore <- matrix(NA,nrow=nboot,ncol=h)
for(i in 1:nboot)
    {    
        index <- as.integer(runif(n-p, min=1, max=nrow(eb)))
        es <- eb[index,1]
        xs <- ysb(x, bb, es)
        bs <- Shaman.Stine(xs,p,h)$coef
        fore[i,] <- AR.ForeB(xs,bs,h,ef,length(bs)-1)
    }

Interval <- matrix(NA,nrow=h,ncol=length(prob),dimnames=list(1:h,prob))
for( i in 1:h)
Interval[i,] <- quantile(fore[,i],probs=prob)
return(list(PI=Interval,forecast=BC$forecast))
}
