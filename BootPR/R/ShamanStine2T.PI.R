ShamanStine2T.PI <-
function(x,p,h,nboot,prob,pmax)
{
#set.seed(12345)
n <- nrow(x)

BC <- Shaman.StineT(x,p,h)
BCB <- Shaman.StineBT(x,p,h)

bb <- BCB$coef
eb <- sqrt( (n-p) / ( (n-p)-p-1))*BCB$resid
ef <- sqrt( (n-p) / ( (n-p)-p-1))*BC$resid

fore <- matrix(NA,nrow=nboot,ncol=h)
for(i in 1:nboot)
    {    
        index <- as.integer(runif(n-p, min=1, max=nrow(eb)))
        es <- eb[index,1]
        xs <- ysbT(x, bb, es)
        ps <- ART.order(xs,pmax)$ARorder[1]
        bs <- Shaman.StineT(xs,ps,h)$coef
        fore[i,] <- ART.ForeB(xs,bs,h,ef,length(bb)-1)
    }

Interval <- matrix(NA,nrow=h,ncol=length(prob),dimnames=list(1:h,prob))
for( i in 1:h)
Interval[i,] <- quantile(fore[,i],probs=prob)
return(list(PI=Interval,forecast=BC$forecast))
}
