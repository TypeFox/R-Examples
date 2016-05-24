BootAfterBootT.PI <-
function(x,p,h,nboot,prob)
{
#set.seed(12345)
n <- nrow(x)

B <- OLS.ART(x,p,h,prob)
BBC <- BootstrapT(x,p,h,200)
BBCB <- BootstrapTB(x,p,h,200)

bb <- BBCB$coef
eb <- sqrt( (n-p) / ( (n-p)-length(bb)))*BBCB$resid
bias <- B$coef - BBC$coef
ef <- sqrt( (n-p) / ( (n-p)-length(bb)))*BBC$resid

fore <- matrix(NA,nrow=nboot,ncol=h)
for(i in 1:nboot)
    {    
        index <- as.integer(runif(n-p, min=1, max=nrow(eb)))
        es <- eb[index,1]
        xs <- ysbT(x, bb, es)
        bs <- LSMT(xs,p)$coef 
        bsc <- bs-bias
        bsc <- adjust(bs,bsc,p)

        if(sum(bsc) != sum(bs))
        bsc[(p+1):(p+2),] <- RE.LSMT(xs,p,bsc)
    
        fore[i,] <- ART.ForeB(xs,bsc,h,ef,length(bs)-2)
    }

Interval <- matrix(NA,nrow=h,ncol=length(prob),dimnames=list(1:h,prob))
for( i in 1:h)
Interval[i,] <- quantile(fore[,i],probs=prob)
return(list(PI=Interval,forecast=BBC$forecast))
}
