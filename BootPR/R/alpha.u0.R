alpha.u0 <-
function(b,p,n)
{
alphavec <- c(-0.99,seq(-0.9,1,0.1))
stat <- matrix(NA,nrow=length(alphavec),ncol=2)

for( i in 1:length(alphavec) )
{alpha <- alphavec[i]
stat[i,] <- cbind(alpha,alphas0(alpha,b,p,n))}


if( stat[nrow(stat),2] < b[1] ) alphau <- 1
if( stat[1,2] >= b[1]) alphau <- -1
if( stat[1,2] < b[1] & stat[nrow(stat),2] >= b[1])
{
    for (i in 2:length(alphavec))
    {low <- stat[i-1,2]; high <- stat[i,2]
    if(low < b[1] & high > b[1]) stat1 <- stat[(i-1):i,]
    }
    base1 <- stat1[1,2]; base2 <- stat1[1,1]
    slope <- (stat1[2,1]-stat1[1,1])/(stat1[2,2]-stat1[1,2])
    alphau <- (b[1] - base1)*slope+base2
}   

return(alphau)
}
