"discrim" <-
function(x)
{
x <- na.exclude(as.matrix(x))
k <- ncol(x)
N <- nrow(x)
ni <- as.integer(N/3)
TOT <- apply(x, 1, mean)
tmpx <- cbind(x,TOT)[order(TOT),]
tmpxU <- tmpx[(N+1-ni):N,]
tmpxL <- tmpx[1:ni,]
Ui <- apply(tmpxU,2,sum)
Li <- apply(tmpxL,2,sum)
discrim <- (Ui - Li)/ni
return (discrim[1:k])
}

