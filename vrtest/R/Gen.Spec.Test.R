Gen.Spec.Test <-
function(y,B=300)
{
set.seed(12345)
n<- length(y)
e <- y - mean(y)
v <- var(e)
y1 <- y[1:(n-1)]
weiexp <- compweexp(y1)
CvMexp <- 0

for(j in 1:(n-1)) {
   aux2 <- 1/((j*pi)^2)
   aux2 <- aux2/(n-j+1)
   CvMexp <- CvMexp+ aux2* t(e[(1+j):n]) %*% weiexp[1:(n-j),1:(n-j)] %*% e[(1+j):n]
}
CvMexp <- CvMexp/v

CvMexpb <- matrix(0,nrow=B,ncol=2)
for(k in 1:B)
{
    eb <- e * Mammen(n)
    eb <- eb - mean(eb)
    tem <- 0
    for( j in 1:(n-1) ){
    aux2 <- 1/((j*pi)^2)
    aux2 <- aux2/(n-j+1)
    tem <- tem+aux2* t(eb[(1+j):n]) %*% weiexp[1:(n-j),1:(n-j)] %*% eb[(1+j):n]
    }
    CvMexpb[k,] <- cbind(tem/v > CvMexp,tem/v)
    
}    

pboot <- mean(CvMexpb[,1])
Critboot <- quantile(CvMexpb[,2],c(0.9,0.95,0.99))
return(pboot)
}
