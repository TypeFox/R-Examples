VAR.adjust <-
function(b,bias,p,type)
{
k <- nrow(b)
bs1 <- b - bias
delta <- VAR.modul(bs1,p)[1]

if(delta < 1) bs2 <- bs1[,1:(k*p),drop=FALSE] #return(bs1)
if(delta >= 1)
{
delta1 <- 1
while(delta >= 1)
    {
        delta1 <- delta1-0.01
        bias <- delta1*bias
        bs2 <- b[,1:(k*p),drop=FALSE] - bias[,1:(k*p),drop=FALSE]
        if(is.nan(sum(bs2)) | is.infinite(sum(bs2)) )
        {bs2 <- b[,1:(p*k),drop=FALSE]; break}
        delta <- VAR.modul(bs2,p)[1]
    }
}

if(type=="const" | type=="const+trend") bs2 <- cbind(bs2,bs1[,(p*k+1):ncol(b),drop=FALSE])

   #intercept adjustment
   # index <- 1:k + 1
   # sum1 <- 0
   # for(i in 1:p){
   # sum1 <- sum1 + bs2[,index]
   # index <- index+ k}
    #bs2[,1] <- (diag(k)-sum1) %*% matrix(colMeans(x))
return(bs2)
}
