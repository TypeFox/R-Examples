disc2 <-
function (x, k,out=c("symb","num")) 
{
    n = length(x)
    ciclo = ceiling(n/k)
    y = x
    if(out=="num")
     {for (i in 1:(k - 1)) {
        y[order(x)[((i - 1) * ciclo + 1):(i * ciclo)]] = i
    }
    y[order(x)[((k - 1) * ciclo + 1):n]] = k
    return(y)
           }
     else{
cutpoints=(1:k-1)*ciclo
cutpoints=cutpoints[-1]
cutpoints=c(-Inf,.5*(y[order(x)[cutpoints]]+y[order(x)[cutpoints+1]]),Inf)
y=cut(y,cutpoints)
return(y)
}
}
