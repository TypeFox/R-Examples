DLtest <-
function(y,p)
{
ym <- as.matrix(y-mean(y))
n <- nrow(ym); s2 <- sum(ym^2)/(n-p)

sum3 <- numeric(n-p)
sum2 <- 0
for(j in (p+1):n) {
    sum1 <- 0
    for(i in (p+1):n){
    indicate <- 0
    zi <- ym[(i-1):(i-p),1]
    zj <- ym[(j-1):(j-p),1]
    tem1 <- as.numeric(zi <= zj)
    if( prod(tem1) == 1) indicate <- 1
    sum1 <- sum1 + ym[i,1]*indicate
    }
    sum2 <- sum2 + sum1^2
    sum3[j-p] <- abs(sum1/sqrt(n-p))
}

Cp <- sum2/(s2*(n-p)^2)
Kp <- max(sum3)/sqrt(s2)

return(list(Cpstat=Cp,Kpstat=Kp))
}
