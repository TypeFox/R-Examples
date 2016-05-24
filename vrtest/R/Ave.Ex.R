Ave.Ex <-
function(y)
{
T <- length(y)
B <- rep(0,T-1)
SIGMAT <- sum(y^2)/(T-1)

CPI <- seq(-0.8,0.8,0.01)
EXLMINF <- 0
EXLRINF <- 0

for (j in 1:161)
{
    B[1] <- y[1]
    for (k in 2:(T-1))
    {
    B[k] <- y[k] + CPI[j]*B[k-1]
    }

SIGMAH <- SIGMAT - (( sum(y[2:T]*B) )^2)/T/(sum(B^2))

LM <- (sum(y[2:T]*B))^2/T*(1 - CPI[j]^2)/SIGMAT/SIGMAT
LR <- T*log(SIGMAT/SIGMAH)
EXLMINF <- EXLMINF + exp(LM/2)/161;
EXLRINF <- EXLRINF + exp(LR/2)/161;
}

EXLMINF <- log(EXLMINF)
EXLRINF <- log(EXLRINF)

return(list(Ex.LM=EXLMINF,Ex.LR=EXLRINF))
}
