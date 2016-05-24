RappW <-
function(x,res,delta){ 
# computes semi-empirical cdf of rho_0(u)=exp(u)-u at exp(x)-x
n  <- length(delta); FU <- FC <- 0; nu  <- sum(delta); nc <- n-nu
rU <- res[delta==1]
rC <- res[delta==0]
if (nu > 0)         FU <- sum( rU <= x & rU > Izero(x) )
if (nc > 0 & x > 3) FC <- nc
if (nc > 0 & x < 3) { 
 tmp     <- rep(0,nc); i3 <- rC < 3
 num     <- (plweibul(x) - plweibul(pmax(rC,Izero(x))))*(rC <= x)
 den     <- 1-plweibul(rC)
 tmp[i3] <- num[i3]/den[i3]
 tmp[is.na(tmp)] <- 1
 FC    <- sum(tmp) }
Fnz <- (FC+FU)/n
Rap <- Fnz/(plweibul(x)-plweibul(Izero(x)))
c(Fnz,Rap)}

