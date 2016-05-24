RappN <-
function(x,res,delta){ 
# computes semi-empirical cdf of 2*rho_0(u)=u^2 at x^2
n  <- length(delta); FU <- FC <- 0; nu  <- sum(delta); nc <- n-nu
rU <- res[delta==1]
rC <- res[delta==0]
if (nu > 0)         FU <- sum( abs(rU) <= x )
if (nc > 0 & x > 8) FC <- nc
if (nc > 0 & x < 8) { 
 t     <- rep(0,nc); i8 <- rC < 8
 num   <- (pnorm(x) - pnorm(pmax(rC,-x)))*(rC <= x)
 den   <- 1-pnorm(rC)
 t[i8] <- num[i8]/den[i8] 
 FC    <- sum(t) }
Fnz <- (FC+FU)/n
Rap <- Fnz/pchisq(x^2,df = 1)
c(Fnz,Rap)}

