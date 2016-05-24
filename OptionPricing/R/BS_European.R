
BS_EC <- function(T=0.25,K=100,r=0.05,sigma=0.2,S0=100){
# Black-Scholes Formula for European Call option
# K ... strike Price
# r ... riskfree rate
# sigma ... (yearly volatility)
# T ... time horizon
# S0 ... stockprice at start
d1 <- (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T));
d2 <- d1-sigma*sqrt(T);
phid1 <- pnorm(d1)
price <- S0*phid1-K*exp(-r*T)*pnorm(d2);
delta <- phid1
gamma <- delta/(S0*sigma*sqrt(T))
res<- c(price,delta,gamma)
names(res) <- c("price","delta","gamma")
return(res)
}
 
BS_EP <- function(T=0.25,K=100,r=0.05,sigma=0.2,S0=100){
# Black-Scholes Formula for European Put option
# K ... strike Price
# r ... riskfree rate
# sigma ... (yearly volatility)
# T ... time horizon
# S0 ... stockprice at start
d1<- (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T));
d2 <- d1-sigma*sqrt(T);
phimd1 <- pnorm(-d1)
price <- -S0*phimd1 + K*exp(-r*T)*pnorm(-d2);
delta <- -phimd1
phid1 <- 1-phimd1
gamma <- phid1/(S0*sigma*sqrt(T))
res<- c(price,delta,gamma)
names(res) <- c("price","delta","gamma")
return(res)
}




