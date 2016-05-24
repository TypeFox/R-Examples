Haz4 <- function(d, S, gam){
#### This function does one iteration of the baseline hazard Hw calculation.
#### All input should be in the order of y and -d. 
#### (-d as a tie breaker)
#### S is baseline survival prob, vector of length n; 
#### d is censoring indecator.
#### gam = (gamma1i, gamma2i), determined by beta Z. A matrix 
#### of  nx2.   
#### This function is the same as Haz2() ? ! yes.

gam12 <- as.vector(gam[,1] - gam[,2])
   temp <- gam12*S
   gweight <- (1- d*temp)/(gam[,2] + temp)
   ## n <- length(d)                3/2015  MZ
   ## Hw <- rep(0, n)
   ttemp <- cumsumsurv(gweight)    ###  rev( cumsum( rev(gweight) )) 3/2015  MZ
   Hw <- d/ttemp
###################### change Nov. 2013 
#for(k in 1:n) {
#     if(d[k] ==1) Hw[k] <- 1/( sum(gweight[k:n]) )
#              }
###################### get rid of for loop.
list( Hazw=Hw, Su=exp(-cumsum(Hw)) )
}