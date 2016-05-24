condnorm <-
function(xx, mu1=0, sd=1, c1) {
   # this is the pdf of z1, with mean mu1, conditional on abs(z1) > c1
   denom <- 1 - (pnorm(c1, mean=mu1, sd=sd) - pnorm(-c1, mean=mu1, sd=sd) )
   posvals.ind <- {{xx < -c1} | {xx > c1}}
   result <- xx*0
   result[posvals.ind] <- dnorm(xx[posvals.ind], mean=mu1)/denom
   result
   }
