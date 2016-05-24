# hilbe.NBR2.F8.2.r
# Negative binomial regression distributions with 
#    user specified series of mean values for a specified alpha value
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# Figure 8.2-6   default: means at .5, 1, 2, 5, 10 and alpha=.33
#
obs <- 11
mu <- c(0.5,1,2,5,10)
y <- 0:10
alpha <- .33
amu <- mu*alpha
layout(1) 
for (i in 1:length(mu)) { 
ynb2 = exp(
        y*log(amu[i]/(1+amu[i])) 
     - (1/alpha)*log(1+amu[i]) 
     + log( gamma(y +1/alpha) )
     - log( gamma(y+1) ) 
     - log( gamma(1/alpha) ) 
)
if (i==1) { 
plot(y, ynb2, col=i, type='l', lty=i) 
  } else { 
  lines(y, ynb2, col=i, lty=i) 
  } 
}





