# hilbe.NBR2.F8.1.r
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# Negative binomial regression distributions with 
#    user specified series of alpha values for a specified mean value
# Figure 8.2-6   default: means at .5, 1, 2, 5, 10 and alpha=0
#
m<- c(0.5,1,2,5,10)  #mean values 
y<- 0:10             #Observed counts 
layout(1) 
for (i in 1:length(m)) { 
  p<- dpois(y, m[i]) #poisson pmf 
  if (i==1) { 
  plot(y, p, col=i, type='l', lty=i) 
  } else { 
  lines(y, p, col=i, lty=i) 
  } 
}




