# hilbe.NBR2.F6.3.r
# Table 6.16
# Graphic of Poisson distributions with user specified mean values
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# Table 6.16; Figure 6.3  
# 
m<- c(0.5,1,3,5,7,9) #Poisson means 
y<- 0:19             #Observed counts 
layout(1) 
for (i in 1:length(m)) { 
  p<- dpois(y, m[i]) #poisson pmf 
  if (i==1) { 
  plot(y, p, col=i, type='l', lty=i) 
  } else { 
  lines(y, p, col=i, lty=i) 
  } 
}



