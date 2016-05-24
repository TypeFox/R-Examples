# seed the algorithm  
m <- mean(faithful$eruptions)
s <- sd(faithful$eruptions)

mle <- nlmax(loglik,p=c(0.5,m,m,s,s),x=faithful$eruptions)$estimate
mle
