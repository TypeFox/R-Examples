###############################################
# Section 8.3 A One-Sided Test of a Normal Mean
###############################################

library(LearnBayes)

 pmean=170; pvar=25
 probH=pnorm(175,pmean,sqrt(pvar))
 probA=1-probH
 prior.odds=probH/probA
 prior.odds

 weights=c(182, 172, 173, 176, 176, 180, 173, 174, 179, 175)
 xbar=mean(weights)
 sigma2=3^2/length(weights)

 post.precision=1/sigma2+1/pvar
 post.var=1/post.precision

 post.mean=(xbar/sigma2+pmean/pvar)/post.precision
 c(post.mean,sqrt(post.var))

 post.odds=pnorm(175,post.mean,sqrt(post.var))/
  (1-pnorm(175,post.mean,sqrt(post.var)))
 post.odds

 BF = post.odds/prior.odds
 BF

 postH=probH*BF/(probH*BF+probA)
 postH

 z=sqrt(length(weights))*(mean(weights)-175)/3
 1-pnorm(z)

 weights=c(182, 172, 173, 176, 176, 180, 173, 174, 179, 175)
 data=c(mean(weights),length(weights),3)
 prior.par=c(170,1000)
 mnormt.onesided(175,prior.par,data)
