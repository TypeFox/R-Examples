loglikem<- function(theta, theta0, data, design, base.dist, agemin, vec=TRUE)
{

theta[1:2] <- exp(theta[1:2])
theta0[1:2] <- exp(theta0[1:2])
beta.sex <- theta[3]
beta.gen <- theta[4]


time0 <- data$time-agemin
status<- data$status
wt <- data$weight

xbeta1 <- beta.sex*data$gender+beta.gen*1
xbeta0 <- beta.sex*data$gender+beta.gen*0

bhaz <- hazards(base.dist, time0, theta[1:2])
bcumhaz <- cumhaz(base.dist, time0, theta[1:2])

H1 <- bcumhaz*exp(xbeta1)
H0 <- bcumhaz*exp(xbeta0)
logh1 <- log(bhaz) + xbeta1
logh0 <- log(bhaz) + xbeta0


  p1 <- cprob(theta0, data=data, mut=1, base.dist=base.dist, agemin=agemin)
  p0 <- cprob(theta0, data=data, mut=0, base.dist=base.dist, agemin=agemin)
  ex1 <- p1/(p1+p0) #P(x=1|Xp, y)=P(y|x=1)*P(x=1|Xp)/(p1+p0) for EM
  #ex1[!is.na(data$mgene)] <- data$mgene[!is.na(data$mgene)]
  

  loglik <- wt * (- H1 + status*logh1 ) *ex1 + wt * (- H0 + status*logh0 ) * (1-ex1)
  loglik[data$time<=agemin] <- 0
  
# Ascertainment correction by design

ip <- data$proband==1
cagep <- data$currentage[ip]-agemin
xbeta.p <- beta.sex*data$gender[ip]+beta.gen*data$mgene[ip]
bcumhaz.p <- cumhaz(base.dist, cagep, theta[1:2])
wt.p <- data$weight[ip]

slogasc.p <- wt.p*log(1-exp(-bcumhaz.p*exp(xbeta.p))) 

if(design=="cli" | design=="cli+"){
  
  i.m <- data$generation==1 & data$gender==0
  i.f <- data$generation==1 & data$gender==1
  i.s <- data$generation==2 & data$proband==0 & data$status==1
  
cage.m <- data$currentage[i.m]-agemin
xbeta.m0 <- beta.sex*data$gender[i.m]+beta.gen*0
xbeta.m1 <- beta.sex*data$gender[i.m]+beta.gen*1
bcumhaz.m <- cumhaz(base.dist, cage.m, theta[1:2])

cage.f <- data$currentage[i.f]-agemin
xbeta.f0 <- beta.sex*data$gender[i.f]+beta.gen*0
xbeta.f1 <- beta.sex*data$gender[i.f]+beta.gen*1
bcumhaz.f <- cumhaz(base.dist, cage.f, theta[1:2])
                    
cage.s <- data$currentage[i.s]-agemin
xbeta.s0 <- beta.sex*data$gender[i.s]+beta.gen*0
xbeta.s1 <- beta.sex*data$gender[i.s]+beta.gen*1
bcumhaz.s <- cumhaz(base.dist, cage.s, theta[1:2])
                                        
wt.m <- data$weight[i.m]
wt.f <- data$weight[i.f]
wt.s <- data$weight[i.s]

loglik <- wt * (- H1 + status*logh1 ) *ex1 + wt * (- H0 + status*logh0 ) * (1-ex1)

slogasc.m <-  wt.m*log(1-exp(-bcumhaz.m*exp(xbeta.m1)))*ex1[i.m] + wt.p*log(1-exp(-bcumhaz.m*exp(xbeta.m0)))*(1-ex1[i.m])
slogasc.f <-  wt.f*log(1-exp(-bcumhaz.f*exp(xbeta.f1)))*ex1[i.f] + wt.p*log(1-exp(-bcumhaz.f*exp(xbeta.f0)))*(1-ex1[i.f]) 
slogasc.s <-  wt.s*log(1-exp(-bcumhaz.s*exp(xbeta.s1)))*ex1[i.s] + wt.p*log(1-exp(-bcumhaz.s*exp(xbeta.s0)))*(1-ex1[i.s])

loglik[i.m] <- loglik[i.m] - slogasc.m
loglik[i.f] <- loglik[i.f] - slogasc.f
loglik[i.s] <- loglik[i.s] - slogasc.s

}

loglik[ip] <- loglik[ip] - slogasc.p

if(vec) return(-loglik)
else return(-sum(loglik, na.rm=T) )
            
}