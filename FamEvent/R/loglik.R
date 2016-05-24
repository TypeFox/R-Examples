loglik<- function(theta, data, design, base.dist, agemin)
{
#data <- data[data$currentage>agemin,]
lambda   <- exp(theta[1])
rho  <- exp(theta[2])
beta.sex <- theta[3]
beta.gen <- theta[4]


time0 <- data$time-agemin
status<- data$status
wt <- data$weight
xbeta <- beta.sex*data$gender+beta.gen*data$mgene

bhaz <- hazards(base.dist, time0, c(lambda,rho))
bcumhaz <- cumhaz(base.dist, time0, c(lambda,rho))

H <- bcumhaz*exp(xbeta)
logh <- log(bhaz) + xbeta

loglik <-  wt * (- H + status*logh )
loglik[data$time<=agemin] <- 0

# Ascertainment correction by design

if(design=="cli" | design=="cli+"){
  i.m <- data$generation==1 & data$gender==0
  i.f <- data$generation==1 & data$gender==1
  ip <-  data$generation==2  & data$status==1
}
else  ip <- data$proband==1


cagep <- data$currentage[ip]-agemin
xbeta.p <- beta.sex*data$gender[ip]+beta.gen*data$mgene[ip]
bcumhaz.p <- cumhaz(base.dist, cagep,c(lambda,rho))
wt.p <- data$weight[ip]

slogasc <- sum(wt.p*log(1-exp(-bcumhaz.p*exp(xbeta.p)) ), na.rm=T) 

if(design=="cli" | design=="cli+"){
cage.m <- data$currentage[i.m]-agemin
xbeta.m <- beta.sex*data$gender[i.m]+beta.gen*data$mgene[i.m]
bcumhaz.m <- cumhaz(base.dist, cage.m, c(lambda,rho))

cage.f <- data$currentage[i.f]-agemin
xbeta.f <- beta.sex*data$gender[i.f]+beta.gen*data$mgene[i.f]
bcumhaz.f <- cumhaz(base.dist, cage.f,c(lambda,rho))

wt.p <- data$weight[data$proband==1]
slogasc <- slogasc + sum(wt.p*log(1-exp(-bcumhaz.m*exp(xbeta.m) -bcumhaz.f*exp(xbeta.f) ) ),na.rm=T)
}

likelihood  <- try(sum(loglik, na.rm=T) - slogasc)
return(-likelihood)
}