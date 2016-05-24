
varscore.gb2 <- function(x, shape1, scale, shape2, shape3, w=rep(1, length(x)), hs=rep(1,length(x))){
m <- length(x)
Vsc <- matrix(rep(0,16), ncol=4)
for (i in 1:m){
wsci <- w[i]*dlogf.gb2(x[i], shape1, scale, shape2, shape3)
Wsci <- wsci%*%t(wsci)
Vsc <- Vsc + (hs[i]^2)*Wsci
}
return(Vsc)
}

vepar.gb2 <- function(x, Vsc, shape1, scale, shape2, shape3, w=rep(1, length(x)), hs=rep(1,length(x))){
#estimated variance-covariance matrix of af, bf, pf and qf (EVCM)
m <- length(x)
#the left and right side of the sandwich estimator 
# = - (sum_{i=1}^m w_i * n_i * h_i )
#where h_i is the matrix of second derivatives of the log density
#in this case obtained by the function d2logf.gb2

WDM <- matrix(rep(0,16), ncol=4)
for (i in 1:m){
WDM <- WDM + w[i]*hs[i]*d2logf.gb2(x[i], shape1, scale, shape2, shape3)
}
WDM <- -WDM

#sandwich estimator
EVCM <- solve(WDM)%*%Vsc%*%solve(WDM)
return(list(EVCM,WDM))
}

#numerical derivatives of the indicators

derivind.gb2 <- function(shape1, scale, shape2, shape3){

par <- c(shape1, scale, shape2, shape3)

med <- function(par){
main.gb2(0.6,par[1],par[2],par[3],par[4])[1]
}
mean <- function(par){
main.gb2(0.6,par[1],par[2],par[3],par[4])[2]
}
arpr <- function(par){
main.gb2(0.6,par[1],par[2],par[3],par[4])[3]
}
rmpg <- function(par){
main.gb2(0.6,par[1],par[2],par[3],par[4])[4]
}
qsr <- function(par){
main.gb2(0.6,par[1],par[2],par[3],par[4])[5]
}
gini <- function(par){
main.gb2(0.6,par[1],par[2],par[3],par[4])[6]
}

dmed  <- grad(med, par, method= "Richardson", method.args = list())
dmean <- grad(mean, par, method= "Richardson", method.args = list())
darpr <- grad(arpr, par, method= "Richardson", method.args = list())
drmpg <- grad(rmpg, par, method= "Richardson", method.args = list())
dqsr  <- grad(qsr, par, method= "Richardson", method.args = list())
dgini <- grad(gini, par, method= "Richardson", method.args = list())

MFDI <- rbind(dmed, dmean, darpr, drmpg, dqsr, dgini, deparse.level=0)
return(MFDI) 
}

veind.gb2 <- function(Vpar, shape1, scale, shape2, shape3){

# matrix of first derivatives of the indicators
MFDI <- derivind.gb2(shape1, scale, shape2, shape3)
#delta method for variance estimation of the indicators
# variance-covariance matrix of the indicators
IVCM <- matrix(rep(0, 36), ncol=6)
for (i in 1:6){
for (j in 1:6){
IVCM[i,j] <- t(MFDI[i,])%*%Vpar%*%MFDI[j,]
}
}
return(IVCM)
}


