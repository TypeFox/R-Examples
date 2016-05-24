YP4 <- function(y, d, Z, b1, b2, k){
#### y, d, Z  are the data inputs.(times, status and covariates)
#### b1, b2 needs to be given.
#### b1 = (beta1) and b2=(beta2).
#### Z is matrix of covariates, size nx(p+1); Z=(X_i, Z_i)
#### the dim(X_i) = p and also the number of length(alpha)=p
#### k(+1) is the number of iterations we do, in maximizing
#### the baseline.
####  May be add an optional input of alpha? could be from a Cox model initial analysis??? May be a better S0???
####  S0 = WKM()$Surv


yorder <- order(y, -d)
ysort <- y[yorder]
dsort <- d[yorder]
Zsort <- as.matrix(Z[yorder,])  ## could also use Wdataclean5( )

dimZ <- dim(Zsort)
n <- dimZ[1]        #### length of y, d should be also n
S0 <- (n:1)/(n+1)   #### is there a better initial value?
alpha <- rep(0, dimZ[2] -1)   #### if there is no better initial value.

ab1 <- c(alpha, b1)
ab2 <- c(alpha, b2)
gam <- exp( - Zsort %*% cbind(ab1,ab2) )
tempT <- Haz4(d=dsort, S=S0, gam=gam)      #### change made on Nov. 2013. Do Haz4() first, then survreg().

eZb2 <- exp( Zsort[, dimZ[2]] * b2 )
eZb1Mb2 <- exp( Zsort[, dimZ[2]] * (b1-b2) )
Ais <-  eZb2 * log(1+ eZb1Mb2*(1-tempT$Su)/tempT$Su )            #### could do Haz4(  ) first, with alpha =0 etc.
tempTT <- survreg(Surv(Ais, dsort) ~ Zsort[,-dimZ[2]] -1, dist= "expo")
alpha <- tempTT$coef
names(alpha) <- NULL
alpha <- (-alpha)



#tempT <- Haz4(d=dsort, S=temp3$Su, gam=gam)

#####par(new=TRUE)
#####plot(tempT$Su,col="blue", ylim=c(0,1))

i <- 1
err <- 8
while( i <= k & err > 0.0000001 ) {
alphaOld <- alpha
ab1 <- c(alpha, b1)
ab2 <- c(alpha, b2)
gam <- exp( - Zsort %*% cbind(ab1, ab2) )
tempT <- Haz4(d=dsort, S=tempT$Su, gam=gam)
Ais <- eZb2 * log( 1+ eZb1Mb2*(1-tempT$Su)/tempT$Su)
tempTT <- survreg(Surv(Ais, dsort) ~ Zsort[,-dimZ[2]] -1, dist= "expo")
alpha <- tempTT$coef
names(alpha) <- NULL
alpha <- -alpha
i <- i+1
err <- abs(alpha - alphaOld)
}

if( i == (k+1) ) warning("check convergence of alpha")
ab1 <- c(alpha, b1)
ab2 <- c(alpha, b2)
gam <- exp( - Zsort %*% cbind(ab1, ab2) )
tempT <- Haz4(d=dsort, S=tempT$Su, gam=gam)

list(d=dsort, Hazw=tempT$Hazw, Survival=tempT$Su, gam=gam, alpha = alpha)
}
