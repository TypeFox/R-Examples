

profileDG <- function(fit, steps=100, wh=1:p){
  Pnames <- names(B0 <- coefficients(fit))
  pv0 <- t(as.matrix(B0))
  p <- length(Pnames)
  summ <- summary(fit)
  std.err <- summ$coefficients[, "Std. Error"]
  mf <- model.frame(fit)
  n <- length(Y <- model.response(mf))
  Off <- model.offset(mf)
  if (!length(Off)) Off <- rep(0, n)
  W <- model.weights(mf)
  if (length(W) == 0) W <- rep(1, n)
  OriginalDeviance <- deviance(fit)
  DispersionParameter <- summ$dispersion
  X <- model.matrix(fit)
  fam <- family(fit)
  switch(fam$family, binomial = {
    if (!is.null(dim(Y))) {
      n <- n/2
      Off <- Off[1:n]
      Y <- Y[, 1]/(W <- drop(Y %*% c(1, 1)))
    }
    #zmax <- sqrt(qchisq(1 - alpha/2, p))
    profName <- "z"
    }, poisson = , "Negative Binomial" = {
      #zmax <- sqrt(qchisq(1 - alpha/2, p))
      profName <- "z"
    }, gaussian = , quasi = , inverse.gaussian = , quasibinomial = ,
      quasipoisson = , {
      #zmax <- sqrt(p * qf(1 - alpha/2, p, n - p))
      profName <- "tau"
  })
  prof <- vector("list", length = length(wh))
  names(prof) <- Pnames[wh]
  for (i in wh) {
    zi <- 0
    pvi <- pv0
    Xi <- X[, -i, drop = FALSE]
    pi <- Pnames[i]
    for (sgn in c(-1, 1)) {
      z <- 0
      LP <- X %*% fit$coefficients + Off
      bi <- B0[i] + sgn * 3 * std.err[i]
      if (bi < -25) bi <- -25
      if (bi > 25) bi <- 25
      brange <- seq(B0[i], bi, length.out=steps)[-1]
      for (k in 1:length(brange)){
        off <- Off + X[, i] * brange[k]
        fm <- glm.fit(x = Xi, y = Y, weights = W, etastart = LP,
                  offset = off, family = fam, control = fit$control)
        LP <- Xi %*% fm$coefficients + off
        ri <- pv0
        ri[, names(coef(fm))] <- coef(fm)
        ri[, pi] <- brange[k]
        pvi <- rbind(pvi, ri)
        zz <- (fm$deviance - OriginalDeviance)/DispersionParameter
        if (zz > -0.001) zz <- max(zz, 0) else stop("profiling has found a better solution, so original fit had not converged")
        z <- sgn * sqrt(zz)
        zi <- c(zi, z)
      }
    }
    si <- order(zi)
    prof[[pi]] <- structure(data.frame(zi[si]), names = profName)
    prof[[pi]]$par.vals <- pvi[si, ]
  }
  val <- structure(prof, original.fit = fit, summary = summ)
  class(val) <- c("profile.glm", "profile")
  val
}



##########################################


Poisson.ratio<-function(x, y, conf.level=0.95, alternative="two.sided")
{


alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

Y<-c(x,y)

if( any(Y-as.integer(Y) != 0) )
 {warning("There are non-integer values in the data!")}

nx<-length(x)
ny<-length(y)
X<-as.factor(rep( c(1,0), c(nx,ny) ))
dat<-data.frame(Y=Y,X=X)

fit <- glm(Y~X, data=dat, family=poisson(link="log"))

profit <- profileDG(fit, steps=100)


switch(alternative,
two.sided={

lCItemp <- confint(profit, level= conf.level)
lest <- coefficients(fit)[2]

if(lest<=-25){lest<-(-Inf)}
if(lest>=25){lest<-(Inf)}
estimate <- exp(lest)

if(is.na(lCItemp[2,1])){lower<-(-Inf)}
else{
if(lCItemp[2,1]<=-25){lower<-(-Inf)}
else{
if(lCItemp[2,1]>=25){lower<-Inf}
else{lower<-lCItemp[2,1]}
}}

if(is.na(lCItemp[2,2])){upper<-Inf}
else{
if(lCItemp[2,2]<=-25){upper<-(-Inf)}
else{
if(lCItemp[2,2]>=25){upper<-Inf}
else{upper<-lCItemp[2,2]}
}}

conf.int <- exp(c(lower,upper))

},


less={
lCItemp <- confint(profit, level= 1-(1-conf.level)*2)
lest <- coefficients(fit)[2]
if(lest<=-25){lest<-(-Inf)}
if(lest>=25){lest<-(Inf)}
estimate <- exp(lest)

if(is.na(lCItemp[2,2])){upper<-(Inf)}
else{
if(lCItemp[2,2]<=-25){upper<-(-Inf)}
else{
if(lCItemp[2,2]>=25){upper<-Inf}
else{upper<-lCItemp[2,2]}
}}

conf.int <- exp(c(-Inf,upper))

},

greater={

lCItemp <- confint(profit, level= 1-(1-conf.level)*2)
lest <- coefficients(fit)[2]
if(lest<=-25){lest<-(-Inf)}
if(lest>=25){lest<-(Inf)}
estimate <- exp(lest)

if(is.na(lCItemp[2,1])){lower<-Inf}
else{
if(lCItemp[2,1]<=-25){lower<-(-Inf)}
else{
if(lCItemp[2,1]>=25){lower<-Inf}
else{lower<-lCItemp[2,1]}
}}

conf.int <- exp(c(lower, Inf))

})


METHOD <- "Ratio of means assuming Poisson distribution based on a likelihood profile of a GLM fit with log-link."

attr(conf.int, which="methodname")<-METHOD

return(
list(conf.int=conf.int,
estimate=estimate)  
) 

}


#Poisson.ratio(x=c(0,0,0,1,0), y=c(2,5,3,1,1))
#Poisson.ratio(x=c(0,0,0,0,0), y=c(2,5,3,1,1))
#Poisson.ratio(x=c(1), y=c(2))


############################################


Quasipoisson.ratio<-function(x, y, conf.level=0.95, alternative="two.sided")
{

alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

Y<-c(x,y)

if( any(Y-as.integer(Y) != 0) )
 {warning("There are non-integer values in the data!")}

nx<-length(x)
ny<-length(y)

if(any(c(nx,ny)<2))
 {stop("There are less than 2 observations in at least one group: variance can not be estimated")}

X<-as.factor(rep( c(1,0), c(nx,ny) ))
dat<-data.frame(Y=Y,X=X)

fit <- glm(Y~X, data=dat, family=quasipoisson(link="log"))

profit <- profileDG(fit, steps=100)

switch(alternative,
two.sided={

lCItemp <- confint(profit, level= conf.level)
lest <- coefficients(fit)[2]

if(lest<=-25){lest<-(-Inf)}
if(lest>=25){lest<-(Inf)}
estimate <- exp(lest)

if(is.na(lCItemp[2,1])){lower<-(-Inf)}
else{
if(lCItemp[2,1]<=-25){lower<-(-Inf)}
else{
if(lCItemp[2,1]>=25){lower<-Inf}
else{lower<-lCItemp[2,1]}
}}

if(is.na(lCItemp[2,2])){upper<-Inf}
else{
if(lCItemp[2,2]<=-25){upper<-(-Inf)}
else{
if(lCItemp[2,2]>=25){upper<-Inf}
else{upper<-lCItemp[2,2]}
}}

conf.int <- exp(c(lower,upper))

},


less={
lCItemp <- confint(profit, level= 1-(1-conf.level)*2)
lest <- coefficients(fit)[2]
if(lest<=-25){lest<-(-Inf)}
if(lest>=25){lest<-(Inf)}
estimate <- exp(lest)

if(is.na(lCItemp[2,2])){upper<-(Inf)}
else{
if(lCItemp[2,2]<=-25){upper<-(-Inf)}
else{
if(lCItemp[2,2]>=25){upper<-Inf}
else{upper<-lCItemp[2,2]}
}}

conf.int <- exp(c(-Inf,upper))

},

greater={

lCItemp <- confint(profit, level= 1-(1-conf.level)*2)
lest <- coefficients(fit)[2]
if(lest<=-25){lest<-(-Inf)}
if(lest>=25){lest<-(Inf)}
estimate <- exp(lest)

if(is.na(lCItemp[2,1])){lower<-Inf}
else{
if(lCItemp[2,1]<=-25){lower<-(-Inf)}
else{
if(lCItemp[2,1]>=25){lower<-Inf}
else{lower<-lCItemp[2,1]}
}}

conf.int <- exp(c(lower, Inf))

})

METHOD <- "Ratio of means based on a quasi-likelihood profile of a quasipoisson GLM fit with log-link."

attr(conf.int, which="methodname")<-METHOD

return(
list(conf.int=conf.int,
estimate=estimate)  
) 

}


#Quasipoisson.ratio(x=c(0,0,0,1,0), y=c(2,5,3,1,1))
#Quasipoisson.ratio(x=c(0,0,0,0,0), y=c(2,5,3,1,1))
#Quasipoisson.ratio(x=c(1), y=c(2))


#Quasipoisson.ratio(x=c(0,0,0,1,0), y=c(2,5,3,1,1), alternative="less")
#Quasipoisson.ratio(x=c(0,0,0,0,0), y=c(2,5,3,1,1), alternative="less")
#Quasipoisson.ratio(x=c(1), y=c(2), alternative="less")


#Quasipoisson.ratio(x=c(0,0,0,1,0), y=c(2,5,3,1,1), alternative="greater")
#Quasipoisson.ratio(x=c(0,0,0,0,0), y=c(2,5,3,1,1), alternative="greater")
#Quasipoisson.ratio(x=c(1), y=c(2), alternative="greater")


############################


Negbin.ratio<-function(x, y, conf.level=0.95, alternative="two.sided")
{


alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))

Y<-c(x,y)

if( any(Y-as.integer(Y) != 0) )
 {warning("There are non-integer values in the data!")}

nx<-length(x)
ny<-length(y)


if(any(c(nx,ny)<2))
 {stop("There are less than 2 observations in at least one group: variance can not be estimated")}

X<-as.factor(rep( c(1,0), c(nx,ny) ))
dat<-data.frame(Y=Y,X=X)

fit <- glm.nb(Y~X, data=dat)

profit <- profileDG(fit, steps=100)


switch(alternative,
two.sided={

lCItemp <- confint(profit, level= conf.level)
lest <- coefficients(fit)[2]

if(lest<=-25){lest<-(-Inf)}
if(lest>=25){lest<-(Inf)}
estimate <- exp(lest)

if(is.na(lCItemp[2,1])){lower<-(-Inf)}
else{
if(lCItemp[2,1]<=-25){lower<-(-Inf)}
else{
if(lCItemp[2,1]>=25){lower<-Inf}
else{lower<-lCItemp[2,1]}
}}

if(is.na(lCItemp[2,2])){upper<-Inf}
else{
if(lCItemp[2,2]<=-25){upper<-(-Inf)}
else{
if(lCItemp[2,2]>=25){upper<-Inf}
else{upper<-lCItemp[2,2]}
}}

conf.int <- exp(c(lower,upper))

},


less={
lCItemp <- confint(profit, level= 1-(1-conf.level)*2)
lest <- coefficients(fit)[2]
if(lest<=-25){lest<-(-Inf)}
if(lest>=25){lest<-(Inf)}
estimate <- exp(lest)

if(is.na(lCItemp[2,2])){upper<-(Inf)}
else{
if(lCItemp[2,2]<=-25){upper<-(-Inf)}
else{
if(lCItemp[2,2]>=25){upper<-Inf}
else{upper<-lCItemp[2,2]}
}}

conf.int <- exp(c(-Inf,upper))

},

greater={

lCItemp <- confint(profit, level= 1-(1-conf.level)*2)
lest <- coefficients(fit)[2]
if(lest<=-25){lest<-(-Inf)}
if(lest>=25){lest<-(Inf)}
estimate <- exp(lest)

if(is.na(lCItemp[2,1])){lower<-Inf}
else{
if(lCItemp[2,1]<=-25){lower<-(-Inf)}
else{
if(lCItemp[2,1]>=25){lower<-Inf}
else{lower<-lCItemp[2,1]}
}}

conf.int <- exp(c(lower, Inf))

})


METHOD <- "Ratio of means assuming negative binomial distribution based on a likelihood profile of a GLM fit with log-link."

attr(conf.int, which="methodname")<-METHOD

return(
list(conf.int=conf.int,
estimate=estimate)  
) 

}



#Negbin.ratio(x=c(0,0,0,1,0), y=c(2,5,3,1,1))
#Negbin.ratio(x=c(0,0,0,0,0), y=c(2,5,3,1,1))
#Negbin.ratio(x=c(1), y=c(2))



#Negbin.ratio(x=c(0,0,0,1,0), y=c(2,5,3,1,1), alternative="less")
#Negbin.ratio(x=c(0,0,0,0,0), y=c(2,5,3,1,1), alternative="less")
#Negbin.ratio(x=c(1), y=c(2), alternative="less")


#Negbin.ratio(x=c(0,0,0,1,0), y=c(2,5,3,1,1), alternative="greater")
#Negbin.ratio(x=c(0,0,0,0,0), y=c(2,5,3,1,1), alternative="greater")
#Negbin.ratio(x=c(1), y=c(2), alternative="greater")

