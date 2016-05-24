dotdotdot <- function(data, indices, b,pi){
  x = data[indices]
  y = mean(x) + b - pi;
  temp = b 
  return(y)
}


test.dotdotdot <-function(){
 set.seed(27)
 b = boot(discoveries, dotdotdot, R=100, b=33,pi=63)
 set.seed(27)
	a = pboot(discoveries, dotdotdot, R=100, b=33,pi=63)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
 checkEquals(a,b,"Two extra variables passed via ...")
}

ddd <- function(data, indices, mydd){
  x = data[indices]
  return(mean(x) + mean(mydd))
}

test.ddd <- function(){
#  DEACTIVATED()
 set.seed(27)
 a = boot(discoveries, ddd, R=100, mydd = 33)
 set.seed(27)
	b = pboot(discoveries, ddd, R=100, mydd = 33)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
 checkEquals(a,b,"interger passed via ...")
 
 set.seed(27)
 a = boot(discoveries, ddd, R=100, mydd = discoveries) 
 set.seed(27)
	b = pboot(discoveries, ddd, R=100, mydd = discoveries)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
 checkEquals(a,b,"vector passed via ...")
 
 
 set.seed(27)
 a = boot(discoveries, ddd, R=200, mydd = trees$Girth)
 set.seed(27)
	b = pboot(discoveries, ddd, R=200, mydd = trees$Girth)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
 checkEquals(a,b,"label passed via ...")
 
 set.seed(27)
 a = boot(discoveries, ddd, R=100, mydd = seq(50))
 set.seed(27)
	b = pboot(discoveries, ddd, R=100, mydd = seq(50))
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
 checkEquals(a,b,"expression passed via ...")
 
 set.seed(27)
 a = boot(discoveries, ddd, R=100, mydd = trees[,1])
 set.seed(27)
	b = pboot(discoveries, ddd, R=100, mydd = trees[,1])
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
 checkEquals(a,b,"dataframe row passed via ...")
}

ddddataframe <- function(data, indices, b){
  x = data[indices]
  return(x + mean(b[,1]))
}

test.ddddataframe <- function(){
 set.seed(27)
 a = boot(discoveries, ddddataframe, R=100, b=trees)
 set.seed(27)
	b = pboot(discoveries, ddddataframe, R=100, b=trees)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
 checkEquals(a,b,"dataframe passed via ...")
}

#  taken from the boot library documentation
#  In this example we show the use of boot in a prediction from 
#  regression based on the nuclear data.  This example is taken 
#  from Example 6.8 of Davison and Hinkley (1997).  Notice also 
#  that two extra arguments to statistic are passed through boot.
nuke <- nuclear[,c(1,2,5,7,8,10,11)]
nuke.lm <- glm(log(cost)~date+log(cap)+ne+ ct+log(cum.n)+pt, data=nuke)
nuke.diag <- glm.diag(nuke.lm)
nuke.res <- nuke.diag$res*nuke.diag$sd
nuke.res <- nuke.res-mean(nuke.res)


#  We set up a new data frame with the data, the standardized 
#  residuals and the fitted values for use in the bootstrap.
nuke.data <- data.frame(nuke,resid=nuke.res,fit=fitted(nuke.lm))


#  Now we want a prediction of plant number 32 but at date 73.00
new.data <- data.frame(cost=1, date=73.00, cap=886, ne=0,
                       ct=0, cum.n=11, pt=1)
new.fit <- predict(nuke.lm, new.data)


nuke.fun <- function(dat, inds, i.pred, fit.pred, x.pred)
{
     assign(".inds", inds, envir=.GlobalEnv)
     lm.b <- glm(fit+resid[.inds] ~date+log(cap)+ne+ct+
          log(cum.n)+pt, data=dat)
     pred.b <- predict(lm.b,x.pred)
     remove(".inds", envir=.GlobalEnv)
     c(coef(lm.b), pred.b-(fit.pred+dat$resid[i.pred]))
}

test.nuke <- function(){
  set.seed(27) 
  a <- boot(nuke.data, nuke.fun, R=999, m=1, fit.pred=new.fit, x.pred=new.data)
  set.seed(27) 
	b <- boot(nuke.data, nuke.fun, R=999, m=1, fit.pred=new.fit, x.pred=new.data)
# Ignore the calls having different names when testing equality.
	b$call <- a$call 
  checkEquals(a,b,"nuke test from the boot man page")
}
