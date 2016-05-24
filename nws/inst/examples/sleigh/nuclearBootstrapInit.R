library(boot)
data(nuclear)

# take column 1, 2, 5, 7, 8, 10, and 11 of nuclear data frame.
# column 1 -> cost, column 2 -> date
# column 5 -> cap,  column 7 -> ne
# column 8 -> ct,   column 10 -> cum.n
# column 11 -> pt 
nuke<-nuclear[,c(1,2,5,7,8,10,11)]

# glm is used to fit generalized linear models
# formula for this problem set is "log(cost)~date+log(cap)+ne+ct+log(cum.n)+pt
# data is from nuke, the data frame we constructed earlier
nuke.lm<-glm(log(cost)~date+log(cap)+ne+ct+log(cum.n)+pt, data=nuke)

# more generalized linear model diagnostics
nuke.diag<-glm.diag(nuke.lm)   
nuke.res<-nuke.diag$res*nuke.diag$sd
nuke.res<-nuke.res-mean(nuke.res)

# construct a data frame based on the diagnostics
nuke.data<-data.frame(nuke,resid=nuke.res,fit=fitted(nuke.lm))

# create a data frame with the following variables: 
# cost=1, date=73.00, cap=886, ne=0, ct=0, cum.n=11, pt=1
# these values are used in nuke.fun defined later  
new.data<-data.frame(cost=1,date=73.00, cap=886, ne=0, ct=0, cum.n=11, pt=1)

# nuke.lm is the model object for which prediction is desired.
# new.data contain variables used in the nuke.lm model 
new.fit<-predict(nuke.lm, new.data)


nuke.fun<-function(dat, inds, i.pred, fit.pred, x.pred) {
	assign(".inds", inds, envir=.GlobalEnv)
	lm.b<-glm(fit+resid[.inds] ~date+log(cap)+ne+ct+log(cum.n)+pt, data=dat)
	pred.b <-predict(lm.b, x.pred)
	remove(".inds", envir=.GlobalEnv)
	c(coef(lm.b), pred.b-(fit.pred+dat$resid[i.pred]))
}
