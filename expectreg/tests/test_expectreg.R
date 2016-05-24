require("expectreg")

set.seed(9484470)

###compare 0.5-expectile regression with lm in appropiate cases

ma <- matrix(runif(200)*10,nrow=100)

testfunction <- function(m1,m2) {
    stopifnot(max(abs(fitted(m1)-fitted(m2)))<0.2)
}

#univariate linear model
dat <- data.frame(regressand = 3*ma[,1],covariate.1=ma[,2])
m1 <- expectreg.ls(regressand~ covariate.1,data=dat,smooth="f",expectiles=c(0.5))
#m2 <- expectile.sheets(regressand~covariate.1,data=dat,smooth="f",expectiles=c(0.5))

m3 <- expectreg.ls(regressand~ covariate.1,data=dat,smooth="f",expectiles=c(0.5),estimate="restricted")
m4 <- expectreg.ls(regressand~ covariate.1,data=dat,smooth="f",expectiles=c(0.5),estimate="bundle")

m5 <- expectreg.boost(regressand~ bols(covariate.1),data=dat,expectiles=c(0.5),mstop=112,cv=F)

mcomp <- lm(regressand~covariate.1,data=dat)

testfunction(m1,mcomp)
#testfunction(m2,mcomp)
testfunction(m3,mcomp)
testfunction(m4,mcomp)
testfunction(m5,mcomp)

# univariate model involving sin transformation

dat <- data.frame(regressand = sin(ma[,1]),covariate=ma)
m1 <- expectreg.ls(regressand~rb( covariate.1,"pspline"),data=dat,smooth="schall",expectiles=c(0.5))

#m2 <- expectile.sheets(regressand~rb( covariate.1,"pspline"),data=dat,smooth="acv",expectiles=c(0.5))

m3 <- expectreg.ls(regressand~rb( covariate.1,"pspline"),data=dat,smooth="schall",expectiles=c(0.5),estimate="restricted")
m4 <- expectreg.ls(regressand~rb( covariate.1,"pspline"),data=dat,smooth="schall",expectiles=c(0.5),estimate="bundle")

m5 <- expectreg.boost(regressand~bbs( covariate.1),data=dat,expectiles=c(0.5),mstop=500,cv=F)


mcomp <- lm(regressand~sin(covariate.1),data=dat)

testfunction(m1,mcomp)
#testfunction(m2,mcomp)
testfunction(m3,mcomp)
testfunction(m4,mcomp)
testfunction(m5,mcomp)


# bivariate model: linear and sin
dat <- data.frame(regressand = sin(ma[,1])+3*ma[,2],covariate=ma)

m1 <- expectreg.ls(regressand~rb( covariate.1,"pspline")+rb( covariate.2,"pspline"),data=dat,smooth="schall",expectiles=c(0.5))
#m2 <- expectile.sheets(regressand~rb( covariate.1,"pspline")+base( covariate.2,"pspline"),data=dat,smooth="acv",expectiles=c(0.5))

m3 <- expectreg.ls(regressand~rb( covariate.1,"pspline")+rb( covariate.2,"pspline"),data=dat,smooth="schall",expectiles=c(0.5),estimate="restricted")
m4 <- expectreg.ls(regressand~rb( covariate.1,"pspline")+rb( covariate.2,"pspline"),data=dat,smooth="schall",expectiles=c(0.5),estimate="bundle")

m5 <- expectreg.boost(regressand~bbs( covariate.1)+bbs(covariate.2),data=dat,expectiles=c(0.5),mstop=1200,cv=F)

mcomp <- lm(regressand~sin(covariate.1)+covariate.2,data=dat)

testfunction(m1,mcomp)
#testfunction(m2,mcomp)
testfunction(m3,mcomp)
testfunction(m4,mcomp)
testfunction(m5,mcomp)


data(dutchboys)
dutchb <- dutchboys[sample(1:length(dutchboys[,1]),150),]

expect <- c(0.05,0.5,0.95)


###Values


mLaws <- expectreg.ls(hgt~rb(age,"pspline")+rb(wgt,"pspline"),data=dutchb,smooth="schall",expectiles=expect)
mSheets <- expectreg.ls(hgt~rb(age,"pspline")+rb(wgt,"pspline"),data=dutchb,smooth="f",lambda=1,expectiles=expect,estimate="sheets")
mNoncross <-  expectreg.qp(hgt~rb(age,"pspline"),data=dutchb,smooth="schall",expectiles=expect)
mBundle <- expectreg.ls(hgt~rb(age,"pspline")+rb(wgt,"pspline"),data=dutchb,smooth="schall",expectiles=expect,estimate="bundle")
mRestricted <- expectreg.ls(hgt~rb(age,"pspline")+rb(wgt,"pspline"),data=dutchb,smooth="schall",expectiles=expect,estimate="restricted")
mBoost <- expectreg.boost(hgt~bbs(age)+bbs(wgt),data=dutchb,expectiles=expect,mstop=400)

###INTERCEPTS
stopifnot(length(mLaws$intercepts) == length(mLaws$asymmetries))
#stopifnot(length(mSheets$intercepts) == length(mSheets$asymmetries))
stopifnot(length(mNoncross$intercepts) == length(mNoncross$asymmetries))
stopifnot(length(mBundle$intercepts) == length(mBundle$asymmetries))
stopifnot(length(mRestricted$intercepts) == length(mRestricted$asymmetries))

###COEFFICIENTS A matrix of all the coefficients, for each base element a row and for each expectile a column. 
stopifnot(length(mLaws$coefficients)==2)
#stopifnot(length(mSheets$coefficients)==2)
stopifnot(length(mNoncross$coefficients)==1)
stopifnot(length(mBundle$coefficients)==2)
stopifnot(length(mRestricted$coefficients)==2)

for(i in 1:2) {
    stopifnot(nrow(mLaws$coefficients[[i]])== 20)
    stopifnot(ncol(mLaws$coefficients[[i]])== length(mLaws$asymmetries))
}
#for(i in 1:2) {
#   stopifnot(ncol(mSheets$coefficients[[i]])==length(mSheets$asymmetries))
#}


stopifnot(ncol(mNoncross$coefficients[[1]])== length(mNoncross$asymmetries))


for(i in 1:2) {
    stopifnot(nrow(mBundle$coefficients[[i]])== 20)
    stopifnot(ncol(mBundle$coefficients[[i]])== length(mBundle$asymmetries))
}
for(i in 1:2) {
    stopifnot(nrow(mRestricted$coefficients[i])== 20)
    stopifnot(ncol(mRestricted$coefficients[i])== length(mRestricted$asymmetries))
}

###VALUES

###RESPONSE
stopifnot(dutchb$hgt==mLaws$response)
stopifnot(dutchb$hgt==mSheets$response)
stopifnot(dutchb$hgt==mNoncross$response)
stopifnot(dutchb$hgt==mBundle$response)
stopifnot(dutchb$hgt==mRestricted$response)
stopifnot(dutchb$hgt==mBoost$response)

###COVARIATES
stopifnot(dutchb$age==mLaws$covariates$age && dutchb$wgt==mLaws$covariates$wgt )
stopifnot(dutchb$age==mSheets$covariates$age && dutchb$wgt==mSheets$covariates$wgt )
stopifnot(dutchb$age==mNoncross$covariates$age )
stopifnot(dutchb$age==mBundle$covariates$age && dutchb$wgt==mBundle$covariates$wgt)
stopifnot(dutchb$age==mRestricted$covariates$age && dutchb$wgt==mRestricted$covariates$wgt)

###EXPECTILES
stopifnot(mLaws$asymmetries==expect)
stopifnot(mSheets$asymmetries==expect)
stopifnot(mNoncross$asymmetries==expect)
stopifnot(mBundle$asymmetries==expect)
stopifnot(mRestricted$asymmetries==expect)
stopifnot(mBoost$asymmetries==expect)

###EFFECTS
stopifnot(length(mLaws$covariates)==length(mLaws$effects))
stopifnot(length(mSheets$covariates)==length(mSheets$effects))
stopifnot(length(mNoncross$covariates)==length(mNoncross$effects))
stopifnot(length(mBundle$covariates)==length(mBundle$effects))
stopifnot(length(mRestricted$covariates)==length(mRestricted$effects))
stopifnot(length(mBoost$covariates)==length(mBoost$effects))


#Methods

###predict

stopifnot(isTRUE(all.equal((predict(mLaws,newdata = dutchb))$fitted, fitted(mLaws))))
#stopifnot(isTRUE(all.equal((predict(mSheets,newdata = dutchb))$fitted, fitted(mSheets))))
#stopifnot(isTRUE(all.equal((predict(mNoncross,newdata = dutchb))$fitted, fitted(mNoncross))))

stopifnot(isTRUE(all.equal((predict(mBundle,newdata = dutchb))$fitted, fitted(mBundle))))
stopifnot(isTRUE(all.equal((predict(mRestricted,newdata = dutchb))$fitted, fitted(mRestricted))))

stopifnot(isTRUE(all.equal((predict(mBoost,newdata = dutchb))$fitted, fitted(mBoost))))

###resid
stopifnot(isTRUE(all.equal(resid(mLaws)$fitted, mLaws$response-mLaws$fitted)))
stopifnot(isTRUE(all.equal(resid(mSheets)$fitted, mSheets$response-mSheets$fitted)))
stopifnot(isTRUE(all.equal(resid(mNoncross)$fitted, mNoncross$response-mNoncross$fitted)))

stopifnot(isTRUE(all.equal(resid(mBundle)$fitted, mBundle$response-mBundle$fitted)))
stopifnot(isTRUE(all.equal(resid(mRestricted)$fitted, mRestricted$response-mRestricted$fitted)))

stopifnot(isTRUE(all.equal(resid(mBoost)$fitted, mBoost$response-mBoost$fitted)))
