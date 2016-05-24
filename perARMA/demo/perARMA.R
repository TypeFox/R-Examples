data(volumes)
pgram(t(volumes),length(volumes),datastr='volumes')

############

m=32
win=matrix(1/m,1,m)
dev.set(which=1) 
scoh(t(volumes),m,win,datastr='volumes')
###########
dev.set(which=1)
permest(t(volumes),24, 0.05, NaN,'volumes')

############
dev.set(which=1)
persigest(t(volumes),24,.05, NaN,'volumes')

######### rhoci, rho.constatnt.test, rho.equal.test
peracf(t(volumes),24,seq(1,12),NaN,'volumes')

########
Bcoeff(volumes,24,seq(0,12),NaN,'volumes')

########
Bcoeffa(volumes,24,seq(0,12),NaN,'volumes')

######### nancorr,phth2ab

 PePACF_out<-perpacf(t(volumes),24,12,NaN)   
 ppa=PePACF_out$ppa

 ppfplot(ppa,41,.05,'volumes')
 pc<-ppfcoeffab(ppa,41,1,'volumes')

############ acfpacf.acf, acfpacf.pacf 
pmean<-permest(t(volumes),24, 0.05, NaN,'volumes', pp=0)
xd=pmean$xd
acfpacf(xd,24,24,'volumes') 

##############
perYW(xd,24,2,NaN)

#############  parmafil, R_w_ma
estimators<-perYW(xd,24,2,NaN)
estvec=c(estimators$phi[,1],estimators$phi[,2],estimators$del)
loglikec(estvec,xd,c(24,2,0,0))

##############
z<-parmaresid(xd,0,estimators$del,estimators$phi)
resids=z$resids
#parma_ident(t(resids), 24,NaN, 'residuals_PAR2', outdir='PAR2_ident', details=1) 

##############
dev.set(which=1)
predictperYW(xd,24,2,NaN,2)
dev.set(which=1)
data(volumes.sep)
predseries(volumes.sep,t(volumes),24,2)

