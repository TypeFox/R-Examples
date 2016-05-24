#source: LASSOProstate.R
#Best model, Cp, AIC, BIC, CV
#
library(lars)
library(ElemStatLearn)
data(prostate)
#center and standardize inputs
Xprostate<-scale(as.matrix.data.frame(prostate[,1:8]))
#center output
lpsa<- prostate$lpsa - mean(prostate$lpsa)
#
#Use lars algorithm to fit LASSO
ansLASSO<-lars(Xprostate,lpsa,type="lasso")
#plot lasso trace plot
ansLASSO<-lars(Xprostate,lpsa,type="lasso")
plot(ansLASSO, xvar="df", lwd=3)
#All possible subset LASSO models
coef(ansLASSO)
#
#Best Cp model
summary(LASSO)
#
#Best AIC model
n <- nrow(Xprostate)
LL<- (-n/2)*log((summary(ansLASSO)$Rss)/n)
AICLASSO<- -2*LL + 2*(0:p)
names(AICLASSO)<-paste("df=",0:p,sep="")
indAICModel<-which.min(AICLASSO)
bLARS<-coef(ansLASSO)[indAICModel,]
names(indAICModel)
#
#Best BIC Model
BICLASSO<- -2*LL + log(n)*(0:p)
names(BICLASSO)<-paste("df=",0:p,sep="")
indBICModel<-which.min(BICLASSO)
bLARS<-coef(ansLASSO)[indBICModel,]
names(indBICModel)
#
#Best CV Model
#cv.lars generates plot for K-fold CV
ans <- cv.lars(Xprostate,lpsa, K=10, type="lasso")
#ans contains components: 'cv', 'cv.error', 'fraction'
#Now add DF to plot and show model selected by one-sd rule
indMin <- which.min(ans$cv)
fMin <- (ans$fraction)[indMin]
cutOff <- (ans$cv.error)[indMin] + (ans$cv)[indMin]
indBest <- min(which((ans$cv)<cutOff)) 
fBest<-(ans$fraction)[indBest]
#tablef - first column is 
TotalAbsBeta<-apply(coef(ansLASSO), MARGIN=1, function(x) sum(abs(x)))
p<-ncol(Xprostate)
tablef<-matrix(c((0:p),TotalAbsBeta/TotalAbsBeta[p]),ncol=2)
#tablef
#
#label the df on the plot
axis(side=3, at=tablef[,2], labels=tablef[,1])
mtext("DF (DF=0 corresponds to intercept only)", side=3, line=2)
abline(v=c(fMin,fBest), col="blue", lwd=3)
abline(h=cutOff, col="red", lwd=2, lty=2)
title(sub="`Best` CV and Smallest CV shown")
pHatF <- round(approx(x=tablef[,2],y=tablef[,1], xout=fBest)$y,1)
text(0.85, ans$cv[5], labels=bquote(hat(p)==.(pHatF))) 
#Note that running this script again will produce a different
# looking graph but using the one-sd rule produces an estimate
# of p close to 2
#
#Best model using 10-fold CV with one-sd rule
bHat<-coef(ansLASSO)[1+round(pHatF),]
bHat[bHat!=0]
