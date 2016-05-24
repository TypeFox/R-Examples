### Testing the additional predictive value of high-dimensional data
###
### Copyright 2009-09 Anne-Laure Boulesteix 
###
### 
###
###
### This file is part of the `globalboosttest' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

globalboosttest<-function(X,Y,Z=NULL,nperm=1000,mstop=1000,mstopAIC=FALSE,pvalueonly=TRUE,plot=FALSE,...)
{
n<-nrow(X)

X<-as.matrix(X)

if (is.null(Z))
 {
 Z<-rep(1,n)
 }

if (!is.Surv(Y)&!is.factor(Y)&!is.numeric(Y))
 {
 stop("Y must be a factor, a numeric or a Surv object!")
 }

#  surv
if (is.Surv(Y))
 {
 datasurv<-data.frame(Y=Y,Z)
 linpred<-predict(coxph(Y~.,data=datasurv))
 loss<-CoxPH()
 if (mstopAIC==TRUE)
  {
  warning("AIC-optimized mstop can be applied for binary Y only")
  mstopAIC<-FALSE
  }
 }

if (is.numeric(Y)&!is.Surv(Y))
 {
 linpred <- predict(lm(Y~.,data=data.frame(Y=Y,Z)))
 loss<-GaussReg()
 }


#  factor
if (is.factor(Y))
 {
 linpred <- predict(glm(Y~.,data=data.frame(Y=Y,Z),family=binomial))
 loss <- Binomial()
 if ( sum(diag(table(Y,sign(linpred)))) == length(Y) ) {
    warning("'Perfect separation' case with clinical variables.")
    glmboostreal <- glmboost(X,Y,family=loss,offset=linpred, control = boost_control(mstop = max(mstop)
                             , nu = 0.1, risk = "inbag"), center=TRUE)
    riskreal<-glmboostreal$risk()
    return(list(riskreal=riskreal,riskperm=NA,mstopAIC=NA,pvalue=NA))
    }
 }


# real
glmboostreal<-glmboost(X,Y,family=loss,offset=linpred, control = boost_control(mstop = max(mstop), nu = 0.1, risk = "inbag"), center=TRUE)
riskreal<-glmboostreal$risk()


# permuted
riskperm<-sapply(as.list(1:nperm),FUN=perm,XX=X,Y=Y,n=n,loss=loss,offset=linpred,mstop=max(mstop))


if (plot==TRUE)
 {
 plot(1:max(mstop),riskreal,type="n",ylim=c(min(riskreal,riskperm),max(riskreal,riskperm)),main="Risk",xlab="mstop",ylab="Risk",col=1,...)
 for (i in 1:nperm) 
  {
  points(1:max(mstop),riskperm[,i],type="l",lty=3,col=8)
  }
 points(1:max(mstop),riskreal,type="l",col=1)
 }

pvalue<-c()

for (m in mstop)
 {
 pvalue<-matrix(c(pvalue,sum(riskperm[m,]<=riskreal[m])/nperm),nrow=1)
 }
colnames(pvalue)<-paste("mstop=",mstop,sep="")
if (mstopAIC==TRUE)
 {
 mstopaic<-mstop(AIC(glmboostreal,method="classical"))
 pvalue<-cbind(pvalue,sum(riskperm[mstopaic,]<=riskreal[mstopaic])/nperm)
 colnames(pvalue)[length(mstop)+1]<-"AIC"
 }

if (pvalueonly)
 {
 list(pvalue=pvalue)
 }

else
 {
 if (mstopAIC==TRUE)
  {
  list(riskreal=riskreal,riskperm=riskperm,mstopAIC=mstopaic,pvalue=pvalue)
  }
 else
  {
  list(riskreal=riskreal,riskperm=riskperm,pvalue=pvalue)
  }
 }

}



#####

perm<-function(i,XX,Y,n,loss,offset,mstop)
{
set.seed(i) 
print(paste("perm=",i,sep=""))
sampi<-sample(n) 
return(glmboost(XX[sampi,],Y,family=loss,offset=offset,control = boost_control(mstop = max(mstop), nu = 0.1, risk = "inbag"), center=TRUE)$risk())
}