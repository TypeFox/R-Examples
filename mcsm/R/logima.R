logima=function(N=10^3){
# logistic modelling of the Pima dataset
like=function(beda){
  mia=mean(Pima.tr$bm)
  prod(exp(beda[1]+(Pima.tr$bm[Pima.tr$t=="Yes"]-mia)*beda[2]))/
  prod(1+exp(beda[1]-(Pima.tr$bm-mia)*beda[2]))
  }

library(MASS)
bmi=Pima.tr$bm-mean(Pima.tr$bm)
glm(Pima.tr$t~bmi,family=binomial)

#Call:  glm(formula = Pima.tr$t ~ bmi, family = binomial)

#Coefficients:
#(Intercept)          bmi
#    -0.7249       0.1048

#Degrees of Freedom: 199 Total (i.e. Null);  198 Residual
#Null Deviance:      256.4
#Residual Deviance: 240  AIC: 244

# checking for the range
be1=seq(-1.2,-.3,le=100)
be2=seq(.02,.18,le=130)
li=matrix(0,ncol=130,nro=100)
for (i in 1:100) for (j in 1:130) li[i,j]=like(c(be1[i],be2[j]))
image(be1,be2,li)
contour(be1,be2,li,add=T)
# leads to sigma1=0.55 & sigma2=0.2121

sim=cbind(rnorm(N,m=-.72,sd=.55),rnorm(N,m=0.10,sd=.212))
weit=apply(sim,1,like)/(dnorm(sim[,1],m=-.72,sd=.55)*dnorm(sim[,2],m=0.10,sd=.2121))
weit=weit/sum(weit)

#control variates
vari1=(1/(1+exp(-sim[,1]-sim[,2]*bmi)))-sum((Pima.tr$t=="Yes"))/length(Pima.tr$bmi)
vari2=(bmi/(1+exp(-sim[,1]-sim[,2]*bmi)))-sum(bmi[Pima.tr$t=="Yes"])/length(Pima.tr$bmi)

resim=sample(1:N,N,rep=TRUE,pro=weit)
reg=as.vector(lm(sim[resim,1]~t(rbind(vari1[resim],vari2[resim]))-1)$coef)

par(mfrow=c(2,1),mar=c(4,2,2,1))
plot(cumsum(weit*sim[,1])/cumsum(weit),type="l",ylim=c(-.8,-.6),lwd=2,xlab="Iterations")
lines(cumsum(weit*(sim[,1]-reg[1]*vari1-reg[2]*vari2))/cumsum(weit),col="sienna",lwd=2,lty=2)

reg=as.vector(lm(sim[resim,2]~t(rbind(vari1[resim],vari2[resim]))-1)$coef)
plot(cumsum(weit*sim[,2])/cumsum(weit),type="l",ylim=c(.08,.12),lwd=2,xlab="Iterations")
lines(cumsum(weit*(sim[,2]-reg[1]*vari1-reg[2]*vari2))/cumsum(weit),col="sienna",lwd=2,lty=2)
}
