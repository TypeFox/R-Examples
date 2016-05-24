`BoxCox.FitAR` <-
function(object, interval=c(-1,1), type="BoxCox", InitLambda="none", ...){
if (InitLambda=="none")
    initLQ <- FALSE
else
    initLQ <- TRUE
if (object$ARModel=="ARp") 
    ycall<-"GetFitARpLS(z=y, p=pBXCX)"
if (object$ARModel=="ARz") 
    ycall<-"GetFitARz(z=y, p=pBXCX)"
zBXCX<-object$z
pBXCX <- object$pvec
if (initLQ)
    zBXCX<-bxcx(zBXCX, InitLambda, type=type, InverseQ=TRUE)
if (min(zBXCX) <= 0){
    cat(" minimum data value <= 0 so -min+0.25 added to all values", fill=TRUE)
    data<- data + 0.25 - min(data)
    }
#now get Jacobian
J <- sum(log(zBXCX))
#definite loglikelihood function
LogL<-function(lam){
    y<-bxcx(zBXCX,lam)
    y<-y-mean(y)
    out<-eval(parse(text=ycall))
    (out$loglikelihood) + (lam-1)*J
    }
#optimize
ans<-optimize(LogL,interval=interval,maximum=TRUE)
lamHat<-ans$maximum
LLlamHat<-ans$objective
#determine right and left ends for plotting RL
RL<-1
rightLam<-lamHat
MaxIt<-10
iter<-0
while (RL>0.01 && iter<MaxIt){
    rightLam<-rightLam+0.1
    RL<-exp(LogL(rightLam)-LLlamHat)
    iter<-iter+1
    }
RL<-1
leftLam<-lamHat
MaxIt<-10
iter<-0
while (RL>0.01 && iter<MaxIt){
    leftLam<-leftLam-0.1
    RL<-exp(LogL(leftLam)-LLlamHat)
    iter<-iter+1
    }
lams<-seq(leftLam,rightLam,length=21)
LL<-numeric(length(lams))
for (i in 1:length(lams)){
    LL[i]<-exp(LogL(lams[i])-LLlamHat)
    }
ans<-spline(lams,LL)
dataTI<-attr(data,"title")
plot(ans$x, ans$y, type="l", xlab=expression(lambda), ylab=expression("R("*lambda*")"), main="Relative Likelihood Analysis\n95% Confidence Interval",
sub=dataTI)
abline(h=0.1465, col="blue", lwd=2)
text((lamHat+rightLam)/1.85,0.8,bquote(hat(lambda)==.(round(lamHat,3))))
}

