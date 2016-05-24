`BoxCox.numeric` <-
function(object, interval=c(-1,1), IIDQ = FALSE, ...){
n<-length(object)
if (n <= 20 || IIDQ)
    p<-0
else 
    p<-SelectModel(object, lag.max=min(30,round(length(object)/4)), Best=1)
z<-as.vector(object)
if (min(z) <= 0){
    cat(" minimum data value <= 0 so -min+0.25 added to all values", fill=TRUE)
    z<- z + 0.25 - min(z)
    }
#now get Jacobian
J <- sum(log(z))
#definite loglikelihood function
LogL<-function(lam){
    y<-bxcx(z,lam)
    y<-y-mean(y)
    out<-GetFitAR(y,p)
    (out$loglikelihood) + (lam-1)*J
    }
#optimize
ans<-optimize(LogL,interval=interval,maximum=TRUE)
lamHat<-ans$maximum
LLlamHat<-ans$objective
#determine right and left ends for plotting RL
RL<-1
rightLam<-lamHat
while (RL>0.01){
    rightLam<-rightLam+0.1
    RL<-exp(LogL(rightLam)-LLlamHat)
    }
RL<-1
leftLam<-lamHat
while (RL>0.01){
    leftLam<-leftLam-0.1
    RL<-exp(LogL(leftLam)-LLlamHat)
    }
lams<-seq(leftLam,rightLam,length=21)
LL<-numeric(length(lams))
for (i in 1:length(lams))
    LL[i]<-exp(LogL(lams[i])-LLlamHat)
ans<-spline(lams,LL)
dataTI<-attr(data,"title")
plot(ans$x, ans$y, type="l", xlab=expression(lambda), ylab=expression("R("*lambda*")"), main="Relative Likelihood Analysis\n95% Confidence Interval",
sub=dataTI)
abline(h=0.1465, col="blue", lwd=2)
text((lamHat+rightLam)/1.85,0.8,bquote(hat(lambda)==.(round(lamHat,3))))
}

