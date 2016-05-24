`BoxCox.Arima` <-
function(object, interval=c(-1,1), type="BoxCox", InitLambda="none", ...){
if (InitLambda=="none")
    initLQ <- FALSE
else
    initLQ <- TRUE
xcall<-deparse(object$call,width.cutoff=500)
#replace x = ? with x=y
ycall<-sub("x = .*, o", "x=y, o", xcall, perl=TRUE)
#get data
test<-gsub("arima(x =","",xcall,fixed=TRUE)
data<-gsub(",.*)","",test)
data<-eval(parse(text=data))
if (initLQ)
    data<-bxcx(data, InitLambda, type=type, InverseQ=TRUE)
if (min(data) <= 0){
    cat(" minimum data value <= 0 so -min+0.25 added to all values", fill=TRUE)
    data<- data + 0.25 - min(data)
    }
#parse object to determine d, ds, s
Gcall<-deparse(object$call,width.cutoff=500)
modelorder<-sub("^.*, order = ", "", Gcall, perl=TRUE)
modelorder<-sub("\\, sea.*", "", modelorder, perl=TRUE )
modelorder<-sub("))",")",modelorder)
modelorder<-eval(parse(text=modelorder))
d<-modelorder[2]
ds<-s<-0
ind<-grep("seasonal", Gcall)
if (length(ind) > 0){
    test<-sub("^.*list","list",Gcall)
    FirstRB<-c(regexpr(")",test))
    SecondRB<-c(regexpr(")",substring(test,FirstRB+1,nchar(test))))
    test<-substring(test,1,FirstRB+SecondRB)
    seaord<-eval(parse(text=test))
    ds<-seaord[[1]][2]
    s<-seaord$period
}
D<-d+s*ds
#now get Jacobian
J <- sum(log(data[(D+1):length(data)]))
# For comparison, this is WRONG: J <- sum(log(data))
#definite loglikelihood function
LogL<-function(lam){
    y<-bxcx(data,lam)
    out.arima<-eval(parse(text=ycall))
    (out.arima$loglik) + (lam-1)*J
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
for (i in 1:length(lams)){
    LL[i]<-exp(LogL(lams[i])-LLlamHat)
    }
ans<-spline(lams,LL)
dataTI<-attr(data,"title")
plot(ans$x, ans$y, type="l", xlab=expression(lambda), ylab=expression("R("*lambda*")"), main="Relative Likelihood Analysis\n95% Confidence Interval",
sub=dataTI)
abline(h=0.1465, col="blue", lwd=2)
text((lamHat+rightLam)/1.8,0.8,bquote(hat(lambda)==.(round(lamHat,3))))
}

