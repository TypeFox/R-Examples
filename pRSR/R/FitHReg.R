FitHReg <-
function (y, t = 1:length(y),nf=150) 
{
    theta <- numeric(2)
    theta[1] <- length(y)
	theta[2] <- nf #number of frequencies
    ansH <- .C("GetHReg", y = as.double(y), t = as.double(t), 
        theta = as.double(theta), PACKAGE = "pRSR")
    LR <- (ansH$theta)[1]
    f <- (ansH$theta)[2]
    fault <- (ansH$theta)[3]
    x1<-cos(2*pi*f*t)
    x2<-sin(2*pi*f*t)
    ansLM<-lm(y~x1+x2)
    co<-coef(ansLM)
    names(co)<-c("mu","A","B")
    a<-summary(ansLM)
    v<-a$cov.unscaled
    Rsq<-a$r.squared
    fstat<- a$fstatistic
    sig<-a$sigma
    ans<-list(coefficients=co, residuals=residuals(ansLM), Rsq=Rsq, fstatistic=fstat, sigma=sig, freq=f, LRStat=LR, fault=fault)
    class(ans)<-"HReg"
    ans
}

