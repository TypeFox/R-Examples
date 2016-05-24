plotdhglm <-
function (OUTPUT,type="mean",random=NULL) {
    type="mean"
    random=NULL
    par(mfrow=c(2,2))
    mu<-OUTPUT[2][[1]]
    StudentResidual<-OUTPUT[1][[1]]
    x<-mu
    y<-StudentResidual
    fit<- smooth.spline(x,y,cv=TRUE)
    plot(x, y, main="Residuals vs Fitted", xlab="mu", ylab="StudentizedResidual", cex=0.5) #plot data point
    lines(fit$x, fit$y) #plot smooth spline fit
    y<-abs(StudentResidual)
    fit<- smooth.spline(x,y,cv=TRUE)
    plot(x, y, main="|Residuals| vs Fitted",xlab="mu", ylab="|StudentizedResidual|", cex=0.5) #plot data point
    lines(fit$x, fit$y) #plot smooth spline fit
    qqnorm(StudentResidual); qqline(StudentResidual) # Q-Q plot
    hist(StudentResidual)
}
