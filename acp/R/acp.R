.poissonlik<-function(theta,y,x){ 
n<-nrow(x)
k<-ncol(x)
beta<-theta[1:k]
mu<-exp(x%*%beta)
logl<-(y*(x%*%beta))-mu-log(factorial(y))
sum(-logl)
}

.poissonEst <- function(x, y)
{
k<-ncol(x)
## perform optimization
op<-optim(array(0.1,dim=c(1,k)),.poissonlik,y=y,x=x,method="BFGS",hessian=T)
## compute coefficients and logl
coef<-c(op$par)
logl<-c(op$value)
## degrees of freedom and standard deciation of residuals 
df <- nrow(x)-ncol(x)
sigma2 <- sum((y - exp(x%*%coef))^2)/df
## compute covariance matrix
vcov<-solve(op$hessian)
colnames(vcov) <- rownames(vcov) <- colnames(x)
BIC = -2*(-logl/nrow(x))+(k+1)*log(nrow(x))/nrow(x)
AIC = -2*(-logl/nrow(x))+2*(k+1)/nrow(x)
list(coefficients = coef, logl = -logl, AIC = AIC, BIC = BIC,
vcov = vcov,
sigma = sqrt(sigma2),
df = df)
}

.poisson.default <- function(x, y)
{
x <- as.matrix(x)
y <- as.numeric(y)
est <- .poissonEst(x, y)
est$fitted.values <- as.vector(exp(x %*% est$coefficients))
est$residuals <- y - est$fitted.values
est$logl <- est$logl
est$BIC <- est$BIC
est$call <- match.call()
class(est) <- "acp"
est
}


.acpliknc<-function(theta, y, pp, qq){ 
n<-nrow(y)
lambda <- matrix(NA,n,1)
lambda[1] <- mean(y)
vylag = Lag(y,1)
vparlag = c(theta[2])
if(pp > 1){
for (i in 1:(pp-1)){ vparlag <- append(vparlag, theta[2+i]) 
vylag <- cbind(vylag,Lag(y,i+1))
lambda[i+1] <- mean(y)
}}
vfilt = c(theta[2+pp])
vinit = c(lambda[1])
if(qq > 1){
for (i in 1:(qq-1)){ vfilt <- append(vfilt, theta[2+pp+i]) 
vinit <- append(vinit,lambda[1])
}}
lambda[(pp+1):n] <- filter(theta[1] + na.omit(vylag%*%vparlag), filter = vfilt, init = vinit, method = "recursive")
mu<-lambda
logl<-(y*log(mu))-mu-lfactorial(y)
sum(-logl)
}

.acplikc<-function(theta, y, x, pp, qq){ 
n<-nrow(x)
k<-ncol(x)
beta<-theta[1:k]
lambda <- matrix(NA,n,1)
lambda[1] <- mean(y)
vylag = Lag(y,1)
vparlag = c(theta[k+2])
if(pp > 1){
for (i in 1:(pp-1)){ vparlag <- append(vparlag, theta[k+2+i]) 
vylag <- cbind(vylag,Lag(y,i+1))
lambda[i+1] <- mean(y)
}}
vfilt = c(theta[k+2+pp])
vinit = c(lambda[1])
if(qq > 1){
for (i in 1:(qq-1)){ vfilt <- append(vfilt, theta[k+2+pp+i]) 
vinit <- append(vinit,lambda[1])
}}
lambda[(pp+1):n] <- filter(theta[k + 1] + na.omit(vylag%*%vparlag), filter = vfilt, init = vinit, method = "recursive")
mu<-exp(x%*%beta)*lambda
logl<-(y*log(mu))-mu-lfactorial(y)
sum(-logl)
}


.acpEst <- function(x, y, p, q, startval, varopt)
{
n<-nrow(x)
if(sum(x)==0){
if(is.null(startval)){
#arma regression to detect initial values for autoregressive parameters
ar<-arma(y, order = c(p, q))
asval <- matrix(NA,1,1+p+q)
for (i in 1:p){ asval[1+i] = abs(coef(ar)[i]) }
for (i in 1:q){ asval[1+p+i] = abs(coef(ar)[p+i]) }
asval[1]<-(mean(y)*abs(1-asval[2]-asval[1+p+1]))
## perform optimization
op<-optim(c(asval),fn=.acpliknc,gr=NULL,y=y,pp=p,qq=q,method="BFGS",hessian=varopt)
}
else{
op<-optim(startval,fn=.acpliknc,gr=NULL,y=y,pp=p,qq=q,method="BFGS",hessian=varopt)
}

## compute coefficients and logl
coef<-op$par
logl<-op$value
## compute conditional mean
lambda <- matrix(NA,n,1)
lambda[1] <- mean(y)
vylag = Lag(y,1)
vparlag = c(coef[2])
if(p > 1){
for (i in 1:(p-1)){ vparlag <- append(vparlag, coef[2+i]) 
vylag <- cbind(vylag,Lag(y,i+1))
lambda[i+1] <- mean(y)
}}
vfilt = c(coef[2+p])
vinit = c(lambda[1])
if(q > 1){
for (i in 1:(q-1)){ vfilt <- append(vfilt, coef[2+p+i]) 
vinit <- append(vinit,lambda[1])
}}
lambda[(p+1):n] <- filter(coef[1] + na.omit(vylag%*%vparlag), filter = vfilt, init = vinit, method = "recursive")
mu<-lambda
## compute standardized residuals
pres<-(y-mu)/sqrt(mu)
## degrees of freedom and standard deviation of residuals 
df <- nrow(y)-p-q-1
sigma2 <- sum(pres^2)/df
if(is.null(op$hessian))
{
vcov<- matrix(NA,1,p+q+1)
## create vector of parameter names
namevec <- matrix(NA,1,p+q+1)
namevec[1]<-"a"
for (i in 1:p){ 
namevec[1+i]<-paste("b",i)
}
for (i in 1:q){ 
namevec[1+p+i]<-paste("c",i)
}
BIC = -2*(-logl/n)+(p+q+1)*log(n)/n
AIC = -2*(-logl/n)+2*(p+q+1)/n
list(coefficients = coef, logl = -logl, vcov = vcov, BIC = BIC, AIC = AIC,
df = df)
}
else
{
## compute covariance matrix
vcov<-solve(op$hessian)
## create vector of parameter names
namevec <- matrix(NA,1,p+q+1)
namevec[1]<-"a"
for (i in 1:p){ 
namevec[1+i]<-paste("b",i)
}
for (i in 1:q){ 
namevec[1+p+i]<-paste("c",i)
}
colnames(vcov) <- rownames(vcov) <- namevec
BIC = -2*(-logl/n)+(p+q+1)*log(n)/n
AIC = -2*(-logl/n)+2*(p+q+1)/n
list(coefficients = coef, logl = -logl, BIC = BIC, AIC = AIC,
vcov = vcov,
sigma = sqrt(sigma2),
df = df, p = p, q = q)
}
}
else{
k<-ncol(x)
if(is.null(startval)){
## poisson regression to detect initial values for covariates
pr<-optim(array(0.1,dim=c(1,k)),.poissonlik,y=y,x=x,method="BFGS",hessian=F)
#arma regression to detect initial values for autoregressive parameters
ar<-arma(y, order = c(p, q))
asval <- matrix(NA,1,1+p+q)
for (i in 1:p){ asval[1+i] = abs(coef(ar)[i]) }
for (i in 1:q){ asval[1+p+i] = abs(coef(ar)[p+i]) }
asval[1]<-(mean(y)*abs(1-asval[2]-asval[1+p+1]))
## perform optimization
op<-optim(c(pr$par,asval),fn=.acplikc,gr=NULL,y=y,x=x,pp=p,qq=q,method="BFGS",hessian=varopt)
}
else{
op<-optim(startval,fn=.acplikc,gr=NULL,y=y,x=x,pp=p,qq=q,method="BFGS",hessian=varopt)
}
## compute coefficients and logl
coef<-op$par
logl<-op$value
## compute conditional mean
lambda <- matrix(NA,n,1)
lambda[1] <- mean(y)
vylag = Lag(y,1)
vparlag = c(coef[k+2])
if(p > 1){
for (i in 1:(p-1)){ vparlag <- append(vparlag, coef[k+2+i]) 
vylag <- cbind(vylag,Lag(y,i+1))
lambda[i+1] <- mean(y)
}}
vfilt = c(coef[k+2+p])
vinit = c(lambda[1])
if(q > 1){
for (i in 1:(q-1)){ vfilt <- append(vfilt, coef[k+2+p+i]) 
vinit <- append(vinit,lambda[1])
}}
lambda[(p+1):n] <- filter(coef[k+1] + na.omit(vylag%*%vparlag), filter = vfilt, init = vinit, method = "recursive")
mu<-exp(x%*%coef[1:k])*lambda
## compute standardized residuals
pres<-(y-mu)/sqrt(mu)
## degrees of freedom and standard deviation of residuals 
df <- nrow(y)-k-p-q-1
sigma2 <- sum(pres^2)/df
if(is.null(op$hessian))
{
vcov<- matrix(NA,1,k+p+q+1)
## create vector of parameter names
namevec <- matrix(NA,1,k+p+q+1)
namevec[,1:k]<-colnames(x)
namevec[k+1]<-"a"
for (i in 1:p){ 
namevec[k+1+i]<-paste("b",i)
}
for (i in 1:q){ 
namevec[k+1+p+i]<-paste("c",i)
}
BIC = -2*(-logl/n)+(k+p+q+1)*log(n)/n
AIC = -2*(-logl/n)+2*(k+p+q+1)/n
list(coefficients = coef, logl = -logl, BIC = BIC, AIC = AIC,
vcov = vcov,
sigma = sqrt(sigma2),
df = df, p = p, q = q)
}
else
{
## compute covariance matrix
vcov<-solve(op$hessian)
## create vector of parameter names
namevec <- matrix(NA,1,k+p+q+1)
namevec[,1:k]<-colnames(x)
namevec[k+1]<-"a"
for (i in 1:p){ 
namevec[k+1+i]<-paste("b",i)
}
for (i in 1:q){ 
namevec[k+1+p+i]<-paste("c",i)
}
colnames(vcov) <- rownames(vcov) <- namevec
BIC = -2*(-logl/n)+(k+p+q+1)*log(n)/n
AIC = -2*(-logl/n)+2*(k+p+q+1)/n
list(coefficients = coef, logl = -logl, BIC = BIC, AIC = AIC,
vcov = vcov,
sigma = sqrt(sigma2),
df = df, p = p, q = q)
}
}
}

acp <- function(x, ...) UseMethod("acp")


acp.default <- function(x, y, p, q, startval, varopt,...)
{
n<-nrow(x)
x <- as.matrix(x)
y <- as.matrix(y)
if(sum(x)==0){
n<-nrow(x)
est <- .acpEst(x, y, p, q ,startval, varopt)
lambda <- matrix(NA,n,1)
lambda[1] <- mean(y)
vylag = Lag(y,1)
vparlag = c(est$coefficients[2])
if(p > 1){
for (i in 1:(p-1)){ vparlag <- append(vparlag, est$coefficients[2+i]) 
vylag <- cbind(vylag,Lag(y,i+1))
lambda[i+1] <- mean(y)
}}
vfilt = c(est$coefficients[2+p])
vinit = c(lambda[1])
if(q > 1){
for (i in 1:(q-1)){ vfilt <- append(vfilt, est$coefficients[2+p+i]) 
vinit <- append(vinit,lambda[1])
}}
lambda[(p+1):n] <- filter(est$coefficients[1] + na.omit(vylag%*%vparlag), filter = vfilt, init = vinit, method = "recursive")
mu<-lambda
est$fitted.values <- as.vector(mu)
est$residuals <- (y - est$fitted.values)/sqrt(est$fitted.values)
est$logl <- est$logl
est$AIC <- est$AIC
est$BIC <- est$BIC 
}
else{
k<-ncol(x)
est <- .acpEst(x, y, p, q, startval, varopt)
lambda <- matrix(NA,n,1)
lambda[1] <- mean(y)
vylag = Lag(y,1)
vparlag = c(est$coefficients[k+2])
if(p > 1){
for (i in 1:(p-1)){ vparlag <- append(vparlag, est$coefficients[k+2+i]) 
vylag <- cbind(vylag,Lag(y,i+1))
lambda[i+1] <- mean(y)
}}
vfilt = c(est$coefficients[k+2+p])
vinit = c(lambda[1])
if(q > 1){
for (i in 1:(q-1)){ vfilt <- append(vfilt, est$coefficients[k+2+p+i]) 
vinit <- append(vinit,lambda[1])
}}
lambda[(p+1):n] <- filter(est$coefficients[k+1] + na.omit(vylag%*%vparlag), filter = vfilt, init = vinit, method = "recursive")
mu<-exp(x%*%est$coefficients[1:k])*lambda
est$fitted.values <- as.vector(mu)
est$residuals <- (y - est$fitted.values)/sqrt(est$fitted.values)
est$logl <- est$logl
est$AIC <- est$AIC
est$BIC <- est$BIC 
}
est$call <- match.call()
class(est) <- "acp"
est
}

print.acp <- function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients)
cat("\nLogl:\n")
print(x$logl)
cat("\nAIC:\n")
print(x$AIC)
cat("\nBIC:\n")
print(x$BIC)
}


summary.acp <- function(object, ...)
{
if(is.na(sum(object$vcov)))
{
cat("Call:\n")
print(object$call)
cat("\nCoefficients:\n")
print(object$coefficients)
cat("\nLogl:\n")
print(object$logl)
cat("\nAIC:\n")
print(object$AIC)
cat("\nBIC:\n")
print(object$BIC)
}
else
{
se <- sqrt(diag(object$vcov))
tval <- coef(object) / se
TAB <- cbind(Estimate = coef(object),
StdErr = se,
t.value = tval,
p.value = 2*pt(-abs(tval), df=object$df))
res <- list(call=object$call,
coefficients=TAB, logl=object$logl, AIC=object$AIC ,BIC=object$BIC)
class(res) <- "summary.acp"
res
}
}

print.summary.acp <- function(x, ...)
{
if(is.na(sum(x$vcov)))
{
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients)
cat("\nLogl:\n")
print(x$logl)
cat("\nAIC:\n")
print(x$AIC)
cat("\nBIC:\n")
print(x$BIC)
}
else
{
cat("Call:\n")
print(x$call)
cat("\n")
printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
cat("\nLogl:\n")
print(x$logl)
cat("\nAIC:\n")
print(x$AIC)
cat("\nBIC:\n")
print(x$BIC)
}
}


acp.formula <- function(formula, data=list(), p = NULL ,q = NULL ,startval=NULL, varopt=T, family="acp",...)
{
mf <- model.frame(formula=formula, data=data)
x <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)
if (family=="acp"){
est <- acp.default(x, y, p, q, startval, varopt,...)
}
else{
est <- .poisson.default(x, y, ...)
}
est$call <- match.call()
est$formula <- formula
est$data <- data
est$family <- family
est
}

predict.acp <- function(object, newydata=NULL , newxdata=NULL,...)
{
if(is.null(newydata)){y <- fitted(object)}
n<-length(newydata)
k<-ifelse(is.null(newxdata),0,ncol(newxdata))
if(k>0){x<-newxdata}
if (object$family=="acp"){
q = object$q
p = object$p
lambda <- matrix(NA,n,1)
lambda[1] <- 1
vylag = Lag(newydata,1)
vparlag = c(coef(object)[k+2])
if(p > 1){
for (i in 1:(p-1)){ vparlag <- append(vparlag, coef(object)[k+2+i]) 
vylag <- cbind(vylag,Lag(newydata,i+1))
lambda[i+1] <- 1
}}
vfilt = c(coef(object)[k+2+p])
vinit = c(lambda[1])
if(q > 1){
for (i in 1:(q-1)){ vfilt <- append(vfilt, coef(object)[k+2+p+i]) 
vinit <- append(vinit,lambda[1])
}}
lambda[(p+1):n] <- filter(coef(object)[k+1] +  na.omit(vylag%*%vparlag), filter = vfilt, init = vinit, method = "recursive")
if (k==0) {
   xb=matrix(1,n,1)
} else {
   xb=exp(as.matrix(x)%*%coef(object)[1:k])
}
y <- as.vector(xb*lambda)
}
else{
y <- as.vector(exp(as.matrix(x)%*%coef(object)[1:k]))
}
y
}


### function for nonrandomized PIT histogram 
###
### input: 
###   x    observed data 
###   Px   CDF at x 
###   Px1  CDF at x-1 

.pit <- function(x, Px, Px1, n.bins=10, y.max=2.75, my.title="PIT Histogram")
  {
  a.mat <- matrix(0,n.bins,length(x))
  k.vec <- pmax(ceiling(n.bins*Px1),1)
  m.vec <- ceiling(n.bins*Px)
  d.vec <- Px-Px1
  for (i in 1:length(x))
      {
      if (k.vec[i]==m.vec[i]) {a.mat[k.vec[i],i]=1}
      else 
        { 
        a.mat[k.vec[i],i]=((k.vec[i]/n.bins)-Px1[i])/d.vec[i]
        if ((k.vec[i]+1)<=(m.vec[i]-1))
           {for (j in ((k.vec[i]+1):(m.vec[i]-1))) {a.mat[j,i]=(1/(n.bins*d.vec[i]))}}
        a.mat[m.vec[i],i]=(Px[i]-((m.vec[i]-1)/n.bins))/d.vec[i]     
        }
      }
  a <- apply(a.mat,1,sum)
  a <- (n.bins*a)/(length(x))
  p <- (0:n.bins)/n.bins
  PIT <- "Probability Integral Transform"
  RF <- "Relative Frequency"
  plot(p, p, ylim=c(0,y.max), type="n", xlab=PIT, ylab=RF, main=my.title) 
  temp1 <- ((1:n.bins)-1)/n.bins
  temp2 <- ((1:n.bins)/n.bins)
  o.vec <- rep(0,n.bins)
  segments(temp1,o.vec,temp1,a)
  segments(temp1,a,temp2,a)
  segments(temp2,o.vec,temp2,a)
  segments(0,0,1,0)
  }
  
evaluation <- function(ydata, yhatdata,...)
{

## compute Pearson standardized residuals
pres<-(ydata-yhatdata)/sqrt(yhatdata)

### parameter settings for computing scores

kk <- 100000                            ### cut-off for summations 
my.k <- (0:kk) - 1                      ### to handle ranked probability score

lambdahat<-yhatdata

pois.Px <- ppois(ydata,lambdahat)                        ### cumulative probabilities
pois.Px1 <- ppois(ydata-1,lambdahat)
pois.px <- dpois(ydata,lambdahat)                        ### probabilities 


pois.logs <- - log(pois.px)
pois.norm <- sum(dpois(my.k,lambdahat)^2) 
pois.qs <- - 2*pois.px + exp(-2*lambdahat)*besselI(2*lambdahat,0)
pois.sphs <- - pois.px / sqrt(exp(-2*lambdahat)*besselI(2*lambdahat,0))
i.cumsum <- cumsum(ppois(my.k,lambdahat)^2)
ii.sum <- sum((ppois(my.k,lambdahat)-1)^2)
ii.cumsum <- cumsum((ppois(my.k,lambdahat)-1)^2)
pois.rps <- (i.cumsum[ydata+1] + ii.sum - ii.cumsum[ydata+1]) 
pois.dss <- (ydata-lambdahat)^2/lambdahat + log(lambdahat)
pois.ses <- (ydata-lambdahat)^2
pois.mae <- abs(ydata-lambdahat)

scores <- matrix(c(round(mean(pois.logs),2),round(mean(pois.qs),3),round(mean(pois.sphs),3),round(mean(pois.rps),2),round(mean(pois.dss),2),round(mean(pois.ses),1),round(mean(pois.mae),1),round(mean(pois.ses)^0.5,1)),ncol=1,byrow=TRUE)
rownames(scores) <- c("logarithmic score","quadratic score","spherical score","ranked probability score","Dawid-Sebastiani score","squared error score","mean absolute error score","root squared error score")
colnames(scores) <- c("Scores")
scores <- as.table(scores)


###create charts

par(mfrow=c(3,3))
plot(lambdahat, type="l",lwd=2, col="red", main="Model fit" )
lines(ydata,lwd=2,col="yellow")
z <-ppois(ydata-1,lambdahat)+(ppois(ydata,lambdahat)-ppois(ydata-1,lambdahat))*runif(length(ydata))
.pit(ydata, ppois(ydata,lambdahat) , ppois(ydata-1,lambdahat))
acf(z,demean=TRUE,main="ACF of PIT")
acf(z^2,demean=TRUE,main="ACF of Squared PIT")
acf(pres,main="ACF of residuals")
acf(pres^2,main="ACF of Squared residuals")
acf(pres^3,main="ACF of Cubic residuals")
scores
}

Berkowitz <- function(ydata, yhatdata, rep, ...)
{
##Perform Berkowitz test
lambdahat<-yhatdata
berkmat <- matrix(NA,rep,1)
for (r in 1:rep) {
z <-ppois(ydata-1,lambdahat)+(ppois(ydata,lambdahat)-ppois(ydata-1,lambdahat))*runif(length(ydata))
u <-qnorm(z)
ar<-arma(u, order = c(1, 0))
arpar <- coef(ar)
res <- residuals(ar)
arssr <- res[2:length(ydata)]%*%res[2:length(ydata)]
s2hat <- (arssr)/ (length(ydata)-2)
s2ut <- s2hat/(1-(arpar[1]^2))
mut <- arpar[2]/(1-arpar[1])
loglR <- -0.5*( u%*%u )
loglUR <- -0.5*log(s2ut)-(((u[1]-mut)^2)/(2*s2ut)) -(length(res)/2)*log(s2hat)-(1/(2*s2hat))*(arssr)
LRberk <- -2*(loglR-loglUR)
berkmat[r] <- LRberk
}
berkowitz <- matrix(c(pchisq(mean(berkmat), df=3, lower.tail=FALSE)),ncol=1,byrow=TRUE)
rownames(berkowitz) <- c("P-Likelihood ratio")
colnames(berkowitz) <- c("Berkowitz test")
berkowitz <- as.table(berkowitz)
berkowitz
}

