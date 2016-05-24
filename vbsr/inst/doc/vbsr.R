### R code from vignette source 'vbsr.rnw'

###################################################
### code chunk number 1: vbsr.rnw:14-26
###################################################
library(vbsr)
set.seed(2)
n <- 100
m <- 95
ntrue <- 10
e <- rnorm(n)
X <- matrix(rnorm(n*m),n,m)
tbeta <- sample(1:m,ntrue)
beta <- rep(0,m)
beta[tbeta]<- rnorm(ntrue,0,2)
y <- X%*%beta+e
res<- vbsr(y,X,family='normal')


###################################################
### code chunk number 2: vbsr.rnw:28-29
###################################################
plot(res$beta,beta)


###################################################
### code chunk number 3: vbsr.rnw:33-35
###################################################
plot(-log10(res$pval),log='y')
lines(c(-10,m+10),c(-log10(0.05/m),-log10(0.05/m)),col='red',lwd=3)


###################################################
### code chunk number 4: vbsr.rnw:39-41
###################################################
cat('True variables:',sort(tbeta),'\n');
cat('Vbsr variables:',which(res$pval<0.05/m),'\n');


###################################################
### code chunk number 5: vbsr.rnw:44-49
###################################################
ols <- lm(y~X);
beta_ols <- summary(ols)$coef[-1,1];
beta_vbsr <- res$beta;
cat('OLS MSE:',mean((beta-beta_ols)^2),'\n');
cat('VBSR MSE:',mean((beta-beta_vbsr)^2),'\n');


###################################################
### code chunk number 6: vbsr.rnw:51-54
###################################################
#barplot(t(cbind(beta[tbeta],summary(ols)$coef[-1,1][tbeta],res$beta[tbeta])),beside=T,col=c('blue','red','green'))
#legend('topleft',c('beta','beta_ols','beta_vbsr'),fill=c('blue','red','green'))
plot((beta-beta_ols)^2,(beta-beta_vbsr)^2)


###################################################
### code chunk number 7: vbsr.rnw:56-57
###################################################
pairs(cbind(beta,beta_ols,beta_vbsr));


###################################################
### code chunk number 8: vbsr.rnw:61-65
###################################################
lmfun <- function(x,y){return(summary(lm(y~x))$coef[2,1]);}
beta_uni <- apply(X,2,lmfun,y);
cat('UNI MSE:',mean((beta-beta_uni)^2),'\n');
cat('VBSR MSE:',mean((beta-beta_vbsr)^2),'\n');


###################################################
### code chunk number 9: vbsr.rnw:71-75
###################################################
g <- rnorm(n);
X <- X+g;
y <- X%*%beta+e
res<- vbsr(y,X,family='normal')


###################################################
### code chunk number 10: vbsr.rnw:77-78
###################################################
plot(res$beta,beta)


###################################################
### code chunk number 11: vbsr.rnw:82-84
###################################################
plot(-log10(res$pval),log='y')
lines(c(-10,m+10),c(-log10(0.05/m),-log10(0.05/m)),col='red',lwd=3)


###################################################
### code chunk number 12: vbsr.rnw:88-90
###################################################
cat('True variables:',sort(tbeta),'\n');
cat('Vbsr variables:',which(res$pval<0.05/m),'\n');


###################################################
### code chunk number 13: vbsr.rnw:93-98
###################################################
ols <- lm(y~X);
beta_ols <- summary(ols)$coef[-1,1];
beta_vbsr <- res$beta;
cat('OLS MSE:',mean((beta-beta_ols)^2),'\n');
cat('VBSR MSE:',mean((beta-beta_vbsr)^2),'\n');


###################################################
### code chunk number 14: vbsr.rnw:100-103
###################################################
#barplot(t(cbind(beta[tbeta],summary(ols)$coef[-1,1][tbeta],res$beta[tbeta])),beside=T,col=c('blue','red','green'))
#legend('topleft',c('beta','beta_ols','beta_vbsr'),fill=c('blue','red','green'))
plot((beta-beta_ols)^2,(beta-beta_vbsr)^2)


###################################################
### code chunk number 15: vbsr.rnw:105-106
###################################################
pairs(cbind(beta,beta_ols,beta_vbsr));


###################################################
### code chunk number 16: vbsr.rnw:110-114
###################################################
lmfun <- function(x,y){return(summary(lm(y~x))$coef[2,1]);}
beta_uni <- apply(X,2,lmfun,y);
cat('UNI MSE:',mean((beta-beta_uni)^2),'\n');
cat('VBSR MSE:',mean((beta-beta_vbsr)^2),'\n');


