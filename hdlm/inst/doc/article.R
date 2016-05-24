### R code from vignette source 'article.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: initial_settings
###################################################
options(stringsAsFactors=FALSE)
options(width=80)
options(prompt="R> ")
library("hdlm")
data("votingRecord")
data("wordDataset")


###################################################
### code chunk number 2: article.Rnw:560-566
###################################################
library('hdlm')
set.seed(1)
x <- matrix(rnorm(100*40),ncol=100)
y <- x[,1] + x[,2] * 0.5 + rnorm(40, sd=0.1)
out <- hdlm(y ~ x, , bootstrap = 2, pval.method="fdr")
out


###################################################
### code chunk number 3: article.Rnw:576-577
###################################################
summary(out)


###################################################
### code chunk number 4: article.Rnw:598-599
###################################################
summary(out, level=2)


###################################################
### code chunk number 5: article.Rnw:609-611
###################################################
output <- capture.output(summary(out, level=3))
cat(output[1:20], sep='\n')


###################################################
### code chunk number 6: fig2
###################################################
par(mfrow=c(2,2))
plot(out)


###################################################
### code chunk number 7: article.Rnw:678-690 (eval = FALSE)
###################################################
## dselector <- function(x,y){
##   n <- nrow(x)
##   p <- ncol(x)
##   lambda <- sqrt(log(ncol(x)+1) * nrow(x)) * sd(as.numeric(y))
## 
##   A <- t(x) %*% x
##   R <- rbind(A, -A)
##   a <- c(as.matrix(t(x) %*% y))
##   r <- c(a-lambda, -a-lambda)
##   beta <- quantreg::rq.fit.fnc(diag(p), rep(0,p), R=R, r=r)$coefficients
##   return(round(beta, 6))
## }


###################################################
### code chunk number 8: article.Rnw:704-705 (eval = FALSE)
###################################################
## hdlm(y ~ x, FITCVFUN = dselector)


###################################################
### code chunk number 9: article.Rnw:743-749 (eval = FALSE)
###################################################
## library("quantreg")
## rq2 <- function(formula) {
##   out <- rq(formula)
##   class(out) <- "rq_new"
##   return(out)
## }


###################################################
### code chunk number 10: article.Rnw:754-759 (eval = FALSE)
###################################################
## summary.rq_new <- function(out) {
##   class(out) <- 'rq'
##   val <- summary.rq(out, se='nid')
##   return(val)
## }


###################################################
### code chunk number 11: article.Rnw:763-764 (eval = FALSE)
###################################################
## out <- hdlm(y ~ x, FUNLM=rq2)


###################################################
### code chunk number 12: article.Rnw:782-785 (eval = FALSE)
###################################################
## LMFUN <- function(x,y) return(glm(y ~ x, family=binomial(link=logit)))
## FUNCVFIT <- function(x,y) return(cv.glmnet(x, y, family='binomial'))
## out <- hdlm(y ~ x, LMFUN = LMFUN, FUNCVFIT = FUNCVFIT)


###################################################
### code chunk number 13: article.Rnw:968-970 (eval = FALSE)
###################################################
## set.seed(1)
## out <- hdlm(Log_Length ~ ., data=wordDataset, bootstrap=5)


###################################################
### code chunk number 14: article.Rnw:973-976
###################################################
set.seed(1)
try_out <- try(out <- hdlm(Log_Length ~ ., data=wordDataset, bootstrap=5), silent=TRUE)
cat(try_out) 


###################################################
### code chunk number 15: article.Rnw:989-993
###################################################
set.seed(1)
out <- hdlm(Log_Length ~ ., data=wordDataset, bootstrap = 5, alpha=1,
    M=20, pval.method = "holm")
summary(out)


###################################################
### code chunk number 16: article.Rnw:1015-1019
###################################################
set.seed(1)
out <- hdglm(Boxer_D_CA ~ ., data=votingRecord, bootstrap=5,
           family='binomial')
summary(out)


###################################################
### code chunk number 17: article.Rnw:1033-1037 (eval = FALSE)
###################################################
## set.seed(1)
## out <- hdglm(Boxer_D_CA ~ ., data=votingRecord, 
##             family='binomial', bayes=TRUE, bayesTune=0.01)
## summary(out)


