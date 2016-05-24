###################################################
### chunk number 1: 
###################################################
library(tgp)
##options(width=65)
seed <- 0; set.seed(seed)


###################################################
### chunk number 2: 
###################################################
# 1-d linear data input and predictive data
X <- seq(0,1,length=50)  # inputs
XX <- seq(0,1,length=99) # predictive locations
Z <- 1 + 2*X + rnorm(length(X),sd=0.25) # responses


###################################################
### chunk number 3: 
###################################################
lin.blm <- blm(X=X, XX=XX, Z=Z)


###################################################
### chunk number 4: blm
###################################################
plot(lin.blm, main='Linear Model,', layout='surf')
abline(1,2,lty=3,col='blue')


###################################################
### chunk number 5: 
###################################################
rl <- readline("press RETURN to continue: ")
graphics.off()


###################################################
### chunk number 6: 
###################################################
lin.gpllm <- bgpllm(X=X, XX=XX, Z=Z)


###################################################
### chunk number 7: gplm
###################################################
plot(lin.gpllm, main='GP LLM,', layout='surf')
abline(1,2,lty=4,col='blue')


###################################################
### chunk number 8: 
###################################################
rl <- readline("press RETURN to continue: ")
graphics.off()


###################################################
### chunk number 9: 
###################################################
lin.gpllm.tr <- bgpllm(X=X, XX=0.5, Z=Z, pred.n=FALSE, trace=TRUE,
                       verb=0)
mla <- mean(lin.gpllm.tr$trace$linarea$la)
mla


###################################################
### chunk number 10: 
###################################################
1-mean(lin.gpllm.tr$trace$XX[[1]]$b1)


