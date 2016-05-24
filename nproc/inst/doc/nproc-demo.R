## ---- eval=FALSE---------------------------------------------------------
#  install.packages("nproc", repos = "http://cran.us.r-project.org")

## ------------------------------------------------------------------------
library(nproc)

## ------------------------------------------------------------------------
n = 1000
set.seed(0)
x = matrix(rnorm(n*2),n,2)
c = 1+3*x[,1]
y = rbinom(n,1,1/(1+exp(-c)))

## ------------------------------------------------------------------------
plot(x[y==1,],col=1,xlim=c(-4,4),xlab='x1',ylab='x2')
points(x[y==0,],col=2,pch=2)
legend("topright",legend=c('Class 1','Class 0'),col=1:2,pch=c(1,2))

## ------------------------------------------------------------------------
fit = npc(x, y, method = "lda", alpha = 0.05)

## ------------------------------------------------------------------------
xtest = matrix(rnorm(n*2),n,2)
ctest = 1+3*xtest[,1]
ytest = rbinom(n,1,1/(1+exp(-ctest)))

## ------------------------------------------------------------------------
pred = predict(fit,xtest)
fit.score = predict(fit,x)
accuracy = mean(pred$pred.label==ytest)
cat("Overall Accuracy: ",  accuracy,'\n')
ind0 = which(ytest==0)
typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
cat('Type I error: ', typeI, '\n')

## ------------------------------------------------------------------------
fit = npc(x, y, method = "logistic", alpha = 0.1)
pred = predict(fit,xtest)
accuracy = mean(pred$pred.label==ytest)
cat("Overall Accuracy: ",  accuracy,'\n')
ind0 = which(ytest==0)
typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
cat('Type I error: ', typeI, '\n')

## ------------------------------------------------------------------------
fit = npc(x, y, method = "logistic", alpha = 0.1, split = 11)
pred = predict(fit,xtest)
accuracy = mean(pred$pred.label==ytest)
cat("Overall Accuracy: ",  accuracy,'\n')
ind0 = which(ytest==0)
typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
cat('Type I error: ', typeI, '\n')

## ------------------------------------------------------------------------
methodlist = c("logistic", "penlog", "svm", "randomforest",
                                        "lda", "nb", "ada")
loc.prob = NULL
for(method in methodlist){
fit = npc(x, y, method = method, alpha = 0.05, loc.prob = loc.prob)
loc.prob = fit$loc.prob #Recycle the loc.prob from the fit to same time for the next fit
pred = predict(fit,xtest)
accuracy = mean(pred$pred.label==ytest)
cat(method, ': Overall Accuracy: ',  accuracy,'\n')
ind0 = which(ytest==0)
typeI = mean(pred$pred.label[ind0]!=ytest[ind0]) #type I error on test set
cat(method, ': Type I error: ', typeI, '\n')
}


## ------------------------------------------------------------------------
fit2 = npc(y = y, score = fit.score$pred.score,
pred.score = pred$pred.score, loc.prob = loc.prob, method = 'custom')

## ------------------------------------------------------------------------
fit = nproc(x, y, method = "svm", loc.prob.lo = loc.prob)
plot(fit)

## ------------------------------------------------------------------------
fit = nproc(x, y, method = "lda", loc.prob.lo = loc.prob)
plot(fit)

## ------------------------------------------------------------------------
fit = nproc(x, y, method = "logistic", conf = TRUE)
plot(fit)

## ------------------------------------------------------------------------
fit = nproc(x, y, method = c('svm','logistic','lda'), conf = T)
plot(fit)

