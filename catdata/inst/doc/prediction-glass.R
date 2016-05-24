### R code from vignette source 'prediction-glass.Rnw'

###################################################
### code chunk number 1: prediction-glass.Rnw:13-15 (eval = FALSE)
###################################################
## options(width=60)
## set.seed(6)


###################################################
### code chunk number 2: prediction-glass.Rnw:18-21 (eval = FALSE)
###################################################
## library(catdata)
## data(glass)
## source("disc_comp_roc.r")


###################################################
### code chunk number 3: prediction-glass.Rnw:25-26 (eval = FALSE)
###################################################
## table(glass$type)


###################################################
### code chunk number 4: prediction-glass.Rnw:29-36 (eval = FALSE)
###################################################
## # Computation with only three response categories:
## data(glass)
## 
## 
## glass <- glass[glass$type == "type1" | glass$type == "type2" | glass$type == "type7", ]
## glass[,10] <- factor(glass[,10], labels = c("type1", "type2", "type7"))
## x <- scale(as.matrix(glass[,1:9]), center = T, scale = T)


###################################################
### code chunk number 5: prediction-glass.Rnw:38-39 (eval = FALSE)
###################################################
## test <- disc.comp(glass$type, x, methods = c("lda", "qda", "rf", "neural", "sv", "knn"), nfold = 5, ntimes = 50)


###################################################
### code chunk number 6: prediction-glass.Rnw:44-47 (eval = FALSE)
###################################################
## par(cex.lab = 1.5, cex.axis = 1.5, mai = c(0.9,0.9,0.5,0.5))
## #misclass error
## boxplot(test$lda$misclass, test$qda$misclass, test$rf$misclass, test$neural$misclass, test$sv$misclass, test$knn$misclass, names = c("lda", "qda", "rf", "nnet", "sv", "knn"), main = "")


###################################################
### code chunk number 7: prediction-glass.Rnw:50-52 (eval = FALSE)
###################################################
## #squared error
## boxplot(test$lda$sqerr, test$qda$sqerr, test$neural$sqerr, test$sv$sqerr, names = c("lda", "qda", "nnet", "sv"), main = "")


###################################################
### code chunk number 8: prediction-glass.Rnw:55-62 (eval = FALSE)
###################################################
## # Computation with binary response:
## data(glass)
## 
## 
## glass <- glass[glass$type == "type1" | glass$type == "type2",]
## glass[,10] <- factor(glass[,10],labels=c("type1", "type2"))
## x <- scale(as.matrix(glass[,1:9]), center = T, scale = T)


###################################################
### code chunk number 9: prediction-glass.Rnw:64-65 (eval = FALSE)
###################################################
## test <- disc.comp(glass$type, x, nfold = 5, ntimes = 50)


###################################################
### code chunk number 10: prediction-glass.Rnw:68-70 (eval = FALSE)
###################################################
## par(cex.lab = 1.5, cex.axis = 1.5, mai = c(0.9,0.9,0.5,0.5))
## roc.curve(test, roc.k = c("lda", "boost", "sv"))


###################################################
### code chunk number 11: prediction-glass.Rnw:73-75 (eval = FALSE)
###################################################
##  par(cex.main = 1.5, cex.lab = 1.4, cex.axis = 1.4)
## boxplot(test$lda$misclass, test$qda$misclass, test$logistic$misclass, test$rf$misclass, test$boost$misclass, test$neural$misclass, test$sv$misclass, test$knn$misclass, names = c("lda", "qda", "log", "rf", "bst", "nnet", "sv", "knn"), main = "")


###################################################
### code chunk number 12: prediction-glass.Rnw:78-79 (eval = FALSE)
###################################################
## boxplot(test$lda$sqerr, test$qda$sqerr, test$logistic$sqerr,  test$boost$sqerr, test$neural$sqerr, test$sv$sqerr,  names = c("lda", "qda", "log",  "boost", "nnet", "sv" ), main = "")


