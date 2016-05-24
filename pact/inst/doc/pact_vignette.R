### R code from vignette source 'pact_vignette.Rnw'

###################################################
### code chunk number 1: pact_vignette.Rnw:45-48
###################################################
options(width=72)
options(continue=" ")
set.seed(10)


###################################################
### code chunk number 2: pact_vignette.Rnw:59-63
###################################################
library("pact")   ### Load the "pact" R-package

data("GSE10846")  ### Load the dataset
GSE10846[1:5,1:5] ### Display a piece of the data


###################################################
### code chunk number 3: pact_vignette.Rnw:68-72
###################################################
Y <- GSE10846[,1:2]   ## Response, survival status
Treatment <- as.factor(GSE10846[,3]) ## Treatment information
Xv <- GSE10846[,-c(1:3)]   ## Covariates
## No Xf. So Xf=NULL, the default


###################################################
### code chunk number 4: pact_vignette.Rnw:81-88
###################################################
### Fit predictive model using univariate variable selection
p1 <- pact.fit(Y=Y, Xv=Xv, Treatment=Treatment, family="cox",  
      varSelect="univar", nsig=5) 

summary(p1) ## Display model coefficients

print(p1)   ## Print the classification function


###################################################
### code chunk number 5: pact_vignette.Rnw:91-97
###################################################
### Model can be used to predict score for new subjects
r <- rnorm(ncol(Xv))   ## Generate dummy covariate data for one new subject
newXv <- Xv[1,]+r
rownames(newXv) <- "New"

predict(p1, newXv)     ## Now predict scores for this subject


###################################################
### code chunk number 6: pact_vignette.Rnw:102-109
###################################################
### Fit predictive model using "lasso" with peanlty.scaling = 2
p2 <- pact.fit(Y=Y, Xv=Xv, Treatment=Treatment, family="cox",  
      varSelect="lasso", penalty.scaling=2)

summary(p2) ## Display coefficients

print(p2)   ## Print classification function


###################################################
### code chunk number 7: pact_vignette.Rnw:116-118
###################################################
### Cross-validate the 'pact' model, p1
cv1 <- pact.cv(p1, nfold=5)


###################################################
### code chunk number 8: figplot1
###################################################
### Evaluate with method="discrete" (Figure 1)
e1 <- eval.pact.cv(cv1, method="discrete", g=log(0.80), perm.test=FALSE)


###################################################
### code chunk number 9: fig1
###################################################
### Evaluate with method="discrete" (Figure 1)
e1 <- eval.pact.cv(cv1, method="discrete", g=log(0.80), perm.test=FALSE)
  e1


###################################################
### code chunk number 10: figplot2
###################################################
### Evaluation with method="continuous". No cut-offs here.
### Plot type 1: KM curves are plotted at the 20th, 40th, 60th and 80th
### percenttiles of the cross-validated treatment scores (Figure 2)
e21 <- eval.pact.cv(cv1, method="cont", plot.score=TRUE, perm.test=FALSE)


###################################################
### code chunk number 11: fig2
###################################################
### Evaluation with method="continuous". No cut-offs here.
### Plot type 1: KM curves are plotted at the 20th, 40th, 60th and 80th
### percenttiles of the cross-validated treatment scores (Figure 2)
e21 <- eval.pact.cv(cv1, method="cont", plot.score=TRUE, perm.test=FALSE)
  e21


###################################################
### code chunk number 12: figplot3
###################################################
### Evaluate with method="continuous". Plot type 2: Prob[surv] beyond user 
### specified landmark time as a function of the predictive score (Figure 3)
e22 <- eval.pact.cv(cv1, method="cont", plot.score=FALSE, plot.time=12)


###################################################
### code chunk number 13: fig3
###################################################
### Evaluate with method="continuous". Plot type 2: Prob[surv] beyond user 
### specified landmark time as a function of the predictive score (Figure 3)
e22 <- eval.pact.cv(cv1, method="cont", plot.score=FALSE, plot.time=12)
  e22


###################################################
### code chunk number 14: pact_vignette.Rnw:181-183
###################################################
data("prostateCancer")
head(prostateCancer)


###################################################
### code chunk number 15: pact_vignette.Rnw:188-192
###################################################
Y <- prostateCancer[,3:4]  ## Survival response
Xf <- prostateCancer[,7:8]  ## Prognostic covariates always in the model
Xv <- prostateCancer[,c(5:6,9)]  ## Covariates for the predictive score
Treatment <- as.factor(prostateCancer[,2])


###################################################
### code chunk number 16: pact_vignette.Rnw:196-203
###################################################
### Fit predictive model, variable selection with "univar"
p11 <- pact.fit(Y=Y, Xf=Xf, Xv=Xv, Treatment=Treatment, family="cox",
                varSelect="univar")
### And display it
summary(p11)
### Print 
print(p11)


###################################################
### code chunk number 17: pact_vignette.Rnw:206-213
###################################################
### Model can be used to predict score for new subjects
### We only need to specify variables in Xv for new subjects

newXv <- data.frame(age=c(60,70), 
                   pf=c("Normal.Activity","Limited.Activity"),
                   ap=c(0.5,0.5))
predict(p11, newXv)


###################################################
### code chunk number 18: pact_vignette.Rnw:219-224
###################################################
data("EORTC10994")
head(EORTC10994, n=4)
Y <- EORTC10994[,4]   ## Response
Xv <- EORTC10994[,c(2,5,6,7)]  ## Variables in Xv
Treatment <- as.factor(EORTC10994[,3])   ## Treatment


###################################################
### code chunk number 19: pact_vignette.Rnw:229-237
###################################################
### Fit predictive model, no variable selection
pbin <- pact.fit(Y=Y, Xv=Xv, Treatment=Treatment, family="binomial", 
               varSelect="none") 

### Evaluate the model using K-fold CV and method="discrete"
cvbin <- pact.cv(pbin, nfold=5)
e3 <- eval.pact.cv(cvbin, method="discrete", g=log(1), perm.test=FALSE)
e3


###################################################
### code chunk number 20: figplot4
###################################################
### Evaluation for binary response with method="continuous". 
### Plot: Probability of response as a function of cross-validated
### predictive score (Figure 4)
e4 <- eval.pact.cv(cvbin, method="continuous", perm.test=FALSE)


###################################################
### code chunk number 21: fig4
###################################################
### Evaluation for binary response with method="continuous". 
### Plot: Probability of response as a function of cross-validated
### predictive score (Figure 4)
e4 <- eval.pact.cv(cvbin, method="continuous", perm.test=FALSE)
  e4


###################################################
### code chunk number 22: perm
###################################################
### Permutation test examples (survival response): method="discrete"
e5 <- eval.pact.cv(cvbin, method="discrete", g=log(1), 
          perm.test=TRUE, nperm=100)
e5  ### (or print(e5))

### Permutation test examples (survival response): method="continuous"
e6 <- eval.pact.cv(cvbin, method="continuous", perm.test=TRUE, nperm=100)
e6  ### (or print(e6))


