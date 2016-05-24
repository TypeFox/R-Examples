# /*
# This is an R script containing R markdown comments.  It can be run as is in R.
# To generate a document containing the formatted R code, R output and markdown 
# click the "Compile Notebook" button in R Studio, or run the command
# rmarkdown::render() - see http://rmarkdown.rstudio.com/r_notebook_format.html
# */


#' ---
#' title: "Analysis of SRBCT Data"
#' output: pdf_document
#' author: ""
#' date: Requires "sda" in version 1.3.2 (January 2014) or later
#' ---

#'
#' # Load "sda" package and create SRBCT data set

library("sda")

#' Load data set from Khan et al. (2001):
data(khan2001)

#' Create data set containing only the SRBCT samples:
del.idx = which( khan2001$y == "non-SRBCT" )
srbct.x = khan2001$x[-del.idx,]
srbct.y = factor(khan2001$y[-del.idx])
dim(srbct.x)
 
#' Four subtypes of cancer:
levels(srbct.y)

#' Divide into training and test data
Xtrain = srbct.x[1:63,]
Ytrain = srbct.y[1:63]
Xtest = srbct.x[64:83,]
Ytest = srbct.y[64:83]


#' \newpage
#'
#' # Diagonal Discriminant Analysis (DDA) 

#' In DDA correlation among predictors is assumed to be zero, i.e. a diagonal 
#' covariance matrix is used.
 
 
#'
#' ## Step 1 - feature ranking

#' As there are more than two groups in the response there are three different
#' ways to obtain a summary test statistic to rank genes:

#' a) ranking by averaged squared t-scores across the four groups
#+ fig.height=5
ra = sda.ranking(Xtrain, Ytrain, fdr=TRUE, plot.fdr=TRUE, diagonal=TRUE, ranking.score="avg")
sum( ra[, "lfdr"]< 0.80)  # 97 genes included in classifier (by FNDR control)
which.max( ra[, "HC"] )  #  145 genes according to HC criterion

#' b) ranking by maximum of squared t-scores across the four groups
#+ fig.height=5
ra = sda.ranking(Xtrain, Ytrain, fdr=TRUE, plot.fdr=TRUE, diagonal=TRUE, ranking.score="max")
sum( ra[, "lfdr"]< 0.80)  # 78 genes included in classifier (by FNDR control)
which.max( ra[, "HC"] )  #  121 genes according to HC criterion

#' c) ranking by mutual information (weighted sum of squared t-scores)
#+ fig.height=5
ra = sda.ranking(Xtrain, Ytrain, fdr=TRUE, plot.fdr=TRUE, diagonal=TRUE, ranking.score="entropy")
sum( ra[, "lfdr"]< 0.80)  # 99 genes included in classifier (by FNDR control)
which.max( ra[, "HC"] )  #  158 genes according to HC criterion


#' here we pick the top 99 genes of option c)
#+ fig.height=8
plot(ra, top=99, main="The 99 Top Ranking Genes", ylab="Gene ID")

#' Select these 99 variables:
idx = ra[1:99,"idx"]
Xtrain2 = Xtrain[,idx]
Xtest2 = Xtest[,idx]

 
#'
#' ## Step 2 - training the classifier

#' Learn DDA predictor:
sda.fit = sda(Xtrain2, Ytrain, diagonal=TRUE)

#'
#' ## Step 3 - prediction

#' Predict class labels from test data and compare with known labels:
dim(Xtest2)
predict(sda.fit, Xtest2) 
ynew = predict(sda.fit, Xtest2)$class

#' Number of missclassified test samples:
sum(ynew != Ytest)
# /* 1 */


#' \newpage
#'
#' # Linear Discriminant Analysis (LDA) 

#' In LDA correlation among predictors is taken into account.
 
#'
#' ## Step 1 - feature ranking

#' As there are more than two groups in the response there are three different
#' ways to obtain a summary test statistic to rank genes:
 
#' a) ranking by averaged squared cat-scores across the four groups
#+ fig.height=5
ra = sda.ranking(Xtrain, Ytrain, fdr=TRUE, plot.fdr=TRUE, ranking.score="avg")
sum( ra[, "lfdr"]< 0.80)  # 93 genes included in classifier (by FNDR control)
which.max( ra[, "HC"] )  #  143 genes according to HC criterion
 
#' b) ranking by maximum of squared cat-scores across the four groups
#+ fig.height=5
ra = sda.ranking(Xtrain, Ytrain, fdr=TRUE, plot.fdr=TRUE, ranking.score="max")
sum( ra[, "lfdr"]< 0.80)  # 156 genes included in classifier (by FNDR control)
which.max( ra[, "HC"] )  #  194 genes according to HC criterion
 
#' c) ranking by mutual information (weighted sum of squared cat-scores)
#+ fig.height=5
ra = sda.ranking(Xtrain, Ytrain, fdr=TRUE, plot.fdr=TRUE, ranking.score="entropy")
sum( ra[, "lfdr"]< 0.80)  # 97 genes included in classifier (by FNDR control)
which.max( ra[, "HC"] )  #  140 genes according to HC criterion
 
#' here we pick the top 97 genes of option c)
#+ fig.height=8
plot(ra, top=97, main="The 97 Top Ranking Genes", ylab="Gene ID")

#' Select these 97 variables:
idx = ra[1:97,"idx"]
Xtrain2 = Xtrain[,idx]
Xtest2 = Xtest[,idx]
 
#'
#' ## Step 2 - training the classifier

#' Learn LDA predictor:
sda.fit = sda(Xtrain2, Ytrain)

#'
#' ## Step 3 - prediction
 
#' Predict class labels from test data and compare with known labels:
dim(Xtest2)
predict(sda.fit, Xtest2)
ynew = predict(sda.fit, Xtest2)$class

#' Number of missclassified test samples:
sum(ynew != Ytest) 
# /* 0 */

 
#' \newpage
#'
#' # Estimate prediction accuracy using crossvalidation 
 
#' Using crossvalidation we can estimate the prediction error
#' from the training data set alone.
 
library("crossval")
  
 
#' Setup prediction function: estimate the accuracy of a predictor with a fixed number of predictors (note
#' this takes into account the uncertainty in estimating the variable ordering).
predfun = function(Xtrain, Ytrain, Xtest, Ytest, numVars, diagonal=FALSE,
                    ranking.score="entropy")
{ 
  # estimate ranking and determine the best numVars variables
  ra = sda.ranking(Xtrain, Ytrain, verbose=FALSE, diagonal=diagonal, 
                   fdr=FALSE, ranking.score=ranking.score)
  selVars = ra[,"idx"][1:numVars]
   
  # fit and predict
  sda.out = sda(Xtrain[, selVars, drop=FALSE], Ytrain, diagonal=diagonal,
                verbose=FALSE)
  ynew = predict(sda.out, Xtest[, selVars, drop=FALSE], verbose=FALSE)$class 
   
  # compute accuracy
  acc = mean(Ytest == ynew)
   
  return(acc)  
}
 
#' Our setup for crossvalidation:
K = 10 # number of folds
B = 20 # number of repetitions
 
 
#'  Crossvalidation estimate of accuracy for
#'  LDA using the top 100 features ranked by CAT scores
#'  (combined across groups using "entropy" for overall ranking):
set.seed(12345)
cv.lda100 = crossval(predfun, Xtrain, Ytrain, K=K, B=B, numVars=100, 
                     diagonal=FALSE, verbose=FALSE)
cv.lda100$stat
# /* 1 */
 
 
#'
#' ## Comparison of LDA / DDA and "entropy" and "max" options
 
#' LDA using the top 10 features ranked by CAT scores
#' (combined across groups using "entropy" for overall ranking):
set.seed(12345)
cv.lda10 = crossval(predfun, Xtrain, Ytrain, K=K, B=B, numVars=10, 
                     diagonal=FALSE, verbose=FALSE)
cv.lda10$stat
# /* 0.9909762 */
 
#' DDA using the top 10 features ranked by t scores
#' (combined across groups using "entropy" for overall ranking):
set.seed(12345)
cv.dda10 = crossval(predfun, Xtrain, Ytrain, K=K, B=B, numVars=10, 
                     diagonal=TRUE, verbose=FALSE)
cv.dda10$stat
# /* 0.9643869 */
 
#' DDA using the top 10 features ranked by t scores,
#' (combined across groups using "max" for overall ranking, as in PAM):
set.seed(12345)
cv.dda10b = crossval(predfun, Xtrain, Ytrain, K=K, B=B, numVars=10, 
                      diagonal=TRUE, ranking.score="max", verbose=FALSE)
cv.dda10b$stat
# /* 0.9585595 */
 
#' **Conclusions**:
#' 
#' 1. LDA/CAT score ranking performs petter than DDA/t-score ranking.
#' 2. "entropy" is better as group summary than "max".
