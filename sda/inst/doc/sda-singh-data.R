# /*
# This is an R script containing R markdown comments.  It can be run as is in R.
# To generate a document containing the formatted R code, R output and markdown 
# click the "Compile Notebook" button in R Studio, or run the command
# rmarkdown::render() - see http://rmarkdown.rstudio.com/r_notebook_format.html
# */


#' ---
#' title: "Analysis of Singh et al. (2002) Prostate Cancer Data"
#' output: pdf_document
#' author: ""
#' date: Requires "sda" in version 1.3.2 (January 2014) or later
#' ---

#'
#' # Load "sda" package and Singh et al. (2002) data set


library("sda")

#' Singh et al. (2002) gene expression data for prostate cancer:
data(singh2002)

Xtrain = singh2002$x
Ytrain = singh2002$y

dim(Xtrain)      # 102 6033 
length(Ytrain)   # 102
levels(Ytrain)


#' \newpage
#'  
#' # Feature ranking with t-scores

#' Corresponds to assuming a diagonal covariance matrix (DDA).
#'  
#' Compute ranking:
ranking.DDA = sda.ranking(Xtrain, Ytrain, diagonal=TRUE)
ranking.DDA[1:10,]

#' Plot t-scores for the top 40 genes:
#+ fig.height=8
plot(ranking.DDA, top=40) 

#' Number of features with local FDR < 0.8 
#' (i.e. features useful for prediction):
sum(ranking.DDA[,"lfdr"] < 0.8) # 166

#' Number of features with local FDR < 0.2 
#' (i.e. significant non-null features):
sum(ranking.DDA[,"lfdr"] < 0.2) # 53

#' Optimal number of features according to Higher Criticism:
#+ fig.height=6
plot(ranking.DDA[,"HC"], type="l")
which.max( ranking.DDA[1:1000,"HC"] ) #129


#' \newpage
#'  
#' # Feature ranking with correlation-adjusted t-scores (CAT scores)

#' Corresponds to assuming a full covariance matrix (LDA).
#'
#' Compute ranking:
ranking.LDA = sda.ranking(Xtrain, Ytrain, diagonal=FALSE)
ranking.LDA[1:10,]

#' Plot cat scores for the top 40 genes:
#+ fig.height=8
plot(ranking.LDA, top=40) 

#' Number of features with local FDR < 0.8 
#' (i.e. features useful for prediction):
sum(ranking.LDA[,"lfdr"] < 0.8) # 131

#' Number of features with local FDR < 0.2 
#' (i.e. significant non-null features):
sum(ranking.LDA[,"lfdr"] < 0.2) # 62

#' Optimal number of features according to Higher Criticism:
#+ fig.height=6
plot(ranking.LDA[,"HC"], type="l")
which.max( ranking.LDA[1:1000,"HC"] ) # 116




#' \newpage
#'  
#' # Estimate prediction accuracy using crossvalidation 
 

library("crossval")

#' Setup prediction function: estimate the accuracy of a predictor with a fixed number of predictors (note
#' this takes into account the uncertainty in estimating the variable ordering).
predfun = function(Xtrain, Ytrain, Xtest, Ytest, numVars, diagonal=FALSE)
{ 
  # estimate ranking and determine the best numVars variables
  ra = sda.ranking(Xtrain, Ytrain, verbose=FALSE, diagonal=diagonal, fdr=FALSE)
  selVars = ra[,"idx"][1:numVars]

  # fit and predict
  sda.out = sda(Xtrain[, selVars, drop=FALSE], Ytrain, diagonal=diagonal, verbose=FALSE)
  ynew = predict(sda.out, Xtest[, selVars, drop=FALSE], verbose=FALSE)$class 

  # count false and true positives/negatives
  negative = levels(Ytrain)[2] # "healthy"
  cm = confusionMatrix(Ytest, ynew, negative=negative) 
  
  return(cm)  
}

#' Our setup for crossvalidation:
K = 10 # number of folds
B = 20 # number of repetitions

#' Estimate accuracy of LDA using the top 120 features ranked by CAT scores:
set.seed(12345)
cv.lda120 = crossval(predfun, Xtrain, Ytrain, K=K, B=B, numVars=120, diagonal=FALSE, verbose=FALSE)
cv.lda120$stat
# /*
#  FP    TP    TN    FN 
#0.110 4.725 4.890 0.475 
# */
diagnosticErrors(cv.lda120$stat)
# /*
#     acc      sens      spec       ppv       npv       lor 
#0.9426471 0.9086538 0.9780000 0.9772492 0.9114632 6.0917753 
# */

#' Estimate accuracy of DDA using the top 120 features ranked by t scores:
set.seed(12345)
cv.dda120 = crossval(predfun, Xtrain, Ytrain, K=K, B=B, numVars=120, diagonal=TRUE, verbose=FALSE)
cv.dda120$stat
# /*
#    FP    TP    TN    FN 
# 0.170 4.715 4.830 0.485
# */
diagnosticErrors(cv.dda120$stat)
# /*
#      acc      sens      spec       ppv       npv       lor 
#0.9357843 0.9067308 0.9660000 0.9651996 0.9087488 5.6211586
# */

#' Same as before but using only the top 10 features:
set.seed(12345)
cv.dda10 = crossval(predfun, Xtrain, Ytrain, K=K, B=B, numVars=10, diagonal=TRUE, verbose=FALSE)
cv.dda10$stat
# /*
#   FP    TP    TN    FN 
#1.370 3.355 3.630 1.845 
# */
diagnosticErrors(cv.dda10$stat)
# /*
#      acc      sens      spec       ppv       npv       lor 
#0.6848039 0.6451923 0.7260000 0.7100529 0.6630137 1.5723944 
# */


