# /*
# This is an R script containing R markdown comments.  It can be run as is in R.
# To generate a document containing the formatted R code, R output and markdown 
# click the "Compile Notebook" button in R Studio, or run the command
# rmarkdown::render() - see http://rmarkdown.rstudio.com/r_notebook_format.html
# */


#' ---
#' title: "Brain Data"
#' output: pdf_document
#' author: ""
#' date: Requires "care" version 1.1.7 (November 2014) or later
#' ---


#' This R script reproduces the analysis of brain gene expression data from 
#' V. Zuber and K. Strimmer. 2011. *High-dimensional regression and variable 
#' selection using CAR scores.* Statist. Appl. Genet. Mol. Biol. **10**: 34
#' (http://dx.doi.org/10.2202/1544-6115.1730)

#'
#' # Load "care" package and brain gene expression data set

library("care")

#' Load Lu et al. (2004) data set:
data(lu2004)
x = lu2004$x
y = lu2004$y # age
dim(x) 
# /* 30 403 */

#' Regularization parameters used to fit linear models:
reg = slm(x,y)$regularization
# /* 0.1373 and 0.0246 */
lambda = reg[1,"lambda"] # correlation shrinkage 
lambda.var = reg[1, "lambda.var"] # variance shrinkage 

#' Compute CAR scores to rank predictors: 
car = carscore(x, y, lambda=lambda)
ocar = order(car^2, decreasing=TRUE)
ocar[1:30]

#' Fit linear models with increasing number of predictors:
numpred =  c(seq(5, 45, by=10), seq(53, 403, by=25))
numpred
car.predlist = make.predlist(ocar, numpred, name="CAR")
car.models = slm(x, y, car.predlist, lambda=lambda, lambda.var=lambda.var)

#' For comparison, use marginal correlations to rank predictors:
marg = carscore(x, y, lambda=lambda, diagonal=TRUE) # shrinking not actually reqired
omarg = order(marg^2, decreasing=TRUE)
marg.predlist = make.predlist(omarg, numpred, name="MARG")
marg.models = slm(x, y, marg.predlist, lambda=lambda, lambda.var=lambda.var)

#' Plot model fitted $R^2$ for all models comparing marginal correlation with CAR scores:
#+ fig.width=8, fig.height=8
ylim = range( c(marg.models$R2, car.models$R2) )
plot(car.models$numpred, car.models$R2, type="b", ylim=ylim,
  xlab="number of included predictors",
  ylab="model fit (R2)",
  main="CAR and Marginal Correlation Models for Brain Data")
points(marg.models$numpred, marg.models$R2, col=2, type="b")
legend("bottomright", c("CAR Score", "Marginal Correlation"), col=c(1,2), lty=c(1,1) )


#' \newpage
#'
#' # Estimate prediction error by crossvalidation

library("crossval")

#' Standardize data following Zuber and Strimmer (2011):
xs = scale(x)
ys = scale(y)

K=5  # number of folds
B=100 # number of repetitions


#' Rank by CAR scores, fit and predict using a specified number of predictors:
predfun = function(Xtrain, Ytrain, Xtest, Ytest, numVars)
{  
  # rank the variables according to squared CAR scores
  car = carscore(Xtrain, Ytrain, verbose=FALSE)
  ocar = order(car^2, decreasing=TRUE)
  selVars = ocar[1:numVars]

  # fit and predict
  slm.fit = slm(Xtrain, Ytrain, predlist=list(selVars), verbose=FALSE)
  Ynew = predict(slm.fit, Xtest, verbose=FALSE)

  # compute squared error risk
  mse = mean( (Ynew - Ytest)^2)  

  return(mse)  
}

#' Compute results from Table 9 in Zuber and Strimmer (2011):
set.seed(12345)
cvp = crossval(predfun, xs, ys, K=K, B=B, numVars = 36, verbose=FALSE)
c(cvp$stat, cvp$stat.se) # 0.3441  0.0074
set.seed(12345)
cvp = crossval(predfun, xs, ys, K=K, B=B, numVars = 60, verbose=FALSE)
c(cvp$stat, cvp$stat.se) # 0.3129  0.0066
set.seed(12345)
cvp = crossval(predfun, xs, ys, K=K, B=B, numVars = 85, verbose=FALSE)
c(cvp$stat, cvp$stat.se) # 0.2978  0.0062 
  

#' This produces Figure 3 in Zuber and Strimmer (2011)
numpred = c(seq(10, 200, 10), 403) # number of predictors
set.seed(12345)
cvsim = lapply(numpred, 
  function(i)
  {
    cat("Number of predictors:", i, "\n")
    cvp = crossval(predfun, xs, ys, K=K, B=B, numVars = i, verbose=FALSE)
    return( cvp$stat.cv )
  }
)
#+ fig.width=8, fig.height=6
boxplot(cvsim, names=numpred,
col=c(rep("grey", 4), rep("white", 16), "grey"),
 main="CAR Models for the Gene Expression Data", xlab="number of included predictors",
       ylab="estimated CV prediction error")



