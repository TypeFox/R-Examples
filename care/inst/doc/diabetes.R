# /*
# This is an R script containing R markdown comments.  It can be run as is in R.
# To generate a document containing the formatted R code, R output and markdown 
# click the "Compile Notebook" button in R Studio, or run the command
# rmarkdown::render() - see http://rmarkdown.rstudio.com/r_notebook_format.html
# */


#' ---
#' title: "Diabetes Data"
#' output: pdf_document
#' author: ""
#' date: Requires "care" version 1.1.7 (November 2014) or later
#' ---


#' This R script reproduces the analysis of diabetes data from 
#' V. Zuber and K. Strimmer. 2011. *High-dimensional regression and variable 
#' selection using CAR scores.* Statist. Appl. Genet. Mol. Biol. **10**: 34
#' (http://dx.doi.org/10.2202/1544-6115.1730)


#'
#' # Load "care" package and diabetes data set

library("care")

#' Diabetes data (442 patients) from Efron et al. 2004. *Least angle regression*
#' Ann. Statist. **32**:407-499.
data(efron2004)
x = efron2004$x
dim(x)
d = ncol(x) # dimension 
n = nrow(x) # samples

#' 10 predictors
xnames = colnames(x)
xnames
# /* "age" "sex" "bmi" "bp"  "s1"  "s2"  "s3"  "s4"  "s5"  "s6" */

#' Response
y = efron2004$y
length(y)

#' \newpage
#'
#' # Comparison of linear regression models

#' Ordering of predictors according to CAR score:
car = carscore(x, y, lambda=0) # no shrinkage estimation needed as n>>d
ocar = order(car^2, decreasing=TRUE)
xnames[ocar]

#' Regression coefficients for models with increasing number of predictors:
car.predlist = make.predlist(ocar, numpred = 1:d, name="CAR")
cm = slm(x, y, car.predlist, lambda=0, lambda.var=0, verbose=FALSE)
bmat= cm$coefficients[,-1]
bmat

#' Plot regression coefficients:
#+ fig.width=7, fig.height=7
plot(1:d, bmat[,1], type="l", 
  ylab="estimated regression coefficients", 
  xlab="number of included predictors", 
  main="CAR Regression Models for Diabetes Data", 
  xlim=c(1,d+1), ylim=c(min(bmat), max(bmat)))

for (i in 2:d) lines(1:d, bmat[,i], col=i, lty=i)
for (i in 1:d) points(1:d, bmat[,i], col=i)
for (i in 1:d) text(d+0.5, bmat[d,i], xnames[i])


#' \newpage
#'
#' # Estimate prediction errors by crossvalidation

library("crossval")

K=10 # number of folds
B=50 # number of repetitions


#' Prediction function used in crossvalidation: Rank by CAR scores, 
#' then fit and predict using a specified number of predictors
#' (note this takes into account the uncertainty in selection and ordering
#' of the predictors)
predfun = function(Xtrain, Ytrain, Xtest, Ytest, numVars)
{  
  # rank the variables according to squared CAR scores
  car = carscore(Xtrain, Ytrain, verbose=FALSE, lambda=0)
  ocar = order(car^2, decreasing=TRUE)
  selVars = ocar[1:numVars]

  # fit and predict
  slm.fit = slm(Xtrain, Ytrain, predlist=list(selVars), verbose=FALSE, 
       lambda=0, lambda.var=0)
  Ynew = predict(slm.fit, Xtest, verbose=FALSE)

  # compute squared error risk
  mse = mean( (Ynew - Ytest)^2)  

  return(mse)  
}

#' Perform crossvalidation:
numpred = 1:10 # number of predictors
set.seed(12345)
cvsim = lapply(numpred, 
  function(i)
  {
    cat("Number of predictors:", i, "\n")
    cvp = crossval(predfun, x, y, K=K, B=B, numVars = i, verbose=FALSE)
    return( cvp$stat.cv )
  }
)

#' Plot results:
boxplot(cvsim, names=numpred,
col=c(rep("grey", 2), rep("white", 8)),
 main="CAR Models for the Diabetes Data", xlab="number of included predictors",
       ylab="estimated CV prediction error")

#' **Conclusion:** after including the three top ranked predictors ("bmi", "s5", "bp") 
#' no further reduction of MSE is seen.




