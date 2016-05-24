# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March. 2016
# Version 1.0
# Licence GPL v3

#-------------
methodInfo <- list(name=c('glmnet','GLMNET','glmelastic','glmlasso'),
                   packages='glmnet',
                   modelTypes = c('pa','pb','ab','n'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(family='binomial',alpha=1),
                   fitFunction = function(formula,data,family,...) {
                     x <- .getData.sdmMatrix(formula,data)
                     y <- .getData.sdmY(formula,data)
                     if (family == 'binomial') m <- cv.glmnet(x,y,family="binomial",type.measure = 'auc')
                     else m <- cv.glmnet(x,y,family=family)
                     glmnet(x=x,y=y,family=family,lambda=m$lambda.1se,...)
                   },
                   settingRules = function(x,fitSettings,predictSettings) {
                     if (x@distribution == 'ab') {
                       fitSettings[['family']] <- 'poisson'
                     } else if (x@distribution == 'n') {
                       fitSettings[['family']] <- 'multinomial'
                     }
                     list(fitSettings=fitSettings,predictSettings=predictSettings)
                   },
                   tuneParams = NULL,
                   predictParams=list(object='model',formula='standard.formula',newx='sdmDataFrame'),
                   predictSettings=list(type='response'),
                   predictFunction=function(object,formula,newx,type) {
                     newx <- .getData.sdmMatrix(formula,newx)
                     predict.glmnet(object,newx,type=type)[,1]
                   },
                   #------ metadata (optional):
                   title='GLM with lasso or elasticnet regularization',
                   creator='Babak Naimi',
                   authors=c('Jerome Friedman, Trevor Hastie, Noah Simon and Rob Tibshirani'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "Regularization Paths for Generalized Linear Models via Coordinate Descent,",
                                          author = as.person("J. Friedman [aut], T. Hastie [aut], R. Tibshirani [aut]"),
                                          year='2008',
                                          journal = "Journal of Statistical Software",
                                          number="33/1",
                                          pages="1-22"
                   )
                   ),
                   description="Fit a generalized linear model via penalized maximum likelihood. The regularization path is computed for the lasso or elasticnet penalty at a grid of values for the regularization parameter lambda."
)