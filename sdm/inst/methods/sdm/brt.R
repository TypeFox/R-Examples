# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March. 2016
# Version 1.0
# Licence GPL v3

#-------------
methodInfo <- list(name=c('brt','BRT','gbm','GBM'),
                  packages='gbm',
                  modelTypes = c('pa','pb','ab','n'),
                  fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                  fitSettings = list(distribution='bernoulli',
                                     n.trees=1000,
                                     interaction.depth=1,
                                     n.minobsinnode = 10,
                                     shrinkage = 0.001,
                                     bag.fraction = 0.5,
                                     train.fraction = 1.0,
                                     cv.folds=0,
                                     keep.data = TRUE,
                                     verbose = "CV",
                                     class.stratify.cv=NULL),
                  fitFunction = 'gbm',
                  settingRules = function(x,fitSettings,predictSettings,userSettings=NULL) {
                    if (!is.null(userSetting)) fitSettings <- .assign(fitSettings,userSettings)
                    else if (x@distribution == 'ab') fitSettings[['distribution']] <- "poisson"
                    else if (x@distribution == 'n') fitSettings[['distribution']] <- "multinomial"
                    
                    list(fitSettings=fitSettings,predictSettings=predictSettings)
                  },
                  tuneParams = NULL,
                  predictParams=list(object='model',newdata='sdmDataFrame'),
                  predictSettings=list(n.trees=1000,type='response'),
                  predictFunction='predict.gbm',
                  #------ metadata (optional):
                  title='Boosted Regression Trees',
                  creator='Babak Naimi',
                  authors=c('Greg Ridgeway'), # authors of the main method
                  email='naimi.b@gmail.com',
                  url='http://r-gis.net',
                  citation=list(bibentry('Article',title = "A working guide to boosted regression trees",
                                         author = as.person("J. Elith [aut], J. R. Leathwick [aut], T. Hastie [aut]"),
                                         year = "2008",
                                         journal = "Journal of Animal Ecology",
                                         number="77",
                                         pages="802-813",
                                         publisher="Wiley Online Library")
                  ),
                  description='Fits Boosting regression trees (BRT), called also generalized boosting regression model (GBM). Boosting is the process of iteratively adding basis functions in a greedy fashion so that each additional basis function further reduces the selected loss function [see the help for gbm function in gbm package]'
)
#------------