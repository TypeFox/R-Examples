# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March. 2016
# Version 1.0
# Licence GPL v3

#-------------
methodInfo <- list(name=c('rf','RF','randomForest','rforest'),
                   packages='randomForest',
                   modelTypes = c('pa','pb','ab','n'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(ntree=1000,
                                      replace=TRUE,
                                      importance=TRUE
                   ),
                   fitFunction = 'randomForest',
                   settingRules = function(x,fitSettings,predictSettings,userSettings=NULL) {
                     if (!is.null(userSetting)) fitSettings <- .assign(fitSettings,userSettings)
                     
                     list(fitSettings=fitSettings,predictSettings=predictSettings)
                   },
                   tuneParams = NULL,
                   predictParams=list(object='model',newdata='sdmDataFrame'),
                   predictSettings=list(type='response'),
                   predictFunction='predict',
                   #------ metadata (optional):
                   title='Random Forest',
                   creator='Babak Naimi',
                   authors=c('Andy Liaw','Matthew Wiener'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "Random Forests",
                                          author = as.person("L. Breiman [aut]"),
                                          year = "2001",
                                          journal = "Machine Learning",
                                          number="45(1)",
                                          pages="5-32"
                   )
                   ),
                   description="implements Breiman's random forest algorithm (based on Breiman and Cutler's original Fortran code) for classification and regression."
)