# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March. 2016
# Version 1.0
# Licence GPL v3

#-------------
methodInfo <- list(name=c('maxlike','MaxLike'),
                   packages=NULL,
                   modelTypes = c('pb'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(link='logit',hessian=TRUE,normalize=TRUE),
                   fitFunction = '.maxlike',
                   settingRules = function(x,fitSettings,predictSettings) {
                     #
                   },
                   tuneParams = NULL,
                   predictParams=list(object='model',newdata='sdmDataFrame'),
                   predictSettings=NULL,
                   predictFunction='predict',
                   #------ metadata (optional):
                   title='Model occurrence probability using presence-only data',
                   creator='Babak Naimi',
                   authors=c('Richard Chandler and Andy Royle'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "Likelihood analysis of species occurrence probability from presence-only data for modelling species distributions",
                                          author = as.person("J.A. Royle [aut], R.B. Chandler [aut], C. Yackulic [aut], J. D. Nichols [aut]"),
                                          year='2012',
                                          journal = "Methods in Ecology and Evolution"
                                          
                   )
                   ),
                   description="Estimates the probability of occurrence using presence-only data and covariates as background data."
)