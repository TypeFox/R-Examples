# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March. 2016
# Version 1.0
# Licence GPL v3

#-------------
methodInfo <- list(name=c('mars','MARS','earth'),
                   packages='earth',
                   modelTypes = c('pa','pb','ab','n'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(weights=NULL,pmethod='none',trace=0,glm=list(family=binomial)),
                   fitFunction = 'earth',
                   settingRules = function(x,fitSettings) {
                     if (x@distribution == 'ab') fitSettings[['glm']] <- list(family=poisson)
                     fitSettings
                   },
                   tuneParams = list(glm=c(NULL,'default')),
                   predictParams=list(object='model',newdata='sdmDataFrame'),
                   predictSettings=list(type='response'),
                   predictFunction='predict',
                   #------ metadata (optional):
                   title='Multivariate Adaptive Regression Splines',
                   creator='Babak Naimi',
                   authors=c('Stephen Milborrow'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "Multivariate Adaptive Regression Splines",
                                          author = as.person("Friedman [aut]"),
                                          year='1991',
                                          journal = "Annals of Statistics",
                                          number="19/1",
                                          pages="1-141"
                   )
                   ),
                   description="Build a regression model using the techniques in Friedman's papers 'Multivariate Adaptive Regression Splines' and 'Fast MARS'."
)