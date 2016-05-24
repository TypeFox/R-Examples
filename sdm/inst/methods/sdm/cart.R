# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March. 2016
# Version 1.0
# Licence GPL v3

#-------------
methodInfo <- list(name=c('cart','CART','tree'),
                   packages='tree',
                   modelTypes = c('pa','pb','ab','n'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(method = "recursive.partition",split = "deviance",x=FALSE,y=FALSE,wts=FALSE,model=FALSE),
                   fitFunction = 'tree',
                   settingRules = function(x,fitSettings) {
                     fitSettings
                   },
                   tuneParams = NULL,
                   predictParams=list(object='model',newdata='sdmDataFrame'),
                   predictSettings=list(type='vector'),
                   predictFunction='predict',
                   #------ metadata (optional):
                   title='Classification or Regression Tree',
                   creator='Babak Naimi',
                   authors=c('B. D. Ripley'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('book',title = "Classification and Regression Trees",
                                          author = as.person("L. Breiman [aut], J.H. Friedman [aut], R.A. Olshen [aut], C.J. Stone [aut]"),
                                          year='1984',
                                          publisher = "Wadsworth"
                   )
                   ),
                   description='A tree is grown by binary recursive partitioning using the response in the specified formula and choosing splits from the terms of the right-hand-side'
)