# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March. 2016
# Version 1.0
# Licence GPL v3

#-------------
methodInfo <- list(name=c('fda','FDA'),
                   packages='mda',
                   modelTypes = c('pa','pb','n'),
                   fitParams = list(formula='standard.formula',data='sdmDataFrame'),
                   fitSettings = list(method = substitute(polyreg),keep.fitted=FALSE),
                   fitFunction = 'fda',
                   settingRules = function(x,fitSettings) {
                     fitSettings
                   },
                   tuneParams = list(method=substitute(c(polyreg,mars,gen.ridge,bruto))),
                   predictParams=list(object='model',newdata='sdmDataFrame'),
                   predictSettings=list(type='posterior'),
                   predictFunction=function(object,newdata,type) {
                     predict(object,newdata,type=type)[,'1']
                   },
                   #------ metadata (optional):
                   title='Flexible Discriminant Analysis',
                   creator='Babak Naimi',
                   authors=c('Trevor Hastie; Robert Tibshirani'), # authors of the main method
                   email='naimi.b@gmail.com',
                   url='http://r-gis.net',
                   citation=list(bibentry('Article',title = "Flexible Disriminant Analysis by Optimal Scoring",
                                          author = as.person("T. Hastie [aut], R. Tibshirani [aut], Buja [aut]"),
                                          year='1994',
                                          journal = "JASA",
                                          pages="1255-1270"
                   )
                   ),
                   description='Flexible discriminant analysis.'
)