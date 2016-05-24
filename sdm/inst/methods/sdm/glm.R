# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March. 2016
# Version 1.0
# Licence GPL v3

#-------------


methodInfo <-list(name=c('glm','GLM','lm'),
           packages='stats',
           modelTypes = c('pa','pb','ab','n'),
           fitParams = list(formula='standard.formula',data='sdmDataFrame'),
           fitSettings = list(family=binomial(link='logit'),weights=NULL,model=FALSE),
           fitFunction = 'glm',
           settingRules = function(x,fitSettings) {
             if (x@distribution == 'ab') fitSettings[['family']] <- poisson
             fitSettings
           },
           tuneParams = NULL,
           predictParams=list(object='model',newdata='sdmDataFrame'),
           predictSettings=list(type='response'),
           predictFunction='predict.glm',
           #------ metadata (optional):
           title='Generalized Linear Model',
           creator='Babak Naimi',
           authors=c('R Core team'), # authors of the main method
           email='naimi.b@gmail.com',
           url='http://r-gis.net',
           citation=list(bibentry('book',title = "Generalized linear models",
                                  author = as.person("P. McCullagh [aut], J. A. Nelder [aut]"),
                                  year = "1989",
                                  publisher = "Chapman and Hall",
                                  address = "London")
           ),
           description='glm is used to fit generalized linear models, specified by giving a symbolic description of the linear predictor and a description of the error distribution [see the help for glm function in stats package]'
)