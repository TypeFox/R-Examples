# Author: Babak Naimi, naimi.b@gmail.com
# Date :  March. 2016
# Version 1.0
# Licence GPL v3

#-------------
methodInfo <-list(name=c('gam','GAM'),
           packages='mgcv',
           modelTypes = c('pa','pb','ab','n'),
           fitParams = list(formula='gam.mgcv.formula',data='sdmDataFrame'),
           fitSettings = list(family=binomial(link='logit'),weights=NULL,subset=NULL,na.action='na.omit',offset=NULL,method='GCV.Cp',optimizer=c("outer","newton"),select=FALSE,knots=NULL,sp=NULL,min.sp=NULL,H=NULL,gamma=1,fit=TRUE,paraPen=NULL,G=NULL),
           fitFunction = 'gam',
           settingRules = function(x='sdmSettings',y='sdmChategories') {
             o <- list()
             if (x@distribution == 'ab') o[['family']] <- poisson
             o
           },
           tuneParams = NULL,
           predictParams=list(object='model',newdata='sdmDataFrame'),
           predictSettings=list(type='response'),
           predictFunction='predict.gam',
           #------ metadata (optional):
           title='Generalized Additive Models with integrated smoothness estimation',
           creator='Babak Naimi',
           authors=c('Simon N. Wood'), # authors of the main method
           email='naimi.b@gmail.com',
           url='http://r-gis.net',
           citation=list(bibentry('book',title = " Generalized Additive Models: An Introduction with R",
                                  author = as.person("S. N. Wood [aut]"),
                                  year = "2006",
                                  publisher = "Chapman and Hall/CRC press")
           ),
           description='Fits a generalized additive model (GAM) to data, the term "GAM" being taken to include any quadratically penalized GLM and a variety of other models estimated by a quadratically penalised likelihood type approach [see the help for gam function in mgcv package]'
)
#------------