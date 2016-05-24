# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Oct. 2012
# Version 1.1
# Licence GPL v3

setMethod ('show' , 'VIF',
           function ( object ) {
             if (length(object@excluded) > 0) {
               cat (length(object@excluded),'variables from the',length(object@variables), 'input variables have collinearity problem:','\n','\n')
               cat (object@excluded,'\n')
             } else cat ('No variable from the',length(object@variables), 'input variables has collinearity problem.','\n')
             cat('\n')
             if (length(object@excluded) > 0) cat('After excluding the collinear variables, the linear correlation coefficients ranges between:','\n')
             else cat('The linear correlation coefficients ranges between:','\n')
             mx <- .minCor(object@corMatrix)
             cat ('min correlation (',mx[1],'~',mx[2],'): ',object@corMatrix[mx[1],mx[2]], '\n')
             mx <- .maxCor(object@corMatrix)
             cat ('max correlation (',mx[1],'~',mx[2],'): ',object@corMatrix[mx[1],mx[2]], '\n')
             cat ('\n')
             cat('---------- VIFs of the remained variables --------','\n')
             print(object@results)
           }
           )

setMethod ('show' , 'speciesLISA', 
           function(object) {
             cat('class                             :' , class(object), '\n')
             cat('LISA statistic                    :' , object@statistic, '\n\n')
             cat('number of species observations    : ' , nrow(object@LISAs), '\n')
             cat('number of predictor variables     : ' , ncol(object@LISAs), '\n')
             cat('min, mean, max of aggregated LISA : ' , min(object@LISA),',' ,mean(object@LISA),',',max(object@LISA), '\n')
           }
)

setMethod ('show' , 'RasterVariogram', 
           function(object) {
             cat('class               :' , class(object), '\n')
             cat('Lag size            :' , object@lag, '\n')
             cat('Number of lags      : ' , object@nlags, '\n')
             cat ('\n')
             cat('------ Variogram data ------','\n')
             print(object@variogram)
           }
)