#' Multinomial Mixed Effects Models
#'
#' The mme package implements three multinomial area level mixed effects models for 
#' small area estimation. The first model (Model 1) is based on the area level multinomial 
#' mixed model with independent random effects for each category of the response 
#' variable (Lopez-Vizcaino et al, 2013). The second model (Model 2) takes advantage 
#' from the availability of survey data from different time periods and uses 
#' a multinomial mixed model with independent random effects for each category of the
#' response variable and with independent time and domain random effects. The 
#' third model (Model 3) is similar to the second one, but with correlated time random effects. To fit the models, we combine the penalized
#' quasi-likelihood (PQL) method, introduced by Breslow and Clayton (1993) for estimating and predicting th fixed
#' and random effects, with the residual maximum likelihood (REML) method for estimating the variance components. 
#' In all models the package use two approaches to estimate the mean square
#' error (MSE), first through an analytical expression and second by bootstrap techniques. 
#'
#' @references Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013). 
#' Multinomial-based small area estimation of labour force indicators. 
#' Statistical Modelling, 13 ,153-178.
#' @references
#' Lopez-Vizcaino, ME, Lombardia, MJ and Morales, D (2013).
#' Small area estimation of labour force indicator under a multinomial mixed model
#' with correlated time and area effects. Submitted for review.
#' @references
#' Breslow, N, Clayton, D (1993).
#' Aproximate inference in generalized linear mixed models. 
#' Journal of the American Statistical Association, 88, 9-25.
#' @name mme-package
#' @aliases mme
#' @docType package
#' @import Matrix
#' @import MASS
NULL

# Data set information
#' Dataset for Model 1
#'
#' Dataset used by the multinomial mixed effects model with one independent random effect in 
#' each category of the response variable (Model 1).  This dataset contains  
#' 15 small areas. The response variable has three categories. The last 
#' is the reference category. 
#' The variables are as follows:
#'
#' \itemize{
#' \item area: area indicator.
#' \item Time: time indicator.
#' \item sample: the sample size of each domain.
#' \item Population: the population size of each domain.
#' \item Y1: the first category of the response variable.
#' \item Y2: the second category of the response variable.
#' \item Y3: the third category of the response variable.
#' \item X1: the covariate for the first category of the response variable.
#' \item X2: the covariate for the second category of the response variable.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name simdata
#' @usage simdata
#' @format A data frame with 15 rows and 9 variables in columns
#' @examples
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata)  #data
#' mod=1 # type of model
#' datar=data.mme(simdata,k,pp,mod)
#' # Model fit
#' result=model(datar$d,datar$t,pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],
#'        datar$n,datar$N,mod)
#'
#' #Analytic MSE
#' mse=msef(pp,datar$X,datar$Z,result,datar$N,datar$n)
#' 
#' B=1    #Bootstrap iterations
#' ss=12345 #SEED
#' set.seed(ss)
#'
#' ##Bootstrap parametric BIAS and MSE
#' mse.pboot=mseb(pp,datar$Xk,datar$X,datar$Z,datar$n,datar$N,result,B,mod)


NULL

#' Dataset for Model 2
#'
#' Dataset used by the multonomial mixed effects model with two independent random effects
#' in each category of the response variable: one domain random effect and another independent time and domain random effect (Model 2).
#' This dataset contains 10 small areas and two periods. The response variable has three categories. The last 
#' is the reference category.
#' The variables are as follows:
#'
#' \itemize{
#' \item area: area indicator.
#' \item Time: time indicator.
#' \item sample: the sample size of each domain.
#' \item Population: the population size of each domain.
#' \item Y1: the first category of the response variable.
#' \item Y2: the second category of the response variable.
#' \item Y3: the third category of the response variable.
#' \item X1: the covariate for the first category of the response variable.
#' \item X2: the covariate for the second category of the response variable.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name simdata2
#' @usage simdata2
#' @format A data frame with 30 rows and 9 variables in columns
#' @examples
#' 
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' data(simdata2)
#' mod=2 #type of model
#' datar=data.mme(simdata2,k,pp,mod)
#'
#' ##Model fit
#' result=model(datar$d,datar$t,pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],
#'        datar$n,datar$N,mod)
#'
#' ##Analytic MSE
#' msef=msef.it(pp,datar$X,result,datar$n,datar$N)
#'
#' B=1    #Bootstrap iterations
#' ss=12345 #SEED
#' set.seed(ss)
#'
#' ##Bootstrap parametric BIAS and MSE
#' mse.pboot=mseb(pp,datar$Xk,datar$X,datar$Z,datar$n,datar$N,result,B,mod)
NULL

#' Dataset for Model 3
#'
#' Dataset used by the multonomial mixed effects model with two independent random effects
#' in each category of the response variable: one domain random effect and another correlated time and domain random effect (Model 3).
#' This dataset contains ten small areas and four periods. The response variable has three categories. The last 
#' is the reference category.
#' The variables are as follows:
#'
#' \itemize{
#' \item area: area indicator.
#' \item Time: time indicator.
#' \item sample: the sample size of each domain.
#' \item Population: the population size of each domain.
#' \item Y1: the first category of the response variable.
#' \item Y2: the second category of the response variable.
#' \item Y3: the third category of the response variable.
#' \item X1: the covariate for the first category of the response variable.
#' \item X2: the covariate for the second category of the response variable.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name simdata3
#' @usage simdata3
#' @format A data frame with 40 rows and 9 variables in columns
#' @examples
#' \dontrun{
#'
#' k=3 #number of categories of the response variable
#' pp=c(1,1) #vector with the number of auxiliary variables in each category
#' mod=3 #type of model
#' data(simdata3) #data
#' datar=data.mme(simdata3,k,pp,mod)
#' ##Model fit
#' result=model(datar$d,datar$t,pp,datar$Xk,datar$X,datar$Z,datar$initial,datar$y[,1:(k-1)],
#'        datar$n,datar$N,mod)
#'
#' ##Analytic MSE
#' msef=msef.ct(pp,datar$X,result,datar$n,datar$N)
#'
#' B=1    #Bootstrap iterations
#' ss=12345 #SEED
#' set.seed(ss)
#'
#' ##Bootstrap parametric BIAS and MSE
#' mse.pboot=mseb(pp,datar$Xk,datar$X,datar$Z,datar$n,datar$N,result,B,mod)}
NULL

