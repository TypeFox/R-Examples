#' Tests the Coefficients of High Dimensional Generalized Linear Models
#'
#' Tests for whole or partial regression coefficient vectors for high dimensional generalized linear models.  
#'  
#' @param Y a vector of observations of length \eqn{n},  where \eqn{n} is the sample size.
#' @param X a design matrix with \eqn{n} rows and \eqn{p} columns, where \eqn{p} is the dimension of the covariates.
#' @param beta_0 a vector with length \eqn{p}. It is the value of regression coefficient under the null hypothesis in global test.  The default
#'               is \eqn{\beta_0=0} and it can be non-zero in the global test. In the test with nuisance coefficients, we only deal with  \eqn{\beta_0^{(2)}=0}. 
#' @param model a character string to describe the model and link function. The default is \code{"gaussian"}, which denotes the linear model using identity link.
#'              The other options are \code{"poisson"}, \code{"logistic"}  and \code{"negative_binomial"} models, where the poisson and negative binomial models using log link.                    
#' @param nuisance an index indicating which coefficients are nuisance parameter. The default is \code{"NULL"} (the global test). 
#' @return An object of class "HDGLM_test" is a list containing the following components:
#' \item{test_stat}{the standardized test statistic}
#' \item{test_pvalue}{pvalue of the test against the null hypothesis} 
#' @note In global test,  the function \code{"HDGLM_test"} can deal with the null hypothesis with non-zero coefficients (\eqn{\beta_0}). However, in test with nuisance coefficient,
#' the function can only deal with the null hypothesis with zero coefficients (\eqn{\beta_0^{(2)}}) in this version.
#' 
#' @author Bin Guo
#' @references Guo, B. and Chen, S. X. (2015). Tests for High Dimensional Generalized Linear Models.
#' @examples
#' ## Example: Linear model 
#' ## Global test: if the null hypothesis is true (beta_0=0)
#' alpha=runif(5,min=0,max=1)
#' ## Generate the data
#' DGP_0=DGP(80,320,alpha) 
#' result=HDGLM_test(DGP_0$Y,DGP_0$X)
#' ## Pvalue
#' result$test_pvalue
#'
#' ## Global test: if the alternative hypothesis is true 
#' ## (the square of the norm of the first 5 nonzero coefficients to be 0.2)
#' ## Generate the data
#' DGP_0=DGP(80,320,alpha,sqrt(0.2),5) 
#' result=HDGLM_test(DGP_0$Y,DGP_0$X) 
#' ## Pvalue
#' result$test_pvalue
#'
#' ## Test with nuisance coefficients: if the null hypothesis is true (beta_0^{(2)}=0)
#' ## The first 10 coefficients to be the nuisance coefficients 
#' betanui=runif(10,min=0,max=1)
#' ## Generate the data
#' DGP_0=DGP(80,320,alpha,0,no=NA,betanui) 
#' result=HDGLM_test(DGP_0$Y,DGP_0$X,nuisance=c(1:10)) 
#' ## Pvalue
#' result$test_pvalue
#'
#' ## Test with nuisance coefficients: if the alternative hypothesis is true 
#' ## (the square of the norm of the first 5 nonzero coefficients in beta_0^{(2)} to be 2)
#' ## The first 10 coefficients to be the nuisance coefficients 
#' betanui=runif(10,min=0,max=1)
#' ## Generate the data
#' DGP_0=DGP(80,330,alpha,sqrt(2),no=5,betanui) 
#' result=HDGLM_test(DGP_0$Y,DGP_0$X,nuisance=c(1:10)) 
#' ## Pvalue
#' result$test_pvalue



HDGLM_test<-function(Y,X,beta_0=NULL,nuisance=NULL,model="gaussian")
{
  
  if(is.null(nuisance)) # Global test
  {
    
    n=as.integer(nrow(X))
    p=as.integer(ncol(X))
    if(is.null(beta_0))
    {
      beta_0=rep(0,p)
    }   
    
    if(model=="gaussian")
    {
      res=Y-X%*%beta_0
    }
    if(model=="poisson")
    {
      res=Y-exp(X%*%beta_0)
    }
    if(model=="logistic")
    { 
      Y1=exp(X%*%beta_0)
      res=Y-Y1/(1+Y1)
    }
    if(model=="negative_binomial")
    {
      res=Y-exp(X%*%beta_0)
    }
    storage.mode(X)<-"double"
    storage.mode(res)<-"double"
    
   # dyn.load("C:/Users/bin/Desktop/Rpackage/Teststat.dll")
    result<-.Fortran("Teststat",X,n,p,res,z=numeric(1))
    stat=result$z
    mypvalue=1-pnorm(stat)  
	output=list(test_stat=stat,test_pvalue=mypvalue)
    return(output)
  }
  if(!is.null(nuisance)) # Nuisance test
  { 
    options( warn = -1 ) 
    n=as.integer(nrow(X))
    p=as.integer(ncol(X))
    p_1=length(nuisance)
    p_2=p-p_1
    Xnui<-X[,nuisance]
    if(model=="gaussian")
    {
      glm1<-glm.fit(Xnui,Y,intercept=FALSE,family=gaussian())  
      fitted<-as.numeric(glm1$fitted.values)
      res<-Y-glm1$fitted.values
    }
    if(model=="poisson")
    {
      glm1<-glm.fit(Xnui,Y,intercept=FALSE,family=poisson())  
      fitted<-as.numeric(glm1$fitted.values)
      res<-Y-glm1$fitted.values
    }
    if(model=="logistic")
    { 
      glm1<-glm.fit(Xnui,Y,intercept=FALSE,family=binomial())  
      fitted<-as.numeric(glm1$fitted.values)
      res<-Y-glm1$fitted.values
    }
    if(model=="negative_binomial")
    {
      glm1<-glm.fit(Xnui,Y,intercept=FALSE,family=poisson()) 
      fitted<-as.numeric(glm1$fitted.values)
      res<-Y-glm1$fitted.values
    }
    Xnew=X[,-nuisance]
    storage.mode(Xnew)<-"double"
    storage.mode(res)<-"double"
    storage.mode(p_2)<-"integer"
    #dyn.load("C:/Users/bin/Desktop/Rpackage/Teststat.dll")
    result<-.Fortran("Teststat",Xnew,n,p_2,res,z=numeric(1))
    stat=result$z
    mypvalue=1-pnorm(stat)  
    output=list(test_stat=stat,test_pvalue=mypvalue)
    return(output)
    options( warn = 0) 
  }
  
}  