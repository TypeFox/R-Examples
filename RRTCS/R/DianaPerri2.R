#' @name DianaPerri2
#' @aliases DianaPerri2
#' @title Diana-Perri-2 model
#' 
#' @description Computes the randomized response estimation, its variance estimation and its confidence interval through the Diana-Perri-2 model.
#' The function can also return the transformed variable. 
#' The Diana-Perri-2 model was proposed by Diana and Perri (2010, page 1879).
#' 
#' @usage DianaPerri2(z,mu,beta,pi,type=c("total","mean"),cl,N=NULL,method="srswr")
#' @param z vector of the observed variable; its length is equal to \eqn{n} (the sample size)
#' @param mu vector with the means of the scramble variables \eqn{W} and \eqn{U}
#' @param beta the constant of weighting
#' @param pi vector of the first-order inclusion probabilities
#' @param type the estimator type: total or mean 
#' @param cl confidence level
#' @param N size of the population. By default it is NULL
#' @param method method used to draw the sample: srswr or srswor. By default it is srswr
#' 
#' @details 
#' In the Diana-Perri-2 model, each respondent is asked to report the scrambled response \eqn{z_i=W(\beta U+(1-\beta)y_i)} where \eqn{\beta \in [0,1)} is a suitable constant 
#' controlled by the researcher and \eqn{W,U} are scramble variables whose distribution is assumed to be known.      
#' 
#' To estimate \eqn{\bar{Y}} a sample of respondents is selected according to simple random sampling with replacement.
#' The transformed variable is
#' \deqn{r_i=\frac{z_i-\beta\mu_W\mu_U}{(1-\beta)\mu_W}}
#' where \eqn{\mu_W,\mu_U} are the means of \eqn{W,U} scramble variables, respectively.                 
#'
#' The estimated variance in this model is
#' \deqn{\widehat{V}(\widehat{\bar{Y}}_R)=\frac{s_z^2}{n(1-\beta)^2\mu_W^2}} 
#' where \eqn{s_z^2=\sum_{i=1}^n\frac{(z_i-\bar{z})^2}{n-1}}.
#' 
#' If the sample is selected by simple random sampling without replacement, the estimated variance is
#' \deqn{\widehat{V}(\widehat{\bar{Y}}_R)=\frac{s_z^2}{n(1-\beta)^2\mu_W^2}\left(1-\frac{n}{N}\right)} 
#' 
#' @return Point and confidence estimates of the sensitive characteristics using the Diana-Perri-2 model. The transformed variable is also reported, if required.
#' 
#' @references  Diana, G., Perri, P.F. (2010).
#' \emph{New scrambled response models for estimating the mean of a sensitive quantitative character.}
#' Journal of Applied Statistics 37 (11), 1875-1890.
#' 
#' @seealso \code{\link{DianaPerri2Data}}
#' @seealso \code{\link{DianaPerri1}}
#' @seealso \code{\link{ResamplingVariance}}
#' 
#' @keywords Randomized_response Quantitative DianaPerri Estimation Variance Transformed_variable Confidence_interval
#' 
#' @examples
#' N=100000
#' data(DianaPerri2Data)
#' dat=with(DianaPerri2Data,data.frame(z,Pi))
#' beta=0.8
#' mu=c(50/48,5/3)
#' cl=0.95
#' DianaPerri2(dat$z,mu,beta,dat$Pi,"mean",cl,N,"srswor")
#' 
#' @export
DianaPerri2=function(z,mu,beta,pi,type=c("total","mean"),cl,N=NULL,method="srswr"){
  if(!is.vector(z)){stop("z must be a vector.")}
  if(any(is.na(z))){stop("There are missing values in z.")}  
  n=length(z)
  
  if(!is.vector(mu)){stop("mu must be a vector.")}
  if(any(is.na(mu))){stop("There are missing values in mu.")}  
  if(length(mu)!=2){stop("The length of the vector mu must be two")}
  
  if((beta<0)|(beta>1)){stop("There are invalid values in beta.")}
  
  if(!is.vector(pi)){stop("pi must be a vector.")}
  if(any(is.na(pi))){stop("There are missing values in pi.")}
  if(any((pi<=0)|(pi>1))){stop("There are invalid values in pi.")}         
  
  if(n!=length(pi)){stop("The lengths of z and pi are different.")}
  
  if(!is.character(type)){stop("type must be a character")}
  if((type!="total")&(type!="mean")){stop("The value of type must be total or mean.")}  
  
  if((cl<=0)|(cl>=1)){stop("The value of cl must be in the range (0,1).")}
  
  if(!is.character(method)){stop("method must be a character")}
  if((method!="srswr")&(method!="srswor")){stop("The value of method must be srswr or srswor.")}  
  
  
  call=match.call()
  
  r=(z-(beta*mu[1]*mu[2]))/((1-beta)*mu[1])

  zalpha=qnorm(1-(1-cl)/2)
  if(type=="total"){
    e=sum(r/pi) 
    if(is.null(N)){
      warning("To calculate the estimated variance is needed the size of the population")
      out=list(Estimation=e)
      return(out)
    }else{
      if(length(N)!=1){stop("N must be a scalar.")}  
      if(N<0){stop("N must be a positive number.")}
      
      ve=N^2*(var(z)/(n*(1-beta)^2*mu[1]^2))
      if(method=="srswor"){
        ve=ve*(1-n/N)
      }
      ci=c((e-zalpha*sqrt(ve)),(e+zalpha*sqrt(ve)))
    }
  }
  
  if(type=="mean"){
    if(is.null(N)){
      N=sum(1/pi)
    }
    if(length(N)!=1){stop("N must be a scalar.")}  
    if(N<0){stop("N must be a positive number.")}
    
    e=(1/N)*sum(r/pi)
    ve=var(z)/(n*(1-beta)^2*mu[1]^2)
    if(method=="srswor"){
      ve=ve*(1-n/N)
    }
    ci=c((e-zalpha*sqrt(ve)),(e+zalpha*sqrt(ve)))
  }
  
  if(ve<0){warning("The variance estimation can not be negative.")}
  
  out1=list(TransformedVariable=r)
  out2=list(Estimation=e,Variance=ve,ConfidenceInterval=ci)
  out=c(Call=call,out1,out2,Name="Diana and Perri",Model="Quantitative",Type=type,Param=list(c("mu"=mu,"beta"=beta)),ConfidenceLevel=cl)
  
  class(out)="RRT"
  return(out)
}