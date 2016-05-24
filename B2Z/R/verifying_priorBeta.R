#######################################################
#This function checks whether the prior distribution  #
#for Beta is valid.                                   #
#######################################################

verifying_prior_Beta<-function(ind,a,b)
  {
  if(ind==1)
    {
    if(!exists("a")){stop("missing lower limit of the uniform prior distribution of Beta.")}
    if(!exists("b")){stop("missing upper limit of the uniform prior distribution of Beta.")}
    if(a<0) {stop("the lower limit of the uniform prior distribution of Beta must be greater than 0.")}
    if(a>=b) {stop(" the lower limit of the uniform prior distribution of Beta must be lower than the upper limit.")}
    }

  if(ind==2)
    {
    if(!exists("a")){stop("missing shape parameter of the gamma prior distribution of Beta.")}
    if(!exists("b")){stop("missing scale paramter  of the gamma prior distribution of Beta.")}
    if(a<=0) {stop("the shape parameter in the gamma prior distribution of Beta must be greater than 0.")}
    if(b<=0) {stop("the scale parameter in the gamma prior distribution of Beta must be greater than 0.")}
    }

  if(ind==3)
    {
    if(!exists("a")){stop("missing scale parameter of the exponential prior distribution of Beta.")}
    if(a<=0) {stop("the scale parameter in the exponential prior distribution of Beta must be greater than 0.")}
    }

  if(ind==4)
    {
    if(!exists("a")){stop("missing mean parameter of the normal prior distribution of Beta.")}
    if(!exists("b")){stop("missing standard deviation paramter  of the normal prior distribution of Beta.")}
    if(b<=0) {stop("the sigma parameter in the normal prior distribution of Beta must be greater than 0.")}
    }

  if(ind==5)
    {
    if(!exists("a")){stop("missing d.f parameter of the t prior distribution of Beta.")}
    if(a<=0) {stop("the d.f. in the t prior distribution of Beta must be greater than 0.")}
    }

  if(ind==6)
    {
    if(!exists("a")){stop("missing shape parameter of the weibull prior distribution of Beta.")}
    if(!exists("b")){stop("missing scale paramter  of the weibull prior distribution of Beta.")}
    if(a<=0) {stop("the shape parameter in the weibull prior distribution of Beta must be greater than 0.")}
    if(b<=0) {stop("the scale parameter in the weibull prior distribution of Beta must be greater than 0.")}
    }

  if(ind==7)
    {
    if(!exists("a")){stop("missing d.f. of the F prior distribution of Beta.")}
    if(!exists("b")){stop("missing d.f. of the F prior distribution of Beta.")}
    if(a<=0) {stop("the d.f. of the F prior distribution of Beta must be greater than 0.")}
    if(b<=0) {stop("the d.f. of the F prior distribution of Beta must be greater than 0.")}
    }

  if(ind==8)
    {
    if(!exists("a")){stop("missing d.f parameter of the chi-squared prior distribution of Beta.")}
    if(a<=0) {stop("the d.f. in the chi-squared prior distribution of Beta must be greater than 0.")}
    }

  if(ind==9)
    {
    if(!exists("a")){stop("missing location parameter of the cauchy prior distribution of Beta.")}
    if(!exists("b")){stop("missing scale paramter  of the cauchy prior distribution of Beta.")}
    if(b<=0) {stop("the scale parameter in the cauchy prior distribution of Beta must be greater than 0.")}
    }


  if(ind==10)
    {
    if(!exists("a")){stop("missing mean parameter of the log-normal prior distribution of Beta.")}
    if(!exists("b")){stop("missing standard deviation paramter  of the log-normal prior distribution of Beta.")}
    if(b<=0) {stop("the standard parameter in the log-normal prior distribution of Beta must be greater than 0.")}
    }

  }



 
  
   
  
  


