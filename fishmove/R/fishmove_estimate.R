fishmove.estimate <- function(data=NA,start=NA,ci=FALSE,rep=100,conf=0.95,...){
  
  # Check if movement data are displacements provided as a single vector (only absolute movement considered)
  if(missing(data)) stop("Please provide displacement/movement data as single numeric vector")
  if(!is.numeric(data)) stop("Please provide displacement/movement data as single numeric vector")
  
  # Get starting values for iteration/optimization
  if(missing(start)) starting_values <- list(sigma_stat=quantile(abs(data),0.1),sigma_mob=quantile(abs(data),0.85),p=0.67)
  else{
    if(!is.list(start)) stop("Please provide named list of starting values for optimization: start=list(sigma_stat=,sigma_mob=,p=)")
    if(!(names(start)[1]=="sigma_stat" && names(start)[2]=="sigma_mob" && names(start)[3]=="p")) stop("Please provide named list of starting values for optimization: start=list(sigma_stat=,sigma_mob=,p=)")
    else  starting_values <- start
  }
  
  # Definition of probability density function based on two superimposed normal distributions
  ddoublenorm <- function(x,sigma_stat,sigma_mob,p){
    dnorm(x,mean=0,sd=sigma_stat)*p+dnorm(x,mean=0,sd=sigma_mob)*(1-p)}
  
  
  if(!ci){
    # Create dataset for negative and positive data
    data_prepared <- c(abs(data),-abs(data))
    res_fit <- fitdistr(x=data_prepared,
                        densfun=ddoublenorm,
                        start=starting_values,
                        method="L-BFGS-B",
                        lower=c(0.0001,0.0001,0.00001),
                        upper=c(Inf,Inf,0.99999))
    
    out <- res_fit
    
  }
  
  
  if(ci){
    warning("Calculation of confidence intervals is currently under developement. Unstable results might occur")
    # Function to fit/optimize leptkurtic dispersal kernel
    fit_doublenorm_boot <- function(data,indices,starting_values){
      
      
      # Create dataset for negative and positive data
      data_prepared <- c(abs(data[indices]),-abs(data[indices]))
      res_fit <- fitdistr(x=data_prepared,
                          densfun=ddoublenorm,
                          start=starting_values,
                          method="L-BFGS-B",
                          lower=c(0.0001,0.0001,0.00001),
                          upper=c(Inf,Inf,0.99999))
      
      res_fit$estimate
    }
    boot_res_fit <- boot(data,fit_doublenorm_boot,starting_values=starting_values,R=rep)
    
    sigma_stat <- boot.ci(boot_res_fit,conf=conf,type="bca",index=1)$bca[c(4,5)]
    sigma_mob <- boot.ci(boot_res_fit,conf=conf,type="bca",index=2)$bca[c(4,5)]
    p <- boot.ci(boot_res_fit,conf=conf,type="bca",index=3)$bca[c(4,5)]
    out <- cbind(sigma_stat,sigma_mob,p)
    out <- rbind(boot_res_fit$t0,out)
    row.names(out) <- c("fit","lwr","upr")
  }

  
  return(out)
}

