distr.fit <- function(object, distr){ 
  if(distr=="poisson"){
    distrcoefs <- NULL
    sigmasq <- 0
  }
  if(distr=="nbinom"){  
    fitval <- object$fitted.values
    ts <- object$response
    n <- object$n_eff
    m <- length(object$coefficients)  
    #Pearson type estimator:
    find_root <- function(v) sum((ts-fitval)^2/(fitval*(1+fitval/v))) - n + m
    root <- try(uniroot(f=find_root, interval=c(0, 1e100)), silent=TRUE)  
    if(class(root)=="try-error"){
      size <- NA
      warning("The parameter of the negative binomial distribution cannot be estimated.\nTry a Poisson distribution with argument 'distr' set to \"poisson\" instead")
    }else{
      size <- root$root 
    }
  #Moment estimator: (does not work reliably)  
  #   v_1 <- ((1/n) * sum(((ts-fitval)^2 - fitval)/(fitval^2)))^(-1)
  #   mom_est <- function(v){
  #     abs(sum((ts-fitval)^2/(fitval*(1+fitval/v))) - n + m)
  #   }  
  #   v_21 <- optim(par=1, fn=mom_est, method="Brent", lower=0, upper=5000)
  #   v_22 <- optimize(f=mom_est, interval=c(0,5000))
    distrcoefs <- c(size=size)
    sigmasq <- 1/size
  }
  result <- list(distr=distr, distrcoefs=distrcoefs, sigmasq=sigmasq)
  return(result)
}
