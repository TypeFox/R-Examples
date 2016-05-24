#'  @title Cox regression for interval censored data using multiple imputation
#'  @author Marc Delord 
#'  @description Uses the multiple imputation approach to compute the regression coefficient and its associated
#'  variance-covariance matrix, and the baseline survival estimates of a Cox proportional hazards regression for interval censorded data
#'  @inheritParams MIICD.crreg
#'  @inheritParams MI.ci
#'  @examples
#'  res <- MIICD.coxph(formula = ~ treatment, k = 5, m = 5, data = bcos, verbose = FALSE)
#'  plot(res)
#'  #diagnostic plot for coefficients end associated standard error
#'  plot(res , type = 'coef' , coef = 1)
#'  plot(res , type = 'sigma' , coef = 1)
#'  @export
#'  @import survival MASS  
#'  @return \code{est} A data frame with estimates
#'  @details
#'  This function uses multiple imputation approach to estimate regression coefficient, its variance-covvariance 
#'   matrix, and baseline survival estimates for a Cox proportional hazards regression for interval censorded data.
#'  
#'  Estimates are computed using Rubin's rules (Rubin (1987)). Estimate of coefficient is computed as the mean of estimates over imputation. #'  The variance-covariance matrix is computed as the within imputation variance and the between imputation variance augmented by an
#'  inflation factor to take into account the finite number of imputation. At each iteration, the baseline survival function is updated
#'  and multiple imputation is performed using updated estimates.
#'  
#'  Print and plot methods are available to handle results.
#'  
#'  The \code{data} must contain at last two columns: \code{left} and \code{right}. For interval censored data, the \code{left} and the
#'  \code{right} columns indicates lower and upper bounds of intervals respectively. \code{Inf} in the right column stands
#'  for right censored observations.
#'  
#'  @references Delord, M. & Genin, E. Multiple Imputation for Competing Risks Regression with Interval Censored Data Journal of Statistical
#'  Computation and Simulation, 2015
#'  @references PAN, Wei. A Multiple Imputation Approach to Cox Regression with Interval-Censored Data. Biometrics, 2000, vol. 56, no 1,
#'   p. 199-203.
#'  @references Rubin, D. B. (1987). Multiple imputation for nonresponse in surveys. 
#'  @references Schenker, N. and Welsh, A. (1988). Asymptotic results for multiple imputation. The Annals of Statistics pages 1550-1566.
#'  @references Tanner, M. A. and Wong, W. H. (1987). An application of imputation to an estimation problem in grouped lifetime analysis.
#'  Technometrics 29, 23-32.
#'  @references Wei, G. C., & Tanner, M. A. (1991). Applications of multiple imputation to the analysis of censored regression data.
#'  Biometrics, 47(4), 1297-1309.
#'  @seealso \link[survival]{Surv}, \link[survival]{survfit}, \link[survival]{coxph}, \link[MASS]{mvrnorm}

MIICD.coxph<-function(formula , k , m , data , method = c('PMDA','ANDA') , verbose = FALSE )
{

if( class( formula ) != "formula"  ) stop('the formula argument must be a formula')
if( !is.numeric(k)  ) stop('k must be an integer')
if( !is.numeric(m)  ) stop('m must be an integer')
if( !is.data.frame(data)  ) stop('data must be a data.frame')
if( !is.logical(verbose)  ) stop('verbose must be logical')

  cl <- match.call()
  if(verbose){
  cat('\nCox regression for interval censored data using data augmentation and multiple imputation\n')
  cat("\nCall:\n", paste(deparse(cl), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat('Interval-censored Response of a Proportional Hazard model:\n\n')
  n<-nrow(data)
  cat('No.Observation:', n , '\n')
  cat('Patern:\n')
  stat<-ifelse(data$right == Inf ,'right-censored', 'event of interest')
  type<-ifelse(data$right == data$left , 'exact' , NA )
  type<-ifelse(data$right!=data$left & data$right!=Inf , 'interval-censored' , type )
  type<-ifelse(data$right == Inf  , 'right-censored' , type )
  print(table('Cause'=stat, type))
  }
  
  method<-match.arg(method)
  
  if(method=='PMDA')
  res<-PMDA.coxph( formula = formula , data = data , k = k , m = m )
  else
  if (method == 'ANDA')
  res<-ANDA.coxph( formula = formula , data = data , k = k , m = m )  
 
  beta  <- tail( t( res$beta ) , 1 )
  sigma <- tail( t( res$sigma ) , 1 ) 
  
  vcov <- res$vcov
  sigma_seq<- res$sigmac^.5
  #Compute the Pvalues
  pv<-1-pchisq((beta/(sigma^.5))**2 , df = 1)
  #print the results to the terminal
  df1<-data.frame('coef' = t(beta) ,'exp(coef)'= t(exp( beta )),'se(coef)'= t(sigma**.5),'z'=t(beta/sigma^.5),pv=t(pv),'.'= '')
  colnames(df1)<-c('coef','exp(coef)' , 'se(coef)' , 'z' , 'p', '')
  
  if(verbose){
  cat("Coefficients:\n")
  print(format(df1  ,  digits  =  max(3L ,  getOption("digits") - 3L))   ,  print.gap  =  3L  ,  quote  =  FALSE)
  cat('\n')}
    
  ret<-list('Coef.' =  beta ,  'vcov'  =  vcov  ,  'Coef_seq' = res$beta , 'Sigma_seq' = sigma_seq , s0 = res$s0 , df = df1 , method = method ,  data = data , call = cl)
  class(ret)<-'MIICD_coxph'
  return(ret)
}

#' @export
print.MIICD_coxph<-function( x , ... ){
    data<-x$data
    cl<-x$call
    n<-nrow(data)
    cat('\nCox regression for interval censored data using data augmentation and multiple imputation\n')
    cat('\nCall:\n', paste(deparse(cl), sep = '\n', collapse = '\n'), '\n\n', sep = '')
    cat('Interval-censored Response of a Proportional Hazard model:\n\n')
    cat('No.Observation:', n , '\n\n')
    cat('Patern:\n\n')
    stat<-ifelse(data$right == Inf ,'right-censored', 'event of interest')
    type<-ifelse(data$right == data$left , 'exact' , NA )
    type<-ifelse(data$right!=data$left & data$right!=Inf , 'interval-censored' , type )
    type<-ifelse(data$right == Inf  , 'right-censored' , type )
    print(table('Cause'=stat, type))
    cat('\nCoefficients:\n')
    print(format( x$df  ,  digits  =  max(3L ,  getOption('digits') - 3L))   ,  print.gap  =  3L  ,  quote  =  FALSE)
    cat('\n')     
  }

#' plot method for MIICD_coxph objects
#' @param x a MIICD_coxph object
#' @param type type of diagnoctic plot to display
#' @param coef An integer: the no of the coefficient to display
#' @inheritParams plot.MI_surv
#' @export
plot.MIICD_coxph<-function ( x , type = c('baseline','coef','sigma') , coef = 1 ,  ylab = 'Survival' , xlab = 'Time' , ... ) 
    {
      type<-match.arg(type)
      if(type == 'baseline'){
      s0<-x$s0
      plot( s0$time , s0$surv , ylab = ylab , xlab = xlab , bty = 'l' , main = x$call ,  type = 's' , ylim =c( 0 , 1))  
      mtext(side = 3 ,  '(final baseline survival)')
    }else if(type == 'coef'){
    plot(as.vector(x$Coef_seq[coef,]) , type ='b' ,xlab='Iteration' , ylab = rownames(x$Coef_seq)[coef])
    }else if(type == 'sigma'){
      plot(as.vector(x$Sigma_seq[coef,]) , type ='b' ,xlab='Iteration' , ylab = paste('sd(',rownames(x$Coef_seq)[coef],')',sep=' '))
    }
}