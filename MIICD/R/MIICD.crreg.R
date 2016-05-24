#'  @title Fine & Gray regression for interval censored competing risks data using multiple imputation
#'  @author Marc Delord \email{<mdelord@@gmail.com>}
#'  @description Uses the multiple imputation approach to compute regression coefficient and its associated
#'  variance-covariance matrix, and baseline cumulative incidence estimates for interval censorded competing risks data
#'  @param formula A formula. The right hand side indicates names of covariables to be found in \code{data}
#'  @param verbose Logical, display the results ? 
#'  @inheritParams MI.ci
#'  @param method Which data augmentation scheme shall be used ? Two algorithms are implemented : \emph{The Poor man's Data Augmentation 
#'  scheme} and the \emph{Asymptotic Normal Data Augmentation scheme} (the later may be preferred).
#'  @export
#'  @import survival MASS mstate
#'  @examples
#'  res <- MIICD.crreg(formula = ~ treatment, k = 5, m = 5, status = 'status',
#'  trans = 1, data = ICCRD,  cens.code = 0, method = 'ANDA', verbose = FALSE )
#'  res
#'  plot(res)
#'  #diagnostic plot for coefficients end associated standard error
#'  plot(res , type = 'coef' , coef = 1)
#'  plot(res , type = 'sigma' , coef = 1)
#'  @return \code{Coef.} Final estimate of the coefficient
#'  @return \code{vcov}  Final estimate of the variance-covariance matrix
#'  @return \code{Coef_seq} Sequence of the coefficient estimate over iterations
#'  @return \code{Sigma_seq} Sequence of the coefficient standard deviation over iterations
#'  @return \code{df} data frame containing the main results 
#'  @return \code{\dots} Other returned values
#'  @details This function uses data augmentation and multiple imputation aproach to estimate regression coefficient, variance-covariance
#'  matrix and baseline cumulative incidence estimates in a competing risks proportional hazards regression model for interval censorded
#'  competing risks data.
#'  
#'  Estimates are computed using Rubin's rules (Rubin (1987)). Estimate of coefficient is computed as the mean of estimates over
#'   imputation. The variance-covariance matrix is computed as the within imputation variance and the between imputation variance
#'   augmented byan inflation factor to take into account the finite number of imputation. At each iteration, the baseline cumulative 
#'   incidence function is updated  and multiple imputation is performed using the updated estimates. Print and plot methods are
#'   available to handle results.
#'  
#'  \code{Print} and \code{plot} methods are available to handle results.
#'  
#'  The \code{data} must contain at last four columns. One named \code{left}, one named \code{right}, the name of the 3^{rd} is indicated
#'  by the \code{status} parameter and one for the covariate to be tested.  For interval censored data, the left and right columns
#'  indicates the lower and the upper bounds of the intervals respectively. \code{Inf} in the right column stands for right censored
#'  observations. When an observation is right censored, the \code{status} column must contain the censor indicator specified by
#'  \code{cens.code}. The transition of interest must be precised by the \code{trans} parameter. 
#'  
#'  @references Delord, M. & Genin, E. Multiple Imputation for Competing Risks Regression with Interval Censored Data Journal of Statistical
#'  Computation and Simulation, 2015
#'  @references Fine JP and Gray RJ (1999) A proportional hazards model for the subdistribution of a competing risk. JASA 94:496-509.
#'  @references PAN, Wei. A Multiple Imputation Approach to Cox Regression with Interval-Censored Data. Biometrics, 2000, vol. 56, no 1,
#'   p. 199-203.
#'  @references Rubin, D. B. (1987). Multiple imputation for nonresponse in surveys. 
#'  @references Schenker, N. and Welsh, A. (1988). Asymptotic results for multiple imputation. The Annals of Statistics pages 1550-1566.
#'  @references Tanner, M. A. and Wong, W. H. (1987). An application of imputation to an estimation problem in grouped lifetime analysis.
#'  Technometrics 29, 23-32.
#'  @references Wei, G. C., & Tanner, M. A. (1991). Applications of multiple imputation to the analysis of censored regression data.
#'  Biometrics, 47(4), 1297-1309.
#'  @seealso \link[survival]{Surv}, \link[survival]{survfit}, \link[riskRegression]{FGR}, \link[MASS]{mvrnorm}
    
MIICD.crreg<-function(formula, k, m, status, trans, cens.code, data, method = c('PMDA','ANDA'), verbose = FALSE )
{
  
if( class( formula ) != "formula"  ) stop('the formula argument must be a formula')
if( !is.numeric(k)  ) stop('k must be an integer')
if( !is.numeric(m)  ) stop('m must be an integer')
if( !is.data.frame(data)  ) stop('data must be a data.frame')
if( !is.logical(verbose)  ) stop('verbose must be logical')
  
  cl <- match.call()
  if(verbose){
  cat('\nFine & Gray regression for competing risks interval censored data using data augmentation and multiple imputation\n')
  cat("\nCall:\n", paste(deparse(cl), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat('Interval-censored Response of a proportional hazard model with competing risks:\n\n')
  x<-nrow(data)
  cat('No.Observation:', x , '\n')
  cat('Patern:\n')
  stat<-ifelse(data[,status]==cens.code,'unknown (right-censored)',as.character(data[,status]))
  type<-ifelse(data$right==data$left , 'exact' , NA )
  type<-ifelse(data$right!=data$left & data$right!=Inf , 'interval-censored' , type )
  type<-ifelse(data[,status]==cens.code  , 'right-censored' , type )
  print(table('Cause'=stat, type))
  }
  
  method <- match.arg(method)
  if( method == 'PMDA' )
  res<-PMDA.crreg( formula , data , k , m , status , trans , cens.code )
  else
  if( method == 'ANDA' )
  res<-ANDA.crreg( formula , data , k , m , status , trans , cens.code )

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
    
  ret<-list('Coef.' =  beta ,  'vcov'  =  vcov  ,  'Coef_seq' = res$beta , 'Sigma_seq' = sigma_seq , CI0 = res$ci0 , df = df1 , method = method , cens.code = cens.code , status = status , data = data , call = cl)
  class(ret)<-'MIICD_crreg'
  return(ret)
}

#' @export
print.MIICD_crreg <- function (x , ... ) 
    {
    data<-x$data
    cl<-x$call
    cens <- x$cens.code
    status <- x$status
    cat('\nFine & Gray regression for competing risks interval censored data using data augmentation and multiple imputation\n')
    cat("\nCall:\n", paste(deparse(cl), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat('Interval-censored Response of a proportional hazard model with competing risks:\n\n')
  n<-nrow(data)
  cat('No.Observation:', n , '\n\n')
  cat('Patern:\n\n')
  stat<-ifelse(data[,status]==cens,'unknown (right-censored)',as.character(data[,status]))
  type<-ifelse(data$right==data$left , 'exact' , NA )
  type<-ifelse(data$right!=data$left & data$right!=Inf , 'interval-censored' , type )
  type<-ifelse(data[,status]==cens  , 'right-censored' , type )
  print(table('Cause'=stat, type))
  cat("\nCoefficients:\n")
  print(format( x$df  ,  digits  =  max(3L ,  getOption("digits") - 3L))   ,  print.gap  =  3L  ,  quote  =  FALSE)
  cat('\n')     
  }

#' plot method for MIICD_crreg objects
#' @param x a MIICD_crreg object
#' @inheritParams plot.MIICD_coxph
#' @export
plot.MIICD_crreg <- function (x, type = c('baseline','coef','sigma'), coef = 1,  ylab = 'Cumulative incidence', xlab = 'Time', ... ) 
    {
      type<-match.arg(type)
      if(type == 'baseline'){
      CI0<-x$CI0
      plot( CI0$time , CI0$est , ylab = ylab , xlab = xlab , bty = 'l' , main = x$call ,  type = 's' , bty = 'l' , ylim = c( 0 , 1) )  
      mtext(side = 3 ,  '(final baseline cumulative incidence)')
    }else if(type == 'coef'){
    plot(as.vector(x$Coef_seq[coef,]) , type ='b', xlab='Iteration', ylab = rownames(x$Coef_seq)[coef])
    }else if(type == 'sigma'){
      plot(as.vector(x$Sigma_seq[coef,]) , type ='b' ,xlab='Iteration' , ylab = paste('sd(',rownames(x$Coef_seq)[coef],')',sep=' '))
    }
    }
