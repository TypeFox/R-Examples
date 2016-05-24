#  @title Survival estimation for interval censored data using multiple imputation
#  @author Marc Delord \email{<mdelord@@gmail.com>}
#  @description  Uses multiple imputation to compute the survival function with interval censored data
#  @inheritParams MI.ci
#  @examples
#   res<-MIICD:::MI.surv_1( m = 10 , data = ICCRD , conf.int = TRUE )
#   res
#   plot( res )
#   #plot( res , fun = 'event')
#  @import survival
#  @return \code{est} A data frame with estimates
#  @details This function uses the multiple imputation aproach to estimate the survival function for interval censored data.
#  
#  Estimates arecomputed using Rubin's rules (Rubin (1987)). Survival function estimate is computed as the mean of survival over
#  imputations. The variance is computed at each point by combining the within imputation variance and the between imputation variance
#  augmented by an inflation factor to take into account the finite number of imputation. If \code{conf.inf} is required, the log-log
#  transformation is used to compute the lower confidence interval.   
#    
#  Print and plot methods are available to handle results.
#  
#  The \code{data} must contain at last two columns: \code{left} and \code{right}. For interval censored data, the \code{left} and the
#  \code{right} columns indicates lower and upper bounds of intervals respectively. \code{Inf} in the \code{right} column stands
#  for right censored observations.
#  
#  @references Delord, M. & Genin, E. Multiple Imputation for Competing Risks Regression with Interval Censored Data Journal of Statistical
#  Computation and Simulation, 2015
#  @references PAN, Wei. A Multiple Imputation Approach to Cox Regression with Interval-Censored Data. Biometrics, 2000, vol. 56, no 1,
#   p. 199-203.
#  @references Rubin, D. B. (1987). Multiple imputation for nonresponse in surveys. 
#  @references Schenker, N. and Welsh, A. (1988). Asymptotic results for multiple imputation. The Annals of Statistics pages 1550-1566.
#  @references Tanner, M. A. and Wong, W. H. (1987). An application of imputation to an estimation problem in grouped lifetime analysis.
#  Technometrics 29, 23-32.
#  @references Wei, G. C., & Tanner, M. A. (1991). Applications of multiple imputation to the analysis of censored regression data.
#  Biometrics, 47(4), 1297-1309.
#  @seealso \link[survival]{Surv}, \link[survival]{survfit}  


MI.surv_1<-function( m ,  data , conf.int = TRUE , alpha = 0.05 ){

if( !is.numeric(m)  ) stop('m must be an integer')
if( !is.data.frame(data)  ) stop('data must be a data.frame')
if( !is.logical(conf.int)  ) stop('conf.int must be logical')
if( alpha <= 0 | alpha >= 1 ) stop('alpha must be in ] 0 , 1 [')

cl<-match.call()  
#Use interval censored data and generate k sets of imputed data
sets<-sapply( 1:m , get.set , data = data  )
#Get and sort single times at wich the cumulative incidence will be estimated
times<-as.vector(sets)
surv<-Surv( time = times , event = rep( data$right != Inf , m ) , type = "right"  )
surv2<-Surv( time = times , event = rep( data$right != Inf , m ) , type = "mstate"  )  
surv2[,2]<-surv[,2]
fitsurv<-survfit( surv2 ~ 1 , weights = rep( 1 , length( times ) ) / m , conf.type = 'none')  
sd <- fitsurv$std.err
pr <- fitsurv$prev
t0 <- fitsurv$time
#get estimated of cumulative incidence and confidence intervals
if(conf.int){
sap<-lapply( 1:m , get_est_mi_surv , imp_sets = sets , data = data )
#obtain data frame of standard errors and point estimates
cis <- sapply(sap,function(x) x[['est']])
t3  <- sapply(sap,function(x) x[['time']])
#get standard errors and point estimates at single times
cis_at_times<-sapply( 1:m , get_values_at_times , values =  cis , times = t3 , at = t0 , list = F )
cis_at_times <- get_z( cis_at_times )
CI <- post_point_est_CI( beta = cis_at_times , sd = sd , times = t0 , conf.int =  conf.int , alpha = alpha )
CI$est<-1-CI$est
CI$lci<-1-CI$lci
CI$uci<-1-CI$uci
colnames(CI)<-c('time','surv','sd','uci','lci')
CI <- unique(replace(CI , is.na(CI) , 1 ))
}else{
CI<-rbind(c(time = 0 ,  surv = 1 ) , data.frame( time = t0 , surv = 1 - pr ))
CI <- unique(replace(CI , is.na(CI) , 1 ) )
}
ret<-list( est = CI , call = cl , data = data , conf.int = conf.int )
class(ret) <- 'MI_surv'  
return(ret)  
}


print.MI_surv_1 <- function (x , ... ) {
  data<-x$data
  cl<-x$call
  cat('\nSurvival estimation for interval censored data using multiple imputation\n')
  cat("\nCall:\n", paste(deparse(cl), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat('Interval-censored response for survival estimate :\n\n')
  n<-nrow(data)
  cat('No.Observation:', n , '\n\n')
  cat('Patern:\n\n')
  stat<-ifelse(data$right == Inf ,'right-censored', 'event of interest')
  type<-ifelse(data$right == data$left , 'exact' , NA )
  type<-ifelse(data$right!=data$left & data$right!=Inf , 'interval-censored' , type )
  type<-ifelse(data$right == Inf  , 'right-censored' , type )
  print(table('Cause'=stat, type))
  cat('\n')  
  cat('$est\n')
  dimest<-paste(dim(x$est)[1] , 'x' , dim(x$est)[2])
  cat(paste('A',dimest,'data frame of required estimates\n'))
  print(head(x$est))
  }

# plot method for MI_surv objects
# @param x A MI_surv object
# @inheritParams plot.MI_surv

plot.MI_surv_1 <- function (x , xlab = 'Time' , ylab = 'Survival' ,  fun = 'surv' , ... )
      {
        data <- x$est
        conf.int <- x$conf.int
        if( fun == 'event' ){
          plot( data$time , 1-data$surv , xlab = xlab , ylab = 'Prevalence' , type = 's' , ylim = c(0,1) ,bty = 'l' )  
          if(conf.int){
            lines( data$time , 1-data$uci , lty = 2 , type = 's' )
            lines( data$time , 1-data$lci , lty = 2 , type = 's' )
            }
            }else{
            plot( data$time , data$surv , xlab = xlab , ylab = ylab , type = 's' , ylim = c(0,1) , bty = 'l' )  
            if(conf.int){
            lines( data$time , data$uci , lty = 2 , type = 's' )
            lines( data$time , data$lci , lty = 2 , type = 's' )
            }
          }
      }