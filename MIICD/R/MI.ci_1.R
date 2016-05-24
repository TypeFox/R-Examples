#  @title Cumulative incidence estimation for interval censored competing risks data using multiple imputation
#  @author Marc Delord \email{<mdelord@@gmail.com>}
#  @description Uses multiple imputation to compute a cumulative incidence function for interval censored competing risks data
#  @inheritParams MI.ci
#  @param ntimes Number of time points where the estimates should be computed (experimental)
#  @examples
#  res <- MIICD:::MI.ci_1(m = 10, status = 'status', trans = 1, data = ICCRD,
#  conf.int = TRUE, cens.code = 0, alpha = 0.05)
#  res
#  plot(res)
#  @import survival
#  @return \code{est} A data frame with estimates
#  @return \code{\dots} Other objects
#
#  @details This function uses the multiple imputation approach to estimate the cumulative incidence function for interval censored
#  competing risks data.
#  
#  Estimates are computed using Rubin's rules (Rubin (1987)). Estimate of the cumulative incidence is computed as the mean of cumulative
#  incidences over imputations. The variance is computed at each point by combining the within imputation variance and the between
#  imputation variance augmented by an inflation factor to take into account the finite number of imputation. If \code{conf.inf} is 
#  required, the log-log transformation is used to compute the lower confidence interval.
#  
#  #Print and plot methods are available to handle results.
#  
#  The \code{data} must contain at last three columns: \code{left},  \code{right} and \code{status}. For interval censored data, the
#  \code{left} and \code{right} columns indicates the lower and the upper bounds of the intervals respectively. \code{Inf} in the
#  \code{right} column stands for right censored observations. When an observation is right censored, the \code{status} column must
#  contain the censor indicator specified by \code{cens.code}. The transition of interest must be specified by the \code{trans}
#  parameter.
#  
#    
#  @references PAN, Wei. A Multiple Imputation Approach to Cox Regression with Interval-Censored Data. Biometrics, 2000, vol. 56, no 1,
#   p. 199-203.
#  @references Rubin, D. B. (1987). Multiple imputation for nonresponse in surveys. 
#  @references Schenker, N. and Welsh, A. (1988). Asymptotic results for multiple imputation. The Annals of Statistics pages 1550-1566.
#  @references Tanner, M. A. and Wong, W. H. (1987). An application of imputation to an estimation problem in grouped lifetime analysis.
#  Technometrics 29, 23-32.
#  @references Wei, G. C., & Tanner, M. A. (1991). Applications of multiple imputation to the analysis of censored regression data.
#  Biometrics, 47(4), 1297-1309.
#  @seealso \link[survival]{Surv}, \link[survival]{survfit}  



MI.ci_1<-function( m , status ,  trans , data , conf.int = TRUE , cens.code , alpha = 0.05 , ntimes = NULL ){

if( !is.numeric(m)  ) stop('m must be an integer')
if( !is.data.frame(data)  ) stop('data must be a data.frame')
if( !is.logical(conf.int)  ) stop('conf.int must be logical')
if( alpha <= 0 | alpha >= 1 ) stop('alpha must be in ] 0 , 1 [')

cl<-match.call()  
#Use interval censored data and generate k sets of imputed data
sets<-sapply( 1:m , get.set , data )
#Get and sort single times at wich the cumulative incidence will be estimated
times<-as.vector(sets)
length(times)
  ###

r2 <- as.character( rep( data[ , status ] , m  ) )
r2<-factor(r2,levels=unique(c(cens.code,unique(r2))))

r1 <- as.character(  data[ , status ] ) 
r1<-factor(r1,levels=unique(c(cens.code,unique(r1))))

fitCI <- survfit( Surv( time = times , event = r2 , type = "mstate"  ) ~ 1  , weights = rep( 1 / m , length(times) ) ,
conf.type = 'none' )

w <- which( fitCI$states == trans )
pr <- fitCI$prev[ , w ]
sd <- fitCI$std.err[ , w ]
t0 <- fitCI$time
#get estimated of cumulative incidence and confidence intervals

if(! is.null(ntimes) ){ 
t1 <- seq( from = range(t0)[1] , to=range(t0)[2] , length = ntimes )
sd_at_times<-sapply( 1 , get_values_at_times , values =  list(sd) , times = list(t0) , at = t1 , list = T )
sd_at_times <- get_z( sd_at_times )
ci_at_times<-sapply( 1 , get_values_at_times , values =  list(pr) , times = list(t0) , at = t1 , list = T )
ci_at_times <- get_z( ci_at_times )
  
}else{
  t1 <- t0
  sd_at_times <- sd
  ci_at_times <- pr  
}
if(conf.int){
sap<-lapply( 1:m , get_est_mi , trans = trans , imp_sets = sets , data = data , cens.code = cens.code ,  r2 = r1 )
#obtain data frame of standard errors and point estimates
cis <- sapply(sap,function(x) x[['est']])
t3  <- sapply(sap,function(x) x[['time']])
#get standard errors and point estimates at single times
cis_at_times<-sapply( 1:m , get_values_at_times , values =  cis , times = t3 , at = t1 , list = F )
cis_at_times <- get_z( cis_at_times )
CI <- post_point_est_CI( beta = cis_at_times , sd = sd_at_times , times = t1 , conf.int =  conf.int , alpha = alpha )
CI <- unique(replace(CI , is.na(CI) , 0 ))
}else{
CI<-rbind(c(time = 0 ,  est = 0 ) , data.frame( time = t0 , est = pr ) )
CI <- unique(replace(CI , is.na(CI) , 0 ) )
}
ret<-list( est = CI , call = cl , data = data , cens.code = cens.code , status = status , conf.int = conf.int )
class(ret) <- 'MI_ci'  
return(ret)
}

print.MI_ci_1 <- function (x , ... ) {
  cat('\nCumulative incidence estimation for interval censored data using multiple imputation\n')
  cat( "\nCall:\n", paste( deparse(x$call) , sep = "\n" , collapse = "\n" ) , "\n\n" , sep = "")
  cat('Interval-censored response for cumulative incidence estimate :\n\n')
  n<-nrow(x$data)
  data<-x$data    
  cens <- x$cens.code
  status<-x$status
  cat('No.Observation:', n , '\n')
  cat('Patern:\n')
  stat<-ifelse(data[,status]==cens,'unknown (right-censored)',as.character(data[,status]))
  type<-ifelse(data$right==data$left , 'exact' , NA )
  type<-ifelse(data$right!=data$left & data$right!=Inf , 'interval-censored' , type )
  type<-ifelse(data[,status]==cens  , 'right-censored' , type )
  print(table('Cause'=stat, type))
  cat('\n')  
  cat('$est\n')
  dimest<-paste(dim(x$est)[1] , 'x' , dim(x$est)[2])
  cat(paste('A',dimest,'data frame of required estimates\n'))
  print(head(x$est))
    }


# plot method for MI_ci objects
# @param x A MI_ci object
# @inheritParams plot.MI_surv

plot.MI_ci_1 <- function (x , xlab = 'Time' , ylab = 'Cumulative incidence' , ... )
      {
      data <- x$est
      conf.int <- x$conf.int
      plot( data$time , data$est , xlab = xlab , ylab = ylab , type = 's' , ylim = c(0,1) , bty ='l')  
      if(conf.int){
      lines( data$time , data$uci , lty = 2 , type = 's' )
      lines( data$time , data$lci , lty = 2 , type = 's' )
      }
      }




