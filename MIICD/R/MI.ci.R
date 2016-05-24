#'  @title Cumulative incidence estimation for interval censored competing risks data using multiple imputation
#'  @author Marc Delord \email{<mdelord@@gmail.com>}
#'  @description Uses multiple imputation to compute the cumulative incidence function for interval censored competing risks data
#'  @param k An integer, indicates the number of iteration to perform
#'  @param m An integer, indicates the number of imputation to perform at each iteration
#'  @param status The name of the column where status are to be found
#'  @param trans Denomination of the event of interest in the status column
#'  @param data The input data (see details)
#'  @param conf.int Logical, computes the confidence interval
#'  @param cens.code Censor indicator in the status column of the data
#'  @param alpha Parametrize the confidence interval width
#'  @examples
#'  res <- MI.ci(k = 5,  m = 5, status = 'status',  trans = 1 , data = ICCRD,
#'  conf.int = TRUE, cens.code = 0 , alpha = 0.05)
#'  res
#'  print(res)
#'  plot(res)
#'  @export
#'  @import survival
#'  @return \code{est} A data frame with estimates
#'  @return \code{\dots} Other objects
#'  @details This function uses a multiple imputation approach to estimate a cumulative incidence function for interval censored competing
#'  risks data.
#'  Estimates are computed using Rubin's rules (Rubin (1987)). The cumulative incidence is computed as the mean of 
#'  cumulative incidences over imputations. The variance is computed at each point by combining the within imputation variance and the
#'  between imputation variance augmented by an inflation factor to take into account the finite number of imputations.
#'  At each iteration, the cumulative incidence is updated and multiple imputation is performed using the updated estimate.
#'  If \code{conf.int} is required, the log-log transformation is used to compute the lower confidence interval.
#'  
#'  Print and plot methods are available to handle results.
#'  
#'  The \code{data} must contain at last three columns: \code{left}, \code{right} and \code{status}. For interval censored data, the
#'  \code{left} and \code{right} columns indicates lower and upper bounds of intervals, respectively. \code{Inf} in the
#'  \code{right} column stands for right censored observations. When an observation is right censored, the \code{status} column must
#'  contain the censor indicator specified by \code{cens.code}. The transition of interest must be specified by the \code{trans}
#'  parameter.
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
#'  @seealso \link[survival]{Surv}, \link[survival]{survfit}
  
MI.ci<-function( k , m , data , status , trans , cens.code , conf.int = F , alpha = 0.05 ){


if( !is.numeric(k)  ) stop('k must be an integer')
if( !is.numeric(m)  ) stop('m must be an integer')
if( !is.data.frame(data)  ) stop('data must be a data.frame')
if( !is.logical(conf.int)  ) stop('conf.int must be logical')
if( alpha <= 0 | alpha >= 1 ) stop('alpha must be in ] 0 , 1 [')
#if(k <= 1) stop('You may consider the MI.ci function')  
  
cl<-match.call()    
prep     <- preproc.crreg( data = data , m = m , trans = trans , status = status , cens.code = cens.code )
data_int <- prep$data2
data_fix <- prep$data1
or       <- prep$or
I        <- prep$I

r2 <- as.character( rep( data[ , status ] , m  ) )
r2 <- replace( r2 , r2 == cens.code ,  0 )

r1 <- as.character(  data[ , status ]  ) 
r1 <- replace( r1 , r1 == cens.code ,  0 )
  
#Multiple Imputation
CI <- MI.ci_1( m = m , status = status , trans = trans , cens.code = cens.code,
data = data , conf.int = F , alpha = alpha , ntimes = NULL )$est
CI$diff <- c(0 , diff( CI$est ) )

for(i in 1:k){

ss1<-apply(data_int , 1 , function(x ) subset( CI , time >=  as.numeric(x['left']) & time <= as.numeric(x['right']) ) )
tk2<-lapply(seq_len(nrow(data_int)) ,function(X)  ss1[[X]]$time)
  
samples<-t( sapply( seq_len(nrow(data_int)) , function(X) {
pk2 <- ss1[[ X ]]$diff  
sapply( 1:m , function(x){
if( sum( pk2 ) & length( pk2 ) > 1 ) sample(  tk2[[ X ]] , size = 1 , prob = pk2 ) 
else   mean( tk2[[ X ]] ) } )  } ) )

samples2<-rbind(samples,data_fix)[or,]
times<-as.vector(samples2)

ci<-Surv( time = times , event = r2 , type = 'mstate')
fitCI<-survfit( ci ~ 1 , weights = rep( 1 , length( times ) ) / m , conf.type = 'none')  
w <- which( fitCI$states == trans )
sd <- fitCI$std.err[ , w ]
pr <- fitCI$prev[ , w ]
t0 <- fitCI$time
CI<-unique(rbind(c(time = 0 ,  est = 0 ) , data.frame( time = t0 , est = pr ) ))
CI$diff <- c(0 , diff( CI$est ) )
}

sap<-lapply( 1:m , get_est_mi , trans = trans , imp_sets = samples2  , data = data , r2 = r1 )
#obtain data frame of standard errors and point estimates
cis <- sapply(sap,function(x) x[['est']])
t3  <- sapply(sap,function(x) x[['time']])
#get standard errors and point estimates at single times
  
cis_at_times<-sapply( 1:m , get_values_at_times , values =  cis , times = t3 , at = t0 , list = is.list(cis) )
cis_at_times <- get_z( cis_at_times )
CI <- post_point_est_CI( beta = cis_at_times , sd = sd , times = t0 , conf.int =  conf.int , alpha = alpha )

if(conf.int){
colnames(CI)<-c('time','prev','sd','uci','lci')
CI <- unique(replace(CI , is.na(CI) , 0 ))
}else{
colnames(CI)<-c('time','prev')
CI <- unique(replace(CI , is.na(CI) , 0 ))
}
ret<-list( est = CI , call = cl , data = data , cens.code = cens.code , status = status , conf.int = conf.int )
class(ret) <- 'MI_ci'  
return(ret)
}

#' @export
print.MI_ci <- function (x , ... ) 
    {
  cat('\nCumulative incidence estimation for interval censored data using data augmentation and multiple imputation\n')
  cat( "\nCall:\n", paste( deparse(x$call) , sep = "\n" , collapse = "\n" ) , "\n\n" , sep = "")
  cat('Interval-censored response for cumulative incidence estimation:\n\n')
  n<-nrow(x$data)
  data<-x$data    
  cens <- x$cens.code
  status <- x$status
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

#' plot method for MI_ci objects
#' @param x A MI_ci object
#' @inheritParams plot.MI_surv
#' @export

plot.MI_ci <- function ( x, xlab = 'Time', ylab = 'Cumulative incidence' , ... )
      {
      data <- x$est
      conf.int = x$conf.int
      plot( data$time , data$prev , xlab = xlab , ylab = ylab , type = 's' , ylim = c(0,1) , bty ='l')  
      if(conf.int){
      lines( data$time , data$uci , lty = 2 , type = 's' )
      lines( data$time , data$lci , lty = 2 , type = 's' )
      }
      }