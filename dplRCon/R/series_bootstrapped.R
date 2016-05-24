#' Series bootstrap
#' 
#' This function calculate the bootstrapped replicated series for calculating concordance  
#' @param	data.boot Data matrix containing standardized ring width indices from which the bootstrapped replicates are calculated. 
#' @param	stat A function giving the statistics that are to be calculated
#' @param	R The number of bootstrapped replicates
#' @param	names.stat the names of the statistics contained in function "stat"
#' @param	aver.by.tree.input This is a True/False. If True then averages tree series, then average tree means. If False then average all series 
#' @return \item{org.stat}{The original input statistics}
#' \item{boot.stat}{The statistics of the bootstrapped replicates, c('mean','median')}
#' \item{bias}{The bias of the boot.stat}
#' \item{std.error}{The standard error of boot.stat}
#' \item{boot.series.mean}{The bootstrapped replicates mean series}	
#' \item{boot.series.std}{The bootstrapped replicates standard deviation series}
#' @examples 
#'\dontrun{#Subset "near-pith" is the material within 0 -20cm from the estimated pith
#' spline200.sub0.20.n   <- TruncSeriesPithoffset( ring.raw, ring.stand, dbh.po.nc, c(1,200))
#' boot.0.20   <-  series.bootstrap( spline200.sub0.20.n$sub.series.stand, stat, 999, 
#'    names.stat, aver.by.tree = FALSE)}
#' @export
series.bootstrap	<- function( data.boot, stat, R, names.stat, aver.by.tree.input ) { 
  #Setup
  N.stat 		<- length( names.stat )
  boot.stat	<- matrix( NA, R, N.stat, dimnames = list( 1:R, names.stat ) )
  boot.series.mean <-matrix( NA, nrow = dim( data.boot )[ 1 ], R )
  boot.series.std	  <- matrix( NA, nrow = dim( data.boot )[ 1 ], R )
  
  #Remove empty columns
  data.boot <- cbind(data.boot, NA)
  data.boot<-data.boot[,-which(apply(data.boot,2,function(x)all(is.na(x))))]
  
  for( i in 1:R){
    n.series	<- dim( data.boot )[ 2 ]
    #Produces a series of random numbers used as an indication as to which series to select
    rand.num<- round( runif( n.series, min = 0.5, max = n.series + 0.49 ), digits = 0 )
    rand.series		<- data.boot[ , sort( rand.num )]
    site.chron.x 		<- site.chron( rand.series, aver.by.tree.input = aver.by.tree.input )
    aver.rand.series		<-site.chron.x$aver.site
    std.rand.series 		<- site.chron.x$std.site
    boot.series.mean[ , i]	<-aver.rand.series
    boot.series.std[ , i]	<-std.rand.series
    boot.stat[ i, ]		<-stat( rand.series )
  }
  
  org.stat		 <- stat( data.boot )
  bias		<- matrix( NA, 1, N.stat, dimnames = list( 1, names.stat) )
  std.error	<- matrix( NA, 1, N.stat, dimnames = list( 1, names.stat) )
  for( k in 1:N.stat){
    bias[ k ]<- 1 / R * sum( boot.stat[ , k] - org.stat[ k ] )
    std.error[ k ]	<- sqrt( 1 /(R-1 ) * sum(( boot.stat[ , k] - mean( boot.stat[ , k])) ^ 2))
  }
  results	 <- list( org.stat = org.stat, boot.stat = boot.stat, bias = bias,std.error = std.error, boot.series.mean = boot.series.mean,boot.series.std = boot.series.std)
  return( results )
}	


#'	Stat
#'  
#'  To be used internally in series.bootstrapped this is the summary that is calculated during bootstrapping. 
#'  @param data data matrix (produced within series.bootstrapped)
#'  @return results required output for series.bootstrapped
#'  @export
stat		<- function(data){
  m<-mean(data, na.rm = T)
  med<-median(data, na.rm = T)
  result<-c(m=m,med=med)
  return(result)
}
#' names.stat
#' 
#' @export
names.stat <- c('mean','median')
