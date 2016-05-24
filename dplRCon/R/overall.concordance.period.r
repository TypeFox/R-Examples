#' Overall Concordance period
#' 
#' Produces the concordance indices for each year over the period of interest and the average of these indices. Trimmed and different levels of concordance indices in the overall concordance, size adjusted (mean concordance indices)  
#' @param x A matrix containing the standardized ring width indices for subset X, column heading tree ID,  row names is time
#' @param y As x, but for subset Y  	
#' @param x.boot the bootstrapped replicates for matrix x. with bootstrapped replicate as row names, and time as the column headings
#' @param y.boot the bootstrapped replicates for matrix y.
#' @param	min.series is a number indicating the minimum number of series for calculating the concordance statistic, recommend 10.
#' @param concordance.indices The function used for calculating the concordance indices for each t.
#' @param	period	Time period for which the concordance is to be calculated, vector of form c(start, end)	
#' @param trim.alpha this is the amount of trimming of the extreme concordance indices, default is 0 (no trimming). Recommend that a 0.005 be used.
#' @param concordance.beta This is used to calculated the proportion of indices above this cut off. Default is proportion above 0.5. 
#' @return \item{pre.in}{Concordance indices.}
#' \item{min.series}{Minimum number of series used to calculate concordance.}
#' \item{num.series.at.t.x}{Number of series in matrix X for each time, t.}  
#' \item{num.series.at.t.y}{Number of series in matrix Y for each time, t.}
#' \item{pre.in.greater.ts.period}{Concordance indices, as time series, when minimum number of series was available.}
#' \item{period}{The time period used to calculate the concordance.}
#' \item{pre.overall.period}{The Concordance Statistic.}
#' @examples
#' #loading data
#' \dontrun{
#' data(ring.raw)
#' data(ring.stand)
#' data(dbh.po.nc)
#' 
#' #Subset near-pith is the material within 0 -20cm from the estimated pith
#' spline200.sub0.20.n   <- TruncSeriesPithoffset( ring.raw, ring.stand, dbh.po.nc, c(1,200))
#' # Subset far-pith is the material further than 20cm from the estimated pith
#' spline200.sub20.2000.n  <- TruncSeriesPithoffset( ring.raw, ring.stand, dbh.po.nc, c(200,200000))
#' 
#' # series.bootstapped
#' boot.0.20   <-  series.bootstrap( spline200.sub0.20.n$sub.series.stand, 
#'    stat, 999, names.stat, aver.by.tree = FALSE)
#' boot.20.2000   <- series.bootstrap(spline200.sub20.2000.n$sub.series.stand, 
#'    stat, 999, names.stat, aver.by.tree = FALSE)
#'
#'overall.precision.HUP <- overall.concordance.period(spline200.sub20.2000.n$sub.series.stand , 
#'    spline200.sub0.20.n$sub.series.stand, 
#'    boot.20.2000$boot.series.mean,  boot.0.20$boot.series.mean ,
#'    1 , concordance.indices, c(1880,1999), trim.alpha=0.005, concordance.beta=0.5)
#'}
#' @export
overall.concordance.period<-function( x, y, x.boot, y.boot, min.series, concordance.indices, period , trim.alpha = 0, concordance.beta = 0.5) {
  # Additional mini functions
  # calculate vector of number of series in x and y for each time, t
  num.series.at.t.x	<- apply( x, 1, SumNotNa )
  num.series.at.t.y	<- apply( y, 1, SumNotNa )
  
  # Sorts the bootstrapped replicates for each time
  sort.boot.fun	<-	function( x ) {
    sort( x, na.last = TRUE ) }
  
  # sums the rows, or column, annoying 'na'
  SumNotNa	<-	function( x ) {
    sum( is.na( x ) == FALSE)
  }
  
  # When period is not provided use the start and end date of matrix X.
  if (missing( period ) ) {
    period 	<- c( tsp( x )[ 1 ], tsp( x )[ 2 ])
  }
  
  # sorting the bootstrapped replicates, sort by rows
  sort.boot.x		<- apply( x.boot, 1, sort.boot.fun )
  sort.boot.y		<- apply( y.boot, 1, sort.boot.fun )
  
  # Calculate concordance indices
  time		<- seq( from = 1, to = dim( x.boot )[1], 1)
  pre.in	<- matrix( NA, nrow = length( time ), ncol = 3, dimnames = list( time, c( "pre.x", "pre.y", "pre" ) ) )
  
  for( i in 1:length( time ) ) {
    pre.in[ i, ]	<- concordance.indices( sort.boot.x, sort.boot.y, time[ i ] , trim.alpha, nx = num.series.at.t.x, ny = num.series.at.t.y)
  }
  
  # Calculating concordance indices where the number of series is greater than min.series
  pre.in.greater	<- matrix( NA, nrow = length( time ), ncol = 3, dimnames = list( time, c( "pre.x", "pre.y", "pre" ) ) )
  
  greater		<- num.series.at.t.x > min.series & num.series.at.t.y > min.series
  for( i in 1:length( time ) ) {
    if (greater[ i ] == TRUE ) {
      pre.in.greater[ i, ] <- pre.in[ i, ]
    }}
  
  # Calculate the overall concordance for the period where there are more than min.series
  pre.in.greater.ts		<- ts( pre.in.greater, start = tsp( x )[ 1 ] )
  pre.in.greater.ts.period	<- window( pre.in.greater.ts, start = period[ 1 ], end = period[ 2 ] )
  greater.period		<- window( ts( greater, start = tsp( x )[ 1 ] ), start = period[ 1 ], end = period[ 2 ] )
  pre.overall.period		<- sum( pre.in.greater.ts.period[ , 3 ] > concordance.beta, na.rm = T ) / sum( greater.period, na.rm = T )
  
  # mean concordance indices
  pre.in.period <- window( ts(pre.in, start = tsp(x)[1] ), start = period[ 1 ], end = period[ 2 ] )
  num.series.at.t.x.period <- window( ts(num.series.at.t.x, start = tsp(x)[1] ), start = period[ 1 ], end = period[ 2 ] )
  num.series.at.t.y.period <- window( ts(num.series.at.t.y, start = tsp(x)[1] ), start = period[ 1 ], end = period[ 2 ] )
  mean.con <- mean(as.numeric(pre.in.period[,3]), na.rm = T)
  se.p  	<- (as.numeric(pre.in.period[,3])*(1 - as.numeric(pre.in.period[,3])))/( num.series.at.t.x.period* num.series.at.t.y.period)
  var.con <- (sum(se.p, na.rm = T))/sum(!is.na(as.numeric(pre.in.period[,3])))^2
  #se.p  	<- (as.numeric(pre.in.period[,3])*( as.numeric(pre.in.period[,3]))-1)
  #var.con <- sum(se.p, na.rm = T)/ sum(!is.na(as.numeric(pre.in.period[,3])))^2
  T.con 	<- (sum(as.numeric(pre.in.period[,3]), na.rm = T) - sum(!is.na(as.numeric(pre.in.period[,3])))) / sqrt(var.con)
  pvalue.con <- pt(T.con, df= mean(num.series.at.t.x.period* num.series.at.t.y.period, na.rm = T))
  ci.con	 <- c(mean.con - qt(0.975, df= mean(num.series.at.t.x.period* num.series.at.t.y.period, na.rm = T))*sqrt(var.con), mean.con + qt(0.975, df= mean(num.series.at.t.x.period* num.series.at.t.y.period, na.rm = T))*sqrt(var.con))
  ci.p <- NULL
  for(k in 1:length(pre.in.period[,3])){
    ci.p[[k]] 	<- exact.ci(as.numeric(pre.in.period[,3]), num.series.at.t.x.period, num.series.at.t.y.period, k) }
  #ci.p <- ts(ci.p, start = period[1], end=period[2])
  ci.p <- ts(ci.p, start = period[1])
  # Produce list of results to be returned
  results	<-	list( pre.in = pre.in, min.series = min.series, num.series.at.t.x = num.series.at.t.x, num.series.at.t.y = num.series.at.t.y, pre.in.greater.ts.period = pre.in.greater.ts.period, pre.in.period = pre.in.period , period = period, pre.overall.period = pre.overall.period, se.p=se.p, ci.p=ci.p,  mean.con = mean.con, var.con = var.con, T.con=T.con, pvalue.con = pvalue.con, ci.con = ci.con )
  return( results )
}
