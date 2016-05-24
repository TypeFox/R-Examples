#' Truncate series by range pithoffsets
#' 
#' The following function truncate the data by a given range from the estimated pith.		
#' @param	file.raw 	data file containing the raw ring widths, in mm
#' @param	file.stand	data file containing the standardized ring widths
#' @param	pithoffset	data set containing the pith offsets for each core, in mm.  
#' @param	range	The distance from the pith use for truncation, given in mm from the core. e.g. range <- c(1,200) truncates values outside this range.  
#' @return \item{sub.series.raw}{A truncated series of raw ring widths (in miromilmeter).}
#' \item{sub.series.stand}{A truncated series of standardized ring widths (in miromilmeter).} 
#' @examples
#' \dontrun{
#' data(ring.raw)
#' data(ring.stand)
#' data(dbh.po.nc)
#' #Subset near-pith is the material within 0 -20cm from the estimated pith
#' spline200.sub0.20.n   <- TruncSeriesPithoffset( ring.raw, ring.stand, dbh.po.nc, c(1,200))
#' # Subset far-pith is the material further than 20cm from the estimated pith
#' spline200.sub20.2000.n  <- TruncSeriesPithoffset( ring.raw, ring.stand, dbh.po.nc, c(200,200000))
#' # Whole dataset, through truncated functions to get in the same formate as the above two datasets
#' spline200.sub0.2000.n  <- TruncSeriesPithoffset( ring.raw, ring.stand, dbh.po.nc, c(00,200000))}

#' @export
TruncSeriesPithoffset	<- function( file.raw,file.stand, pithoffset, range) {
  
  
  dist.pith	<-	pithoffset$pithoffset
  
  # For each series, calculate the distance to the pith for each ring. 
  ts.cum.sum	<- NULL
  
  for( j in 1:dim(file.raw)[ 2 ]){
    n	<- 1
    while( (dimnames(file.raw)[[2]][ j ] != as.character(pithoffset$Series.ID)[ n ]) && 
             (n <= dim ( pithoffset )[ 1 ] ))		{ n<- n+1 }
    if( is.na( dist.pith[ n ]) || all(is.na(file.raw[ ,j] ) )){ 
      ts.cum.sum[[ j ]]	 <- ts( c(NA, NA), start = tsp ( file.raw )[ 1 ] )
    }else{
      series.na.omit		<- na.omit( file.raw[ , j]) 
      cum.sum.1		<- dist.pith[ n ] + series.na.omit[ 1 ]
      cum.sum 		<- c( cum.sum.1 )
      for(i in 1:length( series.na.omit )){
        cum.sum.2	<-cum.sum[ i ]+series.na.omit[ i+1 ]
        cum.sum	<-c( cum.sum, cum.sum.2)
      }
      ts.cum.sum[[ j ]]	<- ts( cum.sum, start = tsp (series.na.omit )[1])
    }
  }
  
  # Produce the truncated matrix
  sub.series.raw		<-NULL
  sub.series.stand		<-NULL
  sub.series.removed		<-NULL
  
  for(j in 1:dim(file.raw)[ 2 ]){
    if(	sum (ts.cum.sum[[ j ]], na.rm = T) == 0 ||
          ts.cum.sum[[ j ]][ 1 ] > range[2] ||
          max( ts.cum.sum[[ j ]], na.rm = T ) < range[1] ){
      sub.series.raw[[ j ]] 		<- ts(NA, start = tsp(file.raw)[1], end = tsp(file.raw)[2])
      sub.series.stand[[ j ]] 		<- ts(NA, start = tsp(file.raw)[1], end = tsp(file.raw)[2])
      sub.series.removed[[ j ]]	 <- file.raw[,j]
      
    }else{ 
      p <- 1
      while(
        ts.cum.sum[[j]][p] < range[1] && 
          p < length(ts.cum.sum[[j]])   ){ p		<-p+1}
      m	<-p
      while(
        ts.cum.sum[[j]][m] < range[2]  && m < length(ts.cum.sum[[j]])   &&  
          is.na( ts.cum.sum[[j]][m] ) != TRUE ) {m	<-m+1}
      if( m ==  length(ts.cum.sum[[j]])  ) { end.m <- tsp(ts.cum.sum[[j]])[2]
      }else{end.m <- tsp(ts.cum.sum[[j]])[1] + (m-1)}
      
      sub.series.raw[[ j ]]	 <- window(na.omit(file.raw[ , j]), start = tsp(na.omit(file.raw[ , j]))[1] + (p-1), 
                                       end = end.m )
      sub.series.stand[[ j ]]	 <- window(na.omit(file.stand[ , j]), start = tsp(na.omit(file.raw[ , j]))[1] + (p-1),
                                         end = end.m )
      
      if ( p == 1 ) {
        sub.series.removed[[ j ]] 	<- ts(NA, start = tsp(file.raw)[1],end = tsp(file.raw)[2])
      }else{
        sub.series.removed[[ j ]]	<- window( na.omit( file.raw[ , j]),
                                             start = tsp( na.omit( file.raw[ , j]))[ 1 ] , 
                                             end = tsp( na.omit( file.raw[ , j]))[ 1 ] + ( p-1 ) ) }
    }	
  }
  
  
  check.sum.raw		<- lapply( sub.series.raw, function( x ){sum( x )})
  
  function.list.to.matrix <- function( x, j){
    x.matrix<- NULL
    for ( k in 1: j ) {
      x.matrix 	<- cbind(x.matrix,x[[k]])  }
    x.matrix
  }
  
  sub.series.raw		<- function.list.to.matrix( sub.series.raw, length( sub.series.raw ) )	
  colnames(sub.series.raw)	<-dimnames(file.raw)[[2]]
  
  sub.series.stand	<- function.list.to.matrix( sub.series.stand, length( sub.series.stand ) )
  colnames(sub.series.stand)	<-dimnames(file.stand)[[2]]
  
  sub.series.removed	<- function.list.to.matrix( sub.series.removed, length( sub.series.removed ) )	
  colnames(sub.series.removed)	<-dimnames(file.raw)[[2]]
  
  #Output results as a list
  results	<- list (check.sum.raw = check.sum.raw, sub.series.removed = sub.series.removed, 
                   sub.series.raw=sub.series.raw, sub.series.stand=sub.series.stand, ts.cum.sum = ts.cum.sum)
  return( results )
}
