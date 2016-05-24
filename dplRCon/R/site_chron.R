#' site.chron
#' 
#'This function calculate the site chronology for the data matrix 
#' @param	data.m	A data matrix contains standardized ring width indices 	
#' @param	aver.by.tree.input If True then averages tree series, then average tree means. If False then average all series 
#' @return	\item{aver.site}{The average for the site.}
#' \item{std.site}{The standard deviation of the site.}
#' \item{aver.tree}{A time series for average of each tree.} 
#' \item{number.trees}{The total number of trees used to produce the site chronology.}
#' \item{core.per.tree}{The number of cores for each tree.}
#' @examples
#' \dontrun{
#' site.full <- site.chron(spline200.sub0.2000.n$sub.series.stand, aver.by.tree=FALSE)
#' }
#' @export
site.chron	<-	function( data.m, aver.by.tree.input ) { 
  names		<- dimnames( data.m )[[ 2 ]]
  names.sub	<- substring( names, 1, 5 )
  
  aver.tree	<- array( NA, dim = dim( data.m), dimnames = dimnames( data.m ) )
  std.tree		<- array( NA, dim = dim( data.m), dimnames = dimnames( data.m ) )
  
  core.per.tree	<- NULL
  
  i	<- 1
  j	<- 1
  while( i < length( names.sub ) ){ 
    j	<- i
    while( names.sub[ i ] == names.sub[ i+1 ] && i < length(names.sub)){
      i = i+1
    }
    tree		<- data.m[ , j : i ]	#extracts cores from a particular tree
    core.per.tree	<-c(core.per.tree, dim( tree )[ 2 ] )		
    
    if( i==j ) {
      aver.tree[ , j]	<- tree
      std.tree[ , j]	<- tree
    }else{
      aver.tree[ , j]	<- apply( tree, 1, function( x )mean( x, na.rm = T))
      #std.tree[ , j]	<- apply( tree, 1,
      #	function( x ){sqrt( var( x, use = "pairwise.complete.obs"))})
      #std.tree[ , j]	<- apply( tree, 1,
      #	function( x ){sqrt( var( na.omit(x)))})
      
    }
    i	<- i + 1
  }	
  if( aver.by.tree.input == TRUE ) { 
    # Performs this step it require cores to be average by tree, before averaging trees
    aver.site 	<- apply( aver.tree, 1, function( x )mean( x, na.rm = T))
    #std.site 		<- apply( std.tree, 1,
    #	function( x ){sqrt( var( na.omit(x)))})
  }else{
    aver.site	<- apply(data.m, 1, function( x )mean( x, na.rm = T ))
    std.site 		<- apply( data.m, 1, 
                         function( x ){sqrt( var( na.omit(x) ))})
  }	
  
  number.trees 	<- sum( core.per.tree )
  aver.site 		<- ts(aver.site, start = tsp(data.m)[1] )
  std.site 		<- ts(std.site, start = tsp(data.m)[1] )
  results 		<- list( aver.site=aver.site, std.site = std.site,
                     aver.tree = aver.tree, number.trees = number.trees,
                     core.per.tree = core.per.tree)
  return(results)
}
