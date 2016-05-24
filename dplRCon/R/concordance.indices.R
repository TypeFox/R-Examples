#' Concordance 
#' 
#' This is an internal function called in the function "Overall Concordance". It is used to calculate the concordance indices for each time, t, using the bootstrapped means. 
#' @param sort.boot.x A matrix of the bootstrapped means from subset X
#' @param sort.boot.y A matrix of the bootstrapped means from subset Y
#' @param t The time index
#' @param trim.alpha The level of outside trimming
#' @param nx A vector providing the number of series at each time for subset X
#' @param ny A vector providing the number of series at each time for subset Y
#' @details The procedure for calculating the concordance is provided in "Concordance: A measure of similarity between matrices of time series with applications in dendroclimatology". 
#' @return A vector containing the a, b and the concordance 
#' @export
concordance.indices   <-	function( sort.boot.x, sort.boot.y, t , trim.alpha, nx, ny) { 
  
  #setup 
  sort.boot	<- cbind( sort.boot.x[ , t ], sort.boot.y[ ,t ] )
  
  SumNotNa	<-	function( x ) {
    sum( is.na( x ) == FALSE)
  }
  
  Xr		<- SumNotNa( sort.boot[ , 1] )
  Yr		<- SumNotNa( sort.boot[ , 2] )
  x.over.y <- NULL
  y.over.x <- NULL
  
  #trims the outside means, based on alpha
  number.trim.x <- round(Xr*trim.alpha,0)
  number.trim.y <- round(Yr*trim.alpha,0)
  Xr <- Xr - number.trim.x
  Yr <- Yr - number.trim.y
  X1 <- 1+ number.trim.x
  Y1 <- 1+ number.trim.y
  
  #size adjusts the bootstrapped means
  sort.boot[,1] = sqrt(nx[t])* (sort.boot[,1]-mean(sort.boot[,1], na.rm = T)) + mean(sort.boot[,1], na.rm = T)
  sort.boot[,2] = sqrt(ny[t])* (sort.boot[,2]-mean(sort.boot[,2], na.rm = T)) + mean(sort.boot[,2], na.rm = T)
  
  # Counts the number of overlapping bootstrapped replicates for each time, t
  if (Xr == 0 | Yr == 0  ){
    pre.x	 <-NA
    pre.y 	<- NA
    pre	<- NA
  }else{
    if (sort.boot[X1 , 1] >= sort.boot[ Y1, 2] & 
          sort.boot[ Xr, 1] <= sort.boot[ Yr, 2]) {
      x.over.y	<- Xr-X1+1
      Yxr		<- sum( sort.boot[ Xr, 1] > sort.boot[ Y1:Yr, 2], na.rm = T )
      Yx1		<- sum( sort.boot[ X1, 1]    > sort.boot[ Y1:Yr, 2], na.rm = T )
      y.over.x	<- Yxr - Yx1
    }
    if (sort.boot[ X1, 1] >= sort.boot[ Y1, 2] & 
          sort.boot[ Xr, 1] >= sort.boot[ Yr, 2] ) {
      Xyr 		<- sum( sort.boot[X1:Xr , 1] < sort.boot[ Yr, 2], na.rm = T )
      x.over.y	<- Xyr
      y.over.x	<- Yr -Y1+1- sum( sort.boot[ X1, 1] > sort.boot[ Y1:Yr, 2], na.rm = T)
    }
    if ( sort.boot[ X1, 1] <= sort.boot[ Y1, 2] &
           sort.boot[ Xr, 1] >= sort.boot[ Yr, 2] ){
      Xy1		<- sum( sort.boot[ X1:Xr, 1] < sort.boot[ Y1, 2], na.rm = T )
      Xyr		<- sum( sort.boot[X1:Xr , 1] > sort.boot[ Yr, 2], na.rm = T )
      x.over.y	<- Xr -X1+1- ( Xy1 + Xyr )
      y.over.x	<- Yr- Y1+1
    }
    if ( sort.boot[ X1, 1] <= sort.boot[ Y1, 2] &
           sort.boot[ Xr, 1] <= sort.boot[ Yr, 2] ) {
      Xy1 		<- sum( sort.boot[X1:Xr , 1] < sort.boot[ Y1, 2], na.rm = T )
      x.over.y	<- Xr -X1+1- Xy1
      y.over.x	<- sum( sort.boot[ Xr, 1] > sort.boot[Y1:Yr , 2], na.rm = T )
    }
    if (sort.boot[ X1, 1] == sort.boot[ Y1, 2] &
          sort.boot[ Xr, 1] == sort.boot[ Yr, 2] ) {
      x.over.y 	<- Xr-X1+1
      y.over.x 	<- Yr-X1+1
    }
    # Summaries the above proportion of overlapping bootstrapped replicates	
    pre.x	<- x.over.y / (Xr- X1+1)
    pre.y 	<- y.over.x / (Yr- Y1 +1)
    pre	<- pre.x * pre.y
  }
  # produces output matrix
  results	<- cbind( pre.x, pre.y, pre )
  return( results )
}
