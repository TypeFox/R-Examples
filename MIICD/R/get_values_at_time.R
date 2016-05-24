

get_values_at_times <- function( j , values , times , at , list = F , unique = T ){
  if(list==F){
    
  sapply( at , function( x ) {
    values[ , j ][ tail( which( sort( unique( times[ , j ] ) )  <= x ) , 1 ) ]
  }
  )
  }else{
    if(unique){      
      sapply( at , function( x ) {
    values[[j]][ tail( which( sort( ( times[[j]] ) )  <= x ) , 1 ) ] 
  }
  )
}else{
      
sapply( at , function( x ) {
    values[[j]][ tail( which( sort( unique( times[[j]] ) )  <= x ) , 1 ) ] 
  }
  )
}
  }
}
