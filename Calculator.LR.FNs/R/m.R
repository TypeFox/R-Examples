m <-
function(M, N)  {
options(warn = -1) 

if ( messages(M) != 1 )  { return( messages(M) ) }
if ( messages(N) != 1 )  { return( messages(N) ) }


if ( M[4] != N[4]) 
  {
   return( noquote( paste0("Production has NOT a closed form of a LR fuzzy number" ) ) )
  } 
 else if (  ( sign(M) == "Positive" )  &  ( sign(N) == "Positive" )  ) 
   {
    a1 = M[1]*N[1]
    a2 = ( M[1]*N[2] ) + ( N[1]*M[2] )
    a3 = ( M[1]*N[3] ) + ( N[1]*M[3] )
    a4 = (M[4]+N[4])/2
    print( noquote( paste0("the result of multiplication is approximately  (core = ", a1, ", left spread = " , a2, ", right spread = " , a3, ")"
     ,  if ( a4 == 0 ) { paste0(" LR" ) }  else if ( a4 == 1 ) { paste0(" RL" ) }  else { paste0(" L" ) }  ) ) )
    return( invisible( c(a1,a2,a3,a4) ) )
   }
 else if (  ( sign(M) == "Negative" )  &  ( sign(N) == "Negative" )  ) 
   {
    a1 = M[1]*N[1]
    a2 = -( M[1]*N[2] ) - ( N[1]*M[2] )
    a3 = -( M[1]*N[3] ) - ( N[1]*M[3] )
    a4 = abs( M[4]-1 )
    print( noquote( paste0("the result of multiplication is approximately  (core = ", a1, ", left spread = " , a2, ", right spread = " , a3, ")"
     ,  if ( a4 == 0 ) { paste0(" LR" ) }  else if ( a4 == 1 ) { paste0(" RL" ) }  else { paste0(" L" ) }  ) ) )
    return( invisible( c(a1,a2,a3,a4) ) )
   }
 else if (  ( sign(M) == "Positive" )  &  ( sign(N) == "Negative" )  ) 
   {
    a1 = M[1]*N[1]
    a2 = ( M[1]*N[2] ) - ( N[1]*M[3] )
    a3 = ( M[1]*N[3] ) - ( N[1]*M[2] )
    a4 = abs( M[4]-1 )
    print( noquote( paste0("the result of multiplication is approximately  (core = ", a1, ", left spread = " , a2, ", right spread = " , a3, ")"
     ,  if ( a4 == 0 ) { paste0(" LR" ) }  else if ( a4 == 1 ) { paste0(" RL" ) }  else { paste0(" L" ) }  ) ) )
   return( invisible( c(a1,a2,a3,a4) ) )
   }
 else if (  ( sign(M) == "Negative" )  &  ( sign(N) == "Positive" )  ) 
   {
    a1 = M[1]*N[1]
    a2 = ( N[1]*M[2] ) - ( M[1]*N[3] )
    a3 = ( N[1]*M[3] ) - ( M[1]*N[2] )
    a4 = abs( N[4]-1 )
    print( noquote( paste0("the result of multiplication is approximately  (core = ", a1, ", left spread = " , a2, ", right spread = " , a3, ")"
     ,  if ( a4 == 0 ) { paste0(" LR" ) }  else if ( a4 == 1 ) { paste0(" RL" ) }  else { paste0(" L" ) }  ) ) )
    return( invisible( c(a1,a2,a3,a4) ) )
   }
 else
  {
   return( noquote( paste0("A regular approximation is not defined for multiplication since at least one of LR fuzzy numbers is non-positive and non-negative fuzzy number")))
  } 

}
