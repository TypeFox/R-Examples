a <-
function(M, N)  {
options(warn = -1) 

if ( messages(M) != 1 )  { return( messages(M) ) }
if ( messages(N) != 1 )  { return( messages(N) ) }

if ( M[4] != N[4] ) 
  {
   return( noquote( paste0("Addition has NOT a closed form of a LR fuzzy number" ) ) )
  } 
 else
   {
    a1 = M[1]+N[1]
    a2 = M[2]+N[2]
    a3 = M[3]+N[3]
    a4 = (M[4]+N[4])/2
    print( noquote( paste0("the result of addition is  (core = ", a1, ", left spread = " , a2, ", right spread = " , a3, ")"
     ,  if ( a4 == 0 ) { paste0(" LR" ) }  else if ( a4 == 1 ) { paste0(" RL" ) }  else { paste0(" L" ) }  ) ) )

    return( invisible( c(a1,a2,a3,a4) ) )
   }
}
