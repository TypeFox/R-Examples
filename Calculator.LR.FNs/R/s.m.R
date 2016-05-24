s.m <-
function(k, N)  {

if ( messages(N) != 1 )  { return( messages(N) ) }
if ( messages(k) != 1 )  { return( messages(k) ) }

if ( length(k)==4 & length(N)==1) {
   zarf = N
   N[1] = k[1]
   N[2] = k[2]
   N[3] = k[3]
   N[4] = k[4]
   k = zarf
   }

 if (k==0) { return( noquote( paste0(" The scalar multiplication is not defined for zero " ) ) ) }
  else {
   a1 = k*N[1]
   a2 = k* (N[2]*(k>0)-N[3]*(k<0))
   a3 = k* (N[3]*(k>0)-N[2]*(k<0))
   a4 = N[4]
   print( noquote( paste0("the result of scalar multiplication is  (core = ", a1, ", left spread = " , a2, ", right spread = " , a3, ")"
     ,  if ( a4 == 0 ) { paste0(" LR" ) }  else if ( a4 == 1 ) { paste0(" RL" ) }  else { paste0(" L" ) }  ) ) )

   return( invisible( c(a1,a2,a3,a4) ) )
 }
}
