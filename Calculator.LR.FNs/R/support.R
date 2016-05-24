support <-
function(M, Left.fun = NULL, Right.fun = NULL){
range1 = M[1]-M[2]-M[3]-100
range2 = M[1]+M[2]+M[3]+100
x = seq(range1, range2, len = 200000)
if ( M[4] == 0 ) { y = Left.fun((M[1]-x)/M[2]) * (x<=M[1]) + Right.fun((x-M[1])/M[3]) * (M[1]<x) }
  else if ( M[4] == 1 ) { y = Right.fun((M[1]-x)/M[2]) * (x<=M[1]) + Left.fun((x-M[1])/M[3]) * (M[1]<x) }
  else if ( M[4] == 0.5 ) { y = Left.fun((M[1]-x)/M[2]) * (x<=M[1]) + Left.fun((x-M[1])/M[3]) * (M[1]<x) }
supp = c()
supp[1] = min(x[0<y & y<1])
supp[2] = max(x[0<y & y<1])
 if ( supp[1] == min(x) ) { supp[1] = -Inf }
 if ( supp[2] == max(x) ) { supp[2] = +Inf }
#print( noquote( paste("The support of fuzzy numver is interval:" ) ) )
return(supp )
  if (Left.fun == Right.fun+100 ) print(2) #Yek jomleye alaki choon CRAN majburet karde bud ke ...
}
