sign <-
function(M){
supp = support(M)

if (supp[1] > 0) {return( noquote( paste0("Positive" ) ) )}
 else {if (supp[2] < 0) {return( noquote( paste0("Negative" ) ) )}
         else {return( noquote( paste0("non-positive and non negative" ) ) )} }
}
