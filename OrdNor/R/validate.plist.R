validate.plist <-
function(plist, no.ord){
warn=FALSE
if (length(plist) != no.ord ) {warning("The value of no.ord and the number of probability vectors in the list do not match."); warn=TRUE} 
for (i in 1:length(plist) ) {
p=plist[[i]]
if ( max(p)>=1) {warning("A probability value cannot be greater than 1. The probability vector ", i, " is not valid."); warn=TRUE }
if ( !(min(p) >0) ) {warning("A probability value cannot be less than 0. The probability vector ", i, " is not valid."); warn=TRUE }

#if ( max(table(p) )!=1 ) {stop("The probability vector ", "i", is not in the form of cumulative probabilities")}
if ( sum( abs( p[order(p)]-p) ) > 0 ) {warning("The probability vector ", i, " is not in the form of cumulative probabilities") ; warn=TRUE }
}
if(warn) {stop("plist is not valid")}

}
