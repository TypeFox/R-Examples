ordinalize <-
function(pvec, z) {

validate.plist(list(pvec),1)
Xord = numeric(length(z))
for (r in 1:length(pvec) ){
if (r !=length(pvec) ) {
t1 = qnorm(pvec[r])
t2 = qnorm(pvec[r+1] )
Xord[(t1<z)&(z<=t2)]= r
} else  {
Xord[z>qnorm(pvec[r])]= r
}
}
return(Xord+1)
}
