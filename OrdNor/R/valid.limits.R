valid.limits <-
function(plist, no.ord, no.norm){

validate.plist(plist, no.ord)


minmat = maxmat= diag(length(plist) + no.norm)

for (r in 2:nrow(minmat) ) {
for (c in 1:(r-1) ){

if(r != c) {

if (r<=length(plist) &  c<=length(plist) ) { 
minmax = Limit_forOO(plist[[r]], plist[[c]])
}
else if (r>length(plist) &  c>length(plist) ) { 
minmax = c(-1,1)
}
else if (r>length(plist) &  c<=length(plist) ){
minmax = Limit_forON(plist[[ c ]])
}
minmat[r,c] = minmax[1]
maxmat[r,c] = minmax[2] 
rm(minmax)
}
}
}

minmat= minmat + t(minmat) ; diag(minmat)=1
maxmat= maxmat + t(maxmat) ; diag(maxmat)=1
return (list(lower=minmat, upper=maxmat) )
}
