roundpretty <- function(kvec,maxdig)
{
if (length(kvec)>1){
  result=rep(NA,length(kvec))
  for (i in 1:length(kvec)){
    result[i]=roundpretty.sub(kvec[i],maxdig)
  }
}
else result=roundpretty.sub(kvec,maxdig)
return(result)
} 

roundpretty.sub <- function(k,maxdig)
{
# round in a pretty way for output
# k ... number to be rounded pretty
# maxdig ... maximum number of digits after the coma
#

if (is.na(k)) {kr=k}
else {
if (k<0.00001) {kr=format(round(k,min(maxdig,8)),nsmall=min(maxdig,8),scientific=FALSE)}
else if (k<0.0001) {kr=format(round(k,min(maxdig,7)),nsmall=min(maxdig,7),scientific=FALSE)}
else if (k<0.001) {kr=format(round(k,min(maxdig,6)),nsmall=min(maxdig,6),scientific=FALSE)}
else if (k<0.01) {kr=format(round(k,min(maxdig,5)),nsmall=min(maxdig,5),scientific=FALSE)}
else if (k<0.1) {kr=format(round(k,min(maxdig,4)),nsmall=min(maxdig,4),scientific=FALSE)}
else if (k<1) {kr=format(round(k,min(maxdig,3)),nsmall=min(maxdig,3),scientific=FALSE)}
else if (k<10) {kr=format(round(k,min(maxdig,2)),nsmall=min(maxdig,2),scientific=FALSE)}
else if (k<100) {kr=format(round(k,min(maxdig,1)),nsmall=min(maxdig,1),scientific=FALSE)}
else {kr=format(round(k,min(maxdig,0)),nsmall=min(maxdig,0),scientific=FALSE)}
}
return(kr)
}
