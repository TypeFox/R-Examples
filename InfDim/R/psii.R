psii <-
function(adjusted_age){
n=length(adjusted_age)

adjust_temp=rep(adjusted_age,each=n)
j_temp=rep(0:(n-1),n)

psii_elements=rep(NA,n*n)
for(i in 1:(n*n)){
psii_elements[i]=p.function(j=j_temp[i],x=adjust_temp[i])
}
psii_temp=matrix(psii_elements,byrow=TRUE,nrow=n,ncol=n)
return(psii_temp)
}

