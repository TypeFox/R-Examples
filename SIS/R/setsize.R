intensure <- function(i,l1,l2,k){
if(length(intersect(l1[1:i],l2[1:i])) >= k)
   return(i)
else
   return(intensure(i+1,l1,l2,k))
}

int.size.k <- function(l1,l2,k){
   iensure = intensure(k,l1=l1,l2=l2,k=k)
   ix01 = l1[1:iensure]
   ix02 = l2[1:iensure]
   ix0 = sort(intersect(ix01,ix02))
   return(ix0)
} 

calculate.nsis <- function(family,varISIS,n,p){
if(varISIS=="aggr") nsis = floor(n/log(n))
else{
  if(family=="gaussian"){nsis = floor(n/log(n))}
  if(family=="binomial"){nsis = floor(n/(4*log(n)))}
  if(family=="poisson"){nsis = floor(n/(2*log(n)))}
  if(family=="cox"){nsis = floor(n/(4*log(n)))}  
}
if(p < n) nsis = p
return(nsis)
}
