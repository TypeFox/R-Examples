
index.prediction=function(res, x){
k=nrow(res)
scor=rep(0,nrow(x))
for(i in 1:k){
 jmax=res[i,"jmax"]
 cutp=res[i,"cutp"]
 maxdir=res[i,"maxdir"]
  if(maxdir==1){ temp=1*(x[,jmax]> cutp)}
  if(maxdir==-1){ temp=1*(x[,jmax]< cutp)}
 scor=scor+temp
}
return(scor)
}


