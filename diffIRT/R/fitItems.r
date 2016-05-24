fitItems=function(nit,order){
  tmp=matrix(0,1,nit)
  c=list()
  for(i in 1:order) c[[i]]=combn(nit,i)
  index=lapply(c,ncol)
  for(i in 1:order)
  for(j in 1:index[[i]]){{
    aa=matrix(0,1,nit)
    aa[1,c[[i]][,j]]=1
  tmp=rbind(tmp,aa)
  }}
  return(tmp)
}

