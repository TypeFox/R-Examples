modetestgauss<-function(lst,n)
{

len<-length(lst$parent)
testvec<-matrix(0,len,1)    #this is output

em<-excmas(lst)

for (i in 1:len){

   if (lst$parent[i]!=0) val<-lst$level[lst$parent[i]]
   else val<-0

   a<-sqrt(n)*em[i]/sqrt(val*lst$volume[i])
   testvec[i]<-2*(1-pnorm(a))
}

return(testvec)
}
