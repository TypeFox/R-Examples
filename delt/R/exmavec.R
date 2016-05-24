exmavec<-function(vol,obsnum,n,lambda)
{

len<-length(vol)
resu<-matrix(0,len,1)

for (i in 1:len){
  
  resu[i]<-exma(vol[i],obsnum[i],n,lambda)

}

return(resu)
}



