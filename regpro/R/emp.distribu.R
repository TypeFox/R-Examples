emp.distribu<-function(arg,dendat)
{
d<-length(arg)
n<-dim(dendat)[1]

if (d>1){
   val<-matrix(0,d,1)
   for (kk in 1:d){
        cateca<-c(dendat[,kk],arg[kk])
        or<-rank(cateca)
        val[kk]<-or[n+1]/(n+2)
   }
}
else{
        cateca<-c(dendat,arg)
        or<-rank(cateca)
        val<-or[n+1]/(n+2) 
}

return(val)
}


