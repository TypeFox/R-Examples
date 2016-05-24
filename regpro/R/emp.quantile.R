emp.quantile<-function(arg,dendat)
{
d<-length(arg)

if (d>1){
   val<-matrix(0,d,1)
   n<-dim(dendat)[1] 
   for (kk in 1:d){
        or<-rank(dendat[,kk])
        uni<-or/(n+1)
        cateca<-c(uni,arg[kk])
        or2<-rank(cateca)
        inter<-or2[length(cateca)]
        if (inter==(n+1)) inter<-n
        so<-sort(dendat[,kk])
        val[kk]<-so[inter]

   }   
}
else{
        n<-length(dendat) 
        or<-rank(dendat)
        uni<-or/(n+1)
        cateca<-c(uni,arg)
        or2<-rank(cateca)
        inter<-or2[length(cateca)]
        if (inter==(n+1)) inter<-n
        val<-sort(dendat)[inter]
}

return(val)
}





