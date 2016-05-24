fs.calc<-function(complex,dendat,h,type="barys")
{
#type =  "barys", "means", "mins", "maxs"
lkm<-dim(complex)[1]
n<-dim(dendat)[1]
d<-dim(dendat)[2]
fs<-matrix(0,lkm,1)

if (type!="barys"){
   f<-matrix(0,n,1)
   for (i in 1:n){
      arg<-dendat[i,]
      f[i]<-kernesti.dens(arg,dendat,h)
   }
   for (i in 1:lkm){
       vs<-complex[i,]
       vals<-f[vs]
       if (type=="means") fs[i]<-mean(vals) 
       if (type=="maxs") fs[i]<-max(vals) 
       if (type=="mins") fs[i]<-min(vals) 
   }
}

if (type=="barys"){
   #barys<-matrix(0,lkm,d)
   for (i in 1:lkm){
      simple<-complex[i,]
      simp<-dendat[simple,]
      arg<-colSums(simp)/(d+1)  #arg<-barys[i,]
      fs[i]<-kernesti.dens(arg,dendat,h)
   }
}

return(fs)
}

