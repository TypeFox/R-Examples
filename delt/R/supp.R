supp<-function(dendat,epsi=0,blown=FALSE){
#Estimates the support of density
#
#dendat on n*xlkm matriisi
#epsi on tekn parametri
#
#Returns xlkm*2-matriisin
#
#kantajaksi estimoidaan [min-epsi,max+epsi]

n<-dim(dendat)[1]
xlkm<-length(dendat[1,])    #dendat matr sarakk lkm on muuttujien lkm
vast<-matrix(0,2*xlkm,1)  

for (i in 1:xlkm){

    minni<-min(dendat[,i])   
    maxxi<-max(dendat[,i])
    if (blown) epsi<-(maxxi-minni)/(2*(n-1))
   
    vast[2*i-1]<-minni-epsi     #sis valien alkupisteet
    vast[2*i]<-maxxi+epsi       #sis valien paatepisteet
}
return(vast)
}

