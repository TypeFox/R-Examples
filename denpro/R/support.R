support<-function(dendat,epsi=0)
{
#estimoi kantajan tih datan perusteella
#dendat on n*xlkm matriisi
#epsi on tekn parametri
#kantajaksi estimoidaan [min-epsi,max+epsi]
#palauttaa xlkm*2-matriisin

xlkm<-length(dendat[1,])    #dendat matr sarakk lkm on muuttujien lkm
vast<-matrix(0,xlkm,2)  
i<-1
while (i<=xlkm){
    vast[i,1]<-min(dendat[,i])-epsi     #sis valien alkupisteet
    vast[i,2]<-max(dendat[,i])+epsi     #sis valien paatepisteet
    i<-i+1
}
return(vast)
}


