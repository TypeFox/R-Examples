leafsfirst.intpol.volu<-function(lst, dendat, rho, M=1000, grid=1)
{
itemnum<-length(lst$volume)
volume<-matrix(0,itemnum,1)
kapat<-matrix(0,itemnum,1)
d<-dim(dendat)[2]

volume.root<-montecarlo.ball(dendat,rho,M)
volume.sum<-itemnum*pi*rho^2
kappa<-volume.root/volume.sum

# 1) lasketaan kapat grid kpl

kapat.lyhyt<-matrix(0,grid,1)
levet.lyhyt<-matrix(0,grid,1)
kapat.lyhyt[1]<-kappa
levet.lyhyt[1]<-1
if (grid>1){
   levstep<-floor(itemnum/grid)
   or<-order(lst$level)
   for (i in 2:grid){
         levlok<-(i-1)*levstep
         note<-or[levlok]
         atomit<-lst$atomlist[note,1:lst$atomnumb[note]]
         pisteet<-matrix(dendat[atomit,],lst$atomnumb[note],d)
         volume.nyt<-montecarlo.ball(pisteet,rho,M)
         volume.sum<-lst$atomnumb[note]*pi*rho^2
         kapat.lyhyt[i]<-volume.nyt/volume.sum
         kapat.lyhyt[i]<-min(kapat.lyhyt[i],kapat.lyhyt[i-1])
         levet.lyhyt[i]<-levlok    
   }
}

# 2) interpoloidaan muut kapat

ra<-rank(lst$level)
for (i in 1:itemnum){
    ranko<-ra[i]
    lohko<-ceiling(grid*ranko/itemnum)
    kapa.ala<-kapat.lyhyt[lohko]
    if (lohko<grid) kapa.yla<-kapat.lyhyt[lohko+1] 
    else            kapa.yla<-kapat.lyhyt[grid]
    leve.ala<-levet.lyhyt[lohko]
    if (lohko<grid) leve.yla<-levet.lyhyt[lohko+1] 
    else            leve.yla<-itemnum
    kappa<-kapa.ala+(ranko-leve.ala)*(kapa.yla-kapa.ala)/(leve.yla-leve.ala)
    kapat[i]<-kappa
    volume.pot<-kappa*lst$atomnumb[i]*pi*rho^2
    #if (lst$parent[i]>0) volume[i]<-min(volume.pot,lst$volume[lst$parent[i]])
    volume[i]<-volume.pot
}

lst$volume<-volume
return(lst)
}


