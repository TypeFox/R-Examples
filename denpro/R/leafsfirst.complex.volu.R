leafsfirst.complex.volu<-function(lst,dendat,complex,rho,vols,M=1000,grid=1,
seed=1)
{
itemnum<-length(lst$volume)
volume<-matrix(0,itemnum,1)
kapat<-matrix(0,itemnum,1)
d<-dim(dendat)[2]

if (grid==0){

for (note in 1:itemnum){
  atomit<-lst$atomlist[note,1:lst$atomnumb[note]]
  pisteet<-matrix(complex[atomit,],lst$atomnumb[note],d+1)
  voltti<-montecarlo.complex(dendat,pisteet,rho,M,seed=seed)
  #if (lst$parent[note]>0) voltti<-min(voltti,lst$volume[lst$parent[note]])
  volume[note]<-voltti
}
lst$volume<-volume

}
else{

volume.root<-montecarlo.complex(dendat,complex,rho,M,seed=seed)
volume.sum<-sum(vols)  #itemnum*pi*rho^2
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
         pisteet<-matrix(complex[atomit,],lst$atomnumb[note],d+1)
         volume.nyt<-montecarlo.complex(dendat,pisteet,rho,M,seed=seed)
         volume.sum<-sum(vols[atomit])    #lst$atomnumb[note]*rho^2/2
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
    atomit<-lst$atomlist[i,1:lst$atomnumb[i]]
    volume.pot<-kappa*sum(vols[atomit]) 
    #volume.pot<-kappa*lst$atomnumb[i]*rho^2/2
    #if (lst$parent[i]>0) volume[i]<-min(volume.pot,lst$volume[lst$parent[i]])
    volume[i]<-volume.pot
}
lst$volume<-volume

}

lst$volume<-volume
return(lst)
}

