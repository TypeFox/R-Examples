profgene<-function(values,recs,frekv=NULL,cvol=TRUE,ccen=TRUE,cfre=FALSE,
outlsets=TRUE,invalue=TRUE)
{

cu<-cumu(values,recs,frekv)
levels<-cu$levels
lsets<-cu$lsets
atoms<-cu$atoms
binfrek<-cu$frekv  #kullekin suorakaiteelle frekvenssi

alkublokki<-200
blokki<-50
links<-toucrec(atoms,alkublokki,blokki)

alkublokki2<-200
blokki2<-50
dentree<-decom(lsets,levels,links,alkublokki2,blokki2)
seplsets<-dentree$lsets
sepval<-dentree$levels
parents<-dentree$parents

if (cfre) nodefrek<-cfrekv(seplsets,binfrek) else nodefrek<-NULL 

if (ccen==TRUE) cvol<-TRUE
if (cvol){
  volum<-cvolum(seplsets,atoms)
  kerroin<-cinte(sepval,volum,parents) 
  sepvalnor<-sepval/kerroin
} 
else{  
  volum<-NULL
  sepvalnor<-NULL
}

if (ccen && cvol) centers<-ccente(seplsets,atoms,volum) else centers<-NULL

if (!(outlsets)) seplsets<-NULL
if (!(invalue)) sepval<-NULL

return(list(parent=parents,level=sepvalnor,invalue=sepval,
volume=volum,center=centers,nodefrek=nodefrek,lsets=seplsets))
#values: normeeratut arvot
#invalues: alkuperaiset frekvenssit/arvot 
#nodefrek: kunkin solmun frekvenssi
}

















