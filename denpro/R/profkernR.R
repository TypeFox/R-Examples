profkernR<-function(kg,dendat,h,N,Q,frekv=NULL,cvol=TRUE,ccen=TRUE,cfre=FALSE){

#set.seed(seed=1)
#dendat<-matrix(rnorm(20),10)
#h<-1 
#N<-c(4,4)
#Q<-3

#kg<-kergrid(dendat,h,N)

nodenumOfDyaker<-length(kg$left)

value<-kg$value
maxval<-max(value)
minval<-min(value)
step<-(maxval-minval)/Q
levseq<-seq(from=minval,to=maxval-step,by=step)

levfrekv<-matrix(0,Q,1)
atomnum<-length(value)
for (i in 1:atomnum){
   for (j in 1:Q){
       if (value[i]>=levseq[j]){
          levfrekv[j]<-levfrekv[j]+1
       }
   }
}
numofall<-sum(levfrekv)
    
dentree<-decomdya(numofall,atomnum,levseq,kg,N,nodenumOfDyaker)
invalue<-dentree$level
parent<-dentree$parent
component<-dentree$component
AtomlistAtom<-dentree$AtomlistAtom
AtomlistNext<-dentree$AtomlistNext

# Tama koodi on jo kergrid:ssa, lasketaan volume of one atom
d<-length(N)
minim<-matrix(0,d,1)  #minim is d-vector of minimums
maxim<-matrix(0,d,1)
for (i in 1:d){
  minim[i]<-min(dendat[,i])
  maxim[i]<-max(dendat[,i])
}
delta<-(maxim-minim+2*h)/(N+1)  
volofatom<-prod(delta)

#if (cfre) nodefrek<-cfrekvdya(seplsets,binfrek) else nodefrek<-NULL 

if (ccen==TRUE) cvol<-TRUE
if (cvol){
  volume<-cvolumdya(volofatom,component,AtomlistNext)
  kerroin<-cinte(invalue,volume,parent) 
  sepvalnor<-invalue/kerroin
} 
else{  
  volume<-NULL
  sepvalnor<-NULL
}

if (ccen && cvol){
  index<-kg$index
  d<-dim(dendat)[2]
  center<-ccentedya(volofatom,component,AtomlistNext,AtomlistAtom,
                    volume,minim,h,delta,index,d)
  }
  else{
      center<-NULL
  }

return(list(parent=parent,level=sepvalnor,invalue=invalue,
volume=volume,center=center))#,nodefrek=nodefrek))
#values: normeeratut arvot
#invalues: alkuperaiset frekvenssit/arvot 
#nodefrek: kunkin solmun frekvenssi
}








