profkernCRC<-function(dendat,h,N,Q,cvol=TRUE,ccen=TRUE,#cfre=FALSE,
kernel="epane",compoinfo=FALSE,trunc=3,threshold=0.0000001,katka=NULL,hw=NULL)
{
#dyn.load("/home/jsk/kerle/kerleCversio")
#pk2<-profkernCRC(dendat,h,N,Q)
#
#set.seed(seed=1)
#dendat<-matrix(rnorm(20),10)
#h<-1 
#N<-c(8,8)
#Q<-3
#
n<-dim(dendat)[1]
d<-length(N)

if (is.null(hw)) weig<-rep(1/n,n) 
else{
   weig<-weightsit(n,hw)

   dendatnew<-dendat
   weignew<-weig
   cumul<-0
   for (i in 1:n){
        if (weig[i]>0){
            cumul<-cumul+1
            dendatnew[cumul,]<-dendat[i,]
            weignew[cumul]<-weig[i] 
        }
   }
   dendat<-dendatnew[1:cumul,]
   weig<-weignew[1:cumul]
   n<-cumul
}
inweig<-matrix(0,n+1,1)
inweig[2:(n+1)]<-weig

hnum<-length(h)
mnn<-maxnodenum(dendat,h,N,n,d)
extMaxnode<-mnn$maxnode
extMaxvals<-mnn$maxpositive
#
if (hnum>1){
 inh<-matrix(0,hnum+1,1)
 inh[2:(hnum+1)]<-h
}
else{
 inh<-h
}
inN<-matrix(0,d+1,1)
inN[2:(d+1)]<-N

if (kernel=="epane") kertype<-1
else kertype<-2  # gaussian

kg<-.C("kergrid",
               as.integer(extMaxnode),
               as.integer(extMaxvals),
               as.double(dendat),
               as.double(inh),
               as.integer(inN),
               as.integer(n),
               as.integer(hnum),
               as.integer(d),
               as.integer(kertype),
               as.double(trunc),
               as.double(threshold),  
               as.double(inweig),        
               ioleft = integer(extMaxnode+1),
               ioright = integer(extMaxnode+1),
               ioparent = integer(extMaxnode+1),
               infopointer = integer(extMaxnode+1),
               iolow = integer(extMaxnode+1),
               ioupp = integer(extMaxnode+1),
               value = double(hnum*extMaxvals),
               index = integer(d*extMaxvals),
               nodefinder = integer(extMaxvals),
               numpositive = integer(1),
               numnode = integer(1),
PACKAGE="denpro")

left<-kg$ioleft[2:(kg$numnode+1)]
right<-kg$ioright[2:(kg$numnode+1)]
parent<-kg$ioparent[2:(kg$numnode+1)]
infopointer<-kg$infopointer[2:(kg$numnode+1)]
iolow<-kg$iolow[2:(kg$numnode+1)]
ioupp<-kg$ioupp[2:(kg$numnode+1)]

value<-kg$value[2:(kg$numpositive+1)]
nodefinder<-kg$nodefinder[2:(kg$numpositive+1)]
vecindex<-kg$index[2:(d*kg$numpositive+1)]
index<-matrix(0,kg$numpositive,d)
for (i in 1:kg$numpositive){
  for (j in 1:d){
     index[i,j]<-vecindex[(i-1)*d+j]
  }
}

nodenumOfDyaker<-length(left)

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
levnum<-length(levseq)
    
inlevseq<-matrix(0,length(levseq)+1,1)
inlevseq[2:(length(levseq)+1)]<-levseq
inN<-matrix(0,d+1,1)
inN[2:(d+1)]<-N
inleft<-matrix(0,length(left)+1,1)
inleft[2:(length(left)+1)]<-left
inright<-matrix(0,length(left)+1,1)
inright[2:(length(left)+1)]<-right
inparent<-matrix(0,length(left)+1,1)
inparent[2:(length(left)+1)]<-parent
invalue<-matrix(0,length(value)+1,1)
invalue[2:(length(value)+1)]<-value
#inindex<-matrix(0,dim(kg$index)[1]+1,dim(kg$index)[2]+1)
#for (i in 1:dim(kg$index)[1]){
#  inindex[i+1,]<-c(0,kg$index[i,])
#}
innodefinder<-matrix(0,length(nodefinder)+1,1)
innodefinder[2:(length(nodefinder)+1)]<-nodefinder

# Tama koodi on jo kergrid:ssa, lasketaan volume of one atom
minim<-matrix(0,d,1)  #minim is d-vector of minimums
maxim<-matrix(0,d,1)
for (i in 1:d){
  minim[i]<-min(dendat[,i])
  maxim[i]<-max(dendat[,i])
}
delta<-(maxim-minim+2*h)/(N+1)  
volofatom<-prod(delta)

inminim<-matrix(0,d+1,1)
inminim[2:(d+1)]<-minim
indelta<-matrix(0,d+1,1)
indelta[2:(d+1)]<-delta

if (!is.null(katka)){
   invalue2<-invalue
   lenni<-length(invalue)
   for (i in 1:lenni){
      if (invalue[i]>=katka) invalue2[i]<-katka
   }
   invalue<-invalue2
}

dentree<-.C("decomdyaC",
               as.integer(numofall),
               as.integer(atomnum),
               as.double(inlevseq),
               as.integer(inN),
               as.integer(d),         
               as.integer(levnum),   
               as.double(volofatom),
               as.double(inminim),
               as.double(h),
               as.double(indelta),
               as.integer(nodenumOfDyaker),
               as.integer(inleft),
               as.integer(inright),
               as.integer(inparent), 
               as.double(invalue),
               as.integer(index),
               as.integer(innodefinder),
               level = double(numofall+1),
               parent = integer(numofall+1),
               component = integer(numofall+1),
               volume = double(numofall+1),
               center = double(d*numofall+1),
               efek = integer(1),
               AtomlistAtom = integer(numofall+1),
               AtomlistNext = integer(numofall+1),
               numOfAtoms = integer(1),
PACKAGE="denpro")

AtomlistAtom<-dentree$AtomlistAtom[2:(dentree$numOfAtoms+1)]
AtomlistNext<-dentree$AtomlistNext[2:(dentree$numOfAtoms+1)]

invalue<-dentree$level[2:(dentree$efek+1)]
parent<-dentree$parent[2:(dentree$efek+1)]
volume<-dentree$volume[2:(dentree$efek+1)]
component<-dentree$component[2:(dentree$efek+1)]
kerroin<-cinte(invalue,volume,parent) 
sepvalnor<-invalue/kerroin
veccenter<-dentree$center[2:(d*dentree$efek+1)]
center<-matrix(0,dentree$efek,d)
for (i in 1:dentree$efek){
  for (j in 1:d){
     center[i,j]<-veccenter[(i-1)*d+j]
  }
}
center<-t(center)

#if (cfre) nodefrek<-cfrekvdya(seplsets,binfrek) else nodefrek<-NULL 

#clus<-F
#if (clus){
#   clustervecs<-cluskern(compo,component,AtomlistAtom,AtomlistNext,kg,dendat,
#   h,N)
#}
#else{
#   clustervecs<-NULL
#}

if (compoinfo)

  return(list(parent=parent,level=sepvalnor,invalue=invalue,
  volume=volume,center=center,
  component=component,
  AtomlistAtom=AtomlistAtom,AtomlistNext=AtomlistNext,index=index))

else

  return(list(parent=parent,level=sepvalnor,invalue=invalue,
  volume=volume,center=center,n=n))

#,nodefrek=nodefrek,clustervec=clustervecs))
#
#values: normeeratut arvot
#invalues: normeeraamattomat arvot 
#nodefrek: kunkin solmun frekvenssi

}












