proftree<-function(tr,
Q=NULL,frekv=NULL,cvol=TRUE,ccen=TRUE,cfre=FALSE)
{

d<-dim(tr$upp)[2]

if (tr$left[1]==0){
  parent=c(0)
  sepvalnor=c(tr$mean[1])
  invalue=c(tr$mean[1])
  volume=c(tr$volume[1])
  rec<-matrix(0,2*d,1)
  for (j in 1:d){
     rec[2*j-1]<-tr$suppo[2*j-1]+tr$low[1,j]*tr$step[j]
     rec[2*j]<-  tr$suppo[2*j-1]+tr$upp[1,j]*tr$step[j]
  }
  center=t(cenone(rec))
}

else{

nodenumOfTree<-length(tr$left)

# make parent
parent<-makeparent(tr$left,tr$right)

mi<-makeinfo(tr$left,tr$right,tr$mean,tr$low,tr$upp)
#infopointer<-mi$infopointer
#terminalnum<-mi$terminalnum
#low<-mi$low
#upp<-mi$upp
#nodefinder<-mi$nodefinder
#value<-mi$value

{
if (!is.null(Q)){
   maxval<-max(mi$value)
   minval<-min(mi$value)
   step<-(maxval-minval)/Q
   levseq<-seq(from=minval,to=maxval-step,by=step)
}
else{
   eppsi<-0        #0.0000001
   levseq<-matrix(0,length(mi$value),1)
   ordu<-order(mi$value)
   ru<-1
   laskuri<-1
   car<-ordu[ru]
   levseq[laskuri]<-mi$value[car]-eppsi
   while (ru < length(mi$value)){
       carnew<-ordu[ru+1]
       if (mi$value[carnew]>mi$value[car]){
          laskuri<-laskuri+1
          levseq[laskuri]<-mi$value[carnew]-eppsi
      }
      ru<-ru+1
   }
   levseq<-levseq[1:laskuri]
   Q<-laskuri
}
}

levfrekv<-matrix(0,Q,1)
atomnum<-length(mi$value)   #=mi$terminalnum
for (i in 1:atomnum){
   for (j in 1:Q){
      if (mi$value[i]>=levseq[j]){
         levfrekv[j]<-levfrekv[j]+1
      }
   }
}
numofall<-sum(levfrekv)

inlevseq<-matrix(0,Q+1,1)
inlevseq[2:(Q+1)]<-levseq
insuppo<-matrix(0,2*d+1,1)
insuppo[2:(2*d+1)]<-tr$suppo
instep<-matrix(0,d+1,1)
sc<-matrix(0,d,1)
for (i in 1:d){
    step[i]<-(tr$support[2*i]-tr$support[2*i-1])/tr$N[i]
}
instep[2:(d+1)]<-sc    #stepcalc(tr$support,tr$N)    #tr$step
inleft<-matrix(0,nodenumOfTree+1,1)
inleft[2:(nodenumOfTree+1)]<-tr$left
inright<-matrix(0,nodenumOfTree+1,1)
inright[2:(nodenumOfTree+1)]<-tr$right
inparent<-matrix(0,nodenumOfTree+1,1)
inparent[2:(nodenumOfTree+1)]<-parent
inval<-matrix(0,nodenumOfTree+1,1)
inval[2:(nodenumOfTree+1)]<-tr$mean  #tr$val
invec<-matrix(0,nodenumOfTree+1,1)
invec[2:(nodenumOfTree+1)]<-tr$direc

for (i in 1:(nodenumOfTree+1)){
  if (is.na(inval[i])){
       inval[i]<-0
       invec[i]<-0
  }
}

ininfopointer<-matrix(0,nodenumOfTree+1,1)
ininfopointer[2:(nodenumOfTree+1)]<-mi$infopointer

invalue<-matrix(0,atomnum+1,1)
invalue[2:(atomnum+1)]<-mi$value
inlow<-matrix(0,atomnum*d+1,1)
inupp<-matrix(0,atomnum*d+1,1)
for (i in 1:atomnum){
   for (j in 1:d){
       inlow[1+(i-1)*d+j]=mi$low[i,j]
       inupp[1+(i-1)*d+j]=mi$upp[i,j]
   }
}
innodefinder<-matrix(0,atomnum+1,1)
innodefinder[2:(atomnum+1)]<-mi$nodefinder

inlowtr<-matrix(0,nodenumOfTree*d+1,1)
inupptr<-matrix(0,nodenumOfTree*d+1,1)
for (i in 1:nodenumOfTree){
   for (j in 1:d){
       inlowtr[1+(i-1)*d+j]=tr$low[i,j]
       inupptr[1+(i-1)*d+j]=tr$upp[i,j]
   }
}

# we have tree with "nodenumOfTree" nodes
# we hae assocoated info with "atomnum" elements => info for each leaf
#   that is, atomnum = number of leaves

dentree<-.C("proftreeC",
               as.integer(numofall),
               as.integer(atomnum),
               as.double(inlevseq),
               as.integer(d),
               as.integer(Q),
               as.double(instep),
               as.double(insuppo), 
               as.integer(nodenumOfTree),
               as.integer(inleft),
               as.integer(inright),
               as.integer(inparent),
               as.integer(inval),
               as.integer(invec),
               as.integer(ininfopointer),
               as.integer(inlowtr),
               as.integer(inupptr),
               as.double(invalue),
               as.integer(inlow),
               as.integer(inupp),
               as.integer(innodefinder),
               level = double(numofall+1),
               parent = integer(numofall+1),
               volume = double(numofall+1),
               center = double(d*numofall+1),
               efek = integer(1),
PACKAGE="denpro")
               #component = integer(numofall+1),
               #AtomlistAtomOut = integer(numofall+1),
               #AtomlistNextOut = integer(numofall+1),
               #numOfAtoms = integer(1))

#               lapu = double(numofall+1))

efek<-dentree$efek
numOfAtoms<-dentree$numOfAtoms
invalue<-dentree$level[2:(efek+1)]
parent<-dentree$parent[2:(efek+1)]
#component<-dentree$component[2:(efek+1)]
#AtomlistAtom<-dentree$AtomlistAtom[2:(numOfAtoms+1)]
#AtomlistNext<-dentree$AtomlistNext[2:(numOfAtoms+1)]

#if (cfre) nodefrek<-cfrekvdya(seplsets,binfrek) else nodefrek<-NULL 
if (ccen==TRUE) cvol<-TRUE
if (cvol){
#  volume<-cvolumbag(component=component,AtomlistAtom=AtomlistAtom,AtomlistNext=AtomlistNext,low=tr$low,upp=tr$upp,steppi=tr$step)
   volume<-dentree$volume[2:(efek+1)]
#  kerroin<-cinte(invalue,volume,parent) 
#  sepvalnor<-invalue/kerroin
   sepvalnor<-invalue
} 
else{  
  volume<-NULL
  sepvalnor<-NULL
}

if (ccen && cvol){
  #center<-ccentebag(component,AtomlistAtom,AtomlistNext,tr$low,tr$upp,volume,
  #                  tr$step,tr$suppo)
  outcenter<-dentree$center[2:(d*efek+1)]
  center<-matrix(0,efek,d)
  for (i in 1:efek){
     for (j in 1:d){
        center[i,j]<-outcenter[(i-1)*d+j]
     }
  }
  }
  else{
      center<-NULL
  }


} #else (tr$left[1]>0)


return(list(parent=parent,level=sepvalnor,invalue=invalue,
volume=volume,center=t(center)))    #nodefrek=nodefrek))
#values: normeeratut arvot
#invalues: alkuperaiset frekvenssit/arvot 
#nodefrek: kunkin solmun frekvenssi
}










