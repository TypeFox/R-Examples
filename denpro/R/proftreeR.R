proftreeR<-function(tr,
Q=NULL,frekv=NULL,cvol=TRUE,ccen=TRUE,cfre=FALSE)
{
#set.seed(seed=1)
#dendat<-matrix(rnorm(20),10)
#h<-1 
#N<-c(4,4)
#Q<-3

d<-dim(tr$upp)[2]

nodenumOfTree<-length(tr$left)

low<-tr$low
upp<-tr$upp
val<-tr$val
#low<-matrix(0,nodenumOfTree,d)
#upp<-matrix(0,nodenumOfTree,d)
#val<-matrix(NA,nodenumOfTree,1)
#for (i in 1:nodenumOfTree){
#  dimu<-tr$vec[i]
#  if (!is.na(dimu) && (dimu>0)) 
#       val[i]<-tr$suppo[2*dimu-1]+tr$val[i]*tr$step[dimu]
#  for (j in 1:d){
#      low[i,j]<-tr$suppo[2*j-1]+tr$low[i,j]*tr$step[j]
#      upp[i,j]<-tr$suppo[2*j-1]+tr$upp[i,j]*tr$step[j]
#  }
#}

# make parent
parent<-matrix(0,length(tr$left),1)
node<-1
while (node<=length(tr$left)){
   if ((!is.na(tr$left[node])) && (tr$left[node]!=0)){
        parent[tr$left[node]]<-node
   }
   if ((!is.na(tr$right[node])) && (tr$left[node]!=0)){
        parent[tr$right[node]]<-node
   }
   node<-node+1
}

mi<-makeinfo(tr$left,tr$right,tr$mean,low,upp)
infopointer<-mi$infopointer
terminalnum<-mi$terminalnum
low<-mi$low
upp<-mi$upp
nodefinder<-mi$nodefinder
value<-mi$value

{
if (!is.null(Q)){
   maxval<-max(value)
   minval<-min(value)
   step<-(maxval-minval)/Q
   levseq<-seq(from=minval,to=maxval-step,by=step)
}
else{
   eppsi<-0        #0.0000001
   levseq<-matrix(0,length(value),1)
   ordu<-order(value)
   ru<-1
   #car<-ordu[ru]
   #while ((ru <= length(value)) && (value[car]==0)){
   #     ru<-ru+1
   #     car<-ordu[ru]
   #}  # we have found first non zero
   laskuri<-1
   car<-ordu[ru]
   levseq[laskuri]<-value[car]-eppsi
   while (ru < length(value)){
       carnew<-ordu[ru+1]
       if (value[carnew]>value[car]){
          laskuri<-laskuri+1
          levseq[laskuri]<-value[carnew]-eppsi
      }
      ru<-ru+1
   }
   levseq<-levseq[1:laskuri]
   Q<-laskuri
}
}

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

dentree<-decombag(numofall,levseq,
tr$left,tr$right,val,tr$vec,infopointer,parent,
nodenumOfTree,terminalnum,
value,low,upp,nodefinder,
d)

invalue<-dentree$level
parent<-dentree$parent
component<-dentree$component
AtomlistAtom<-dentree$AtomlistAtom
AtomlistNext<-dentree$AtomlistNext

#if (cfre) nodefrek<-cfrekvdya(seplsets,binfrek) else nodefrek<-NULL 

if (ccen==TRUE) cvol<-TRUE
if (cvol){
  volume<-cvolumbag(component,AtomlistAtom,AtomlistNext,tr$low,tr$upp,
                    steppi=tr$step)
  kerroin<-cinte(invalue,volume,parent) 
  sepvalnor<-invalue/kerroin
} 
else{  
  volume<-NULL
  sepvalnor<-NULL
}

if (ccen && cvol){
  center<-ccentebag(component,AtomlistAtom,AtomlistNext,tr$low,tr$upp,volume,
                    tr$step,tr$suppo)
  }
  else{
      center<-NULL
  }

return(list(parent=parent,level=sepvalnor,invalue=invalue,
volume=volume,center=center))    #nodefrek=nodefrek))
#values: normeeratut arvot
#invalues: alkuperaiset frekvenssit/arvot 
#nodefrek: kunkin solmun frekvenssi
}









