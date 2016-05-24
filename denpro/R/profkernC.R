profkernC<-function(dendat,h,N,Q,cvol=TRUE,ccen=TRUE,#cfre=FALSE,
numofallR=10000){

#set.seed(seed=1)
#dendat<-matrix(rnorm(20),10)
#h<-1 
#N<-c(8,8)
#Q<-3

n<-dim(dendat)[1]
d<-length(N)
hnum<-length(h)
mnn<-maxnodenum(dendat,h,N,n,d)
extMaxnode<-mnn$maxnode
extMaxvals<-mnn$maxpositive

if (hnum>1){
 inh<-matrix(0,hnum+1,1)
 inh[2:(hnum+1)]<-h
}
else{
 inh<-h
}
inN<-matrix(0,d+1,1)
inN[2:(d+1)]<-N

dentree<-.C("kerprofC",as.integer(extMaxnode),
                  as.integer(extMaxvals),
                  as.double(dendat),
                  as.double(inh),
                  as.integer(inN),
                  as.integer(n),
                  as.integer(hnum),
                  as.integer(d),
                  as.integer(Q),
                  as.integer(numofallR),
                  level = double(numofallR+1),
                  parent = integer(numofallR+1),
                  component = integer(numofallR+1),
                  volume = double(numofallR+1),
                  center = double(d*numofallR+1),
                  efek = integer(1),
PACKAGE="denpro")

invalue<-dentree$level[2:(dentree$efek+1)]
parent<-dentree$parent[2:(dentree$efek+1)]
volume<-dentree$volume[2:(dentree$efek+1)]
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

return(list(parent=parent,level=sepvalnor,invalue=invalue,
volume=volume,center=center))
#,nodefrek=nodefrek,clustervec=clustervecs))
#
#values: normeeratut arvot
#invalues: normeeraamattomat arvot 
#nodefrek: kunkin solmun frekvenssi
}









