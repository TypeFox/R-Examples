exmap<-function(estiseq,mt,hind=NULL,hseq=NULL,
n=NULL,moteslist=NULL,ylist=NULL)
{
#moteslist is a list of alpha values for every node
#not just for the branch nodes, but it might be nonsense for others

pk<-estiseq$lstseq
if (is.null(hseq)) hseq<-mt$hseq
if (is.null(hind)) hind<-c(1:length(mt$hseq))
slis<-mt$hseq[hind]

d<-dim(pk[[1]]$center)[1]

if (is.null(ylist)) ylist<-c(length(slis):1)

hrun<-1
for (i in 1:length(slis)){
   while (hseq[hrun]!=slis[i]){
      hrun<-hrun+1
   }
   if (i==1) treelist<-list(pk[[hrun]])
   else  treelist=c(treelist,list(pk[[hrun]]))
}

parnum<-length(slis)
veclkm<-0

if (d==1){
  crit<-max(treelist[[1]]$center)
}
else{
  crit<-rep(0,d)
}

for (i in 1:parnum){
     scur<-slis[i]

     if (!is.null(ylist))  yhigh<-ylist[i]
     else yhigh<-scur     

     profile<-treelist[[i]]
     
     if (!is.null(moteslist))  motes<-moteslist[[i]]
     else motes<-NULL

     level<-scur
     levelplot<-yhigh

     vecplu<-prof2vecs(profile,levelplot,n,crit,motes=motes)
     vecs<-vecplu$vecs
     depths<-vecplu$depths
     motes<-vecplu$motes
     mlabel<-vecplu$mlabel
     vecnum<-length(depths)
     smoot<-matrix(level,vecnum,1)
     
     # concatenate to big's
     
     veclkmold<-veclkm
     veclkm<-veclkm+vecnum
     if (veclkmold==0){   
        bigvecs<-vecs
        bigdepths<-depths
        bigmotes<-motes
        bigmlabel<-mlabel
        bigsmoot<-smoot
     }
     else{
       temp<-matrix(0,veclkm,4)
       temp[1:veclkmold,]<-bigvecs
       temp[(veclkmold+1):veclkm,]<-vecs
       bigvecs<-temp
       
       tempdep<-matrix(0,veclkm,1)
       tempdep[1:veclkmold]<-bigdepths
       tempdep[(veclkmold+1):veclkm]<-depths
       bigdepths<-tempdep
       
       tempmoo<-matrix(0,veclkm,1)
       tempmoo[1:veclkmold]<-bigmotes
       tempmoo[(veclkmold+1):veclkm]<-motes
       bigmotes<-tempmoo
       
       templab<-matrix(0,veclkm,1)
       templab[1:veclkmold]<-bigmlabel
       templab[(veclkmold+1):veclkm]<-mlabel
       bigmlabel<-templab
       
       tempsmoo<-matrix(0,veclkm,1)
       tempsmoo[1:veclkmold]<-bigsmoot
       tempsmoo[(veclkmold+1):veclkm]<-smoot
       bigsmoot<-tempsmoo
     }
}    
#if (makeplot==T) plotvecs(bigvecs,segme=T,bigdepths) 

return(list(bigvecs=bigvecs,bigdepths=t(bigdepths),motes=t(bigmotes),mlabel=t(bigmlabel),smoot=t(bigsmoot)))
}







