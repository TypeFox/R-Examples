prof2vecs<-function(profile,level,n=NULL,crit,motes=NULL){

parents<-profile$parent
nodenum<-length(parents)
centers<-profile$center
 
nodenum<-length(parents)   
levels<-matrix(level,nodenum,1) #all will be plotted at same lev(=logh)
excma<-excmas(profile)       #instead of volumes, we use excesss mass
                             #to determine the lengths of the vectors
#motes<-mtest(profile,n)

mut<-multitree(parents)

# let us make a vector where modes are labelled with the order, others=0
# later we handle "mlabel" similarily as "motes"
mlabel<-matrix(0,nodenum,1)
mlkm<-moodilkm(parents)      #mlkm$lkm, mlkm$modloc 
for (run in 1:mlkm$lkm){
   alku<-mlkm$modloc[run]
   while ((parents[alku]>0) && 
          (mut$sibling[mut$child[parents[alku]]]==0)){
      alku<-parents[alku]
   }
   mlabel[alku]<-run  
}

mt<-pruneprof(mut)
depths<-depth(mt)
roots<-mt$roots
child<-mt$child
sibling<-mt$sibling

sibord<-siborder(mt,crit,centers)

itemnum<-length(parents)
vecs<-matrix(NA,itemnum,4)  
vecs<-alloroot(vecs,roots,sibord,levels,excma) 
vecs<-plotdata(roots,child,sibling,sibord,levels,excma,vecs)
vecnum<-length(vecs[,1])      #vecs has four columns

#  remove pruned

if (is.null(motes)) motes<-matrix(0,vecnum,1)

tempvecs<-matrix(0,vecnum,4)
tempdepths<-matrix(0,vecnum,1)
tempmotes<-matrix(0,vecnum,1)
tempmlabel<-matrix(0,vecnum,1)
ind<-0
for (i in 1:vecnum){
       if (!(is.na(vecs[i,1]))){
             ind<-ind+1
             tempvecs[ind,]<-vecs[i,]
             tempdepths[ind]<-depths[i]
             tempmotes[ind]<-motes[i]
             tempmlabel[ind]<-mlabel[i]
         }
}
vecs<-tempvecs[1:ind,]
depths<-tempdepths[1:ind]
motes<-tempmotes[1:ind]
mlabel<-tempmlabel[1:ind]

return(list(vecs=vecs,depths=depths,motes=motes,mlabel=mlabel))
}                        





