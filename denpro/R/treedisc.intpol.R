treedisc.intpol<-function(lst, f, ngrid=NULL, r=NULL, lowest="dens")
{
# r is vector of radiuses, we prune shapetree "lst" so that
# its radiuses are given by r

if (lowest=="dens") lowest<-0 else lowest<-min(lst$level)

if (is.null(r)){
      stepsi<-(lst$maxdis-lowest)/(ngrid+1)    
      r<-seq(lowest+stepsi,lst$maxdis-stepsi,stepsi)
}

mt<-multitree(lst$parent)
child<-mt$child
sibling<-mt$sibling

d<-dim(lst$center)[1]
itemnum<-length(lst$parent)

################################################

parent<-matrix(NA,itemnum,1)

rootnum<-length(mt$root)
for (rr in 1:rootnum){

pino<-matrix(0,itemnum,1)
pinoparent<-matrix(0,itemnum,1)
pinorad<-matrix(0,itemnum,1)

pino[1]<-mt$root[rr]
pinoparent[1]<-0
pinorad[1]<-1
pinin<-1
curradind<-1

while (pinin>0){ # && (curradind<=length(r))){
      cur<-pino[pinin]      #take from stack
      curpar<-pinoparent[pinin]
      curradind<-pinorad[pinin]
      pinin<-pinin-1

      # put to the stack
      if (sibling[cur]>0){
            pinin<-pinin+1
            pino[pinin]<-sibling[cur]
            pinoparent[pinin]<-curpar
            pinorad[pinin]<-curradind
      }

      note<-lst$infopointer[cur] #cur
      
      etai<-f[note]

      if (curradind<=length(r)) currad<-r[curradind] else currad<-1000000
      if (etai>currad){
          parent[cur]<-curpar
          curpar<-cur
          curradind<-curradind+1
      }

      # go to left and put right nodes to the stack
      while (child[cur]>0){  # && (curradind<=length(r))){
            cur<-child[cur]

            if (sibling[cur]>0){
                 pinin<-pinin+1
                 pino[pinin]<-sibling[cur]
                 pinoparent[pinin]<-curpar
                 pinorad[pinin]<-curradind
            }
 
            note<-lst$infopointer[cur] #cur

            etai<-f[note]

            if (curradind<=length(r)) currad<-r[curradind] else currad<-1000000
            if (etai>currad){
                parent[cur]<-curpar
                curpar<-cur
                curradind<-curradind+1 
            }
 
      }

}

}

# Prune ##################################

newparent<-matrix(0,itemnum,1)
newcenter<-matrix(0,d,itemnum)
newvolume<-matrix(0,itemnum,1)
newlevel<-matrix(0,itemnum,1)
newpointer<-matrix(0,itemnum,1)
newdistcenter<-matrix(0,d,itemnum)
newproba<-matrix(0,itemnum,1)

i<-1
newlkm<-0
while (i<=itemnum){
    if (!is.na(parent[i])){
         newlkm<-newlkm+1
         newpointer[i]<-newlkm
         if (parent[i]==0)  newparent[newlkm]<-0
         else newparent[newlkm]<-newpointer[parent[i]]
         newcenter[,newlkm]<-lst$center[,i]
         newlevel[newlkm]<-lst$level[i]
         newvolume[newlkm]<-lst$volume[i]
         #newdistcenter[,newlkm]<-lst$distcenter[,i]
         #newproba[newlkm]<-lst$proba[i]
    }
    i<-i+1
}

newparent<-newparent[1:newlkm]
if (newlkm<=1) newcenter<-matrix(newcenter[,1],d,1) 
else newcenter<-newcenter[,1:newlkm]
newvolume<-newvolume[1:newlkm]
newlevel<-newlevel[1:newlkm]
if (newlkm<=1) newdistcenter<-matrix(newdistcenter[,1],d,1) 
else newdistcenter<-newdistcenter[,1:newlkm]
newproba<-newproba[1:newlkm]
newpointer<-newpointer[1:newlkm]

return(list(parent=newparent,level=newlevel,volume=newvolume,center=newcenter,
distcenter=newdistcenter,  #branchradius=newbranchradius,
proba=newproba,
refe=lst$refe,bary=lst$bary,root=1,infopointer=newpointer))

}   

