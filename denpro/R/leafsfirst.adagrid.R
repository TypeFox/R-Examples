leafsfirst.adagrid<-function(pcf)
{
down<-pcf$down
high<-pcf$high
support<-pcf$support
grid<-pcf$grid
value<-pcf$value

d<-dim(down)[2]
lkm<-dim(down)[1]

distat<-pcf$value
infopointer<-seq(1,lkm)

# order the atoms for the level set with level "lev"

ord<-order(distat)
infopointer<-infopointer[ord]

# create tree

parent<-matrix(0,lkm,1)
child<-matrix(0,lkm,1)
sibling<-matrix(0,lkm,1)
volume<-matrix(0,lkm,1)
radius<-matrix(0,lkm,1)
number<-matrix(0,lkm,1)
atomlist<-matrix(0,lkm,lkm)
atomnumb<-matrix(0,lkm,1)
center<-matrix(0,d,lkm)
distcenter<-matrix(0,lkm,d)

highestNext<-matrix(0,lkm,1)    #pointers to the nodes without parent
boundrec<-matrix(0,lkm,2*d) #for each node, the box which bounds all the c:dren

node<-lkm  #ord[lkm]  #the 1st child node is the one with the longest distance
parent[node]<-0
child[node]<-0
sibling[node]<-0

# radius
radius[node]<-distat[ord[node]]

downi<-down[infopointer[node],]
highi<-high[infopointer[node],]
simp<-matrix(0,2*d,1)
for (i in 1:d){
    simp[2*i-1]<-grid[downi[i],i]
    simp[2*i]<-grid[highi[i],i]
}
volume[node]<-massone(simp)
for (dd in 1:d){
   center[dd,node]<-(simp[2*dd-1]+simp[2*dd])/2
}
number[node]<-1
atomlist[node,1]<-infopointer[node]
atomnumb[node]<-1

beg<-node                 #first without parent
highestNext[node]<-0
note<-infopointer[node]   #note<-pcf$nodefinder[infopointer[node]]
for (i in 1:d){
  boundrec[node,2*i-1]<-down[note,i]   
  boundrec[node,2*i]<-high[note,i]  
}

j<-2
while (j<=lkm){
    node<-lkm-j+1   #ord[lkm-j+1]

    # lisaa "node" ensimmaiseksi listaan
    highestNext[node]<-beg  #beg on listan tamanhetkinen ensimmainen
    beg<-node           

    # add node-singleton to boundrec
    rec1<-matrix(0,2*d,1)    #luo sigleton
    note<-infopointer[node]  #note<-pcf$nodefinder[infopointer[node]]
    for (i in 1:d){
         rec1[2*i-1]<-down[note,i]  
         rec1[2*i]<-high[note,i] 
    }
    boundrec[node,]<-rec1

    # radius
    radius[node]<-distat[ord[node]]

    downi<-down[infopointer[node],]
    highi<-high[infopointer[node],]
    simp<-matrix(0,2*d,1)
    for (i in 1:d){
       simp[2*i-1]<-grid[downi[i],i]
       simp[2*i]<-grid[highi[i],i]
    }
    volume[node]<-massone(simp) 
    for (dd in 1:d){
       center[dd,node]<-(simp[2*dd-1]+simp[2*dd])/2
    } 
    number[node]<-1
    atomlist[node,1]<-infopointer[node]
    atomnumb[node]<-1

    curroot<-highestNext[beg]  #node on 1., listassa ainakin 2
    prevroot<-beg
    ekatouch<-0
    while (curroot>0){
       istouch<-touchstep(node,curroot,boundrec,child,sibling,infopointer,down,high,)

        if (istouch==1){

           # paivita parent, child, sibling, volume 
           parent[curroot]<-node           
           if (ekatouch==0) ekatouch<-1 else ekatouch<-0 
           if (ekatouch==1){
              child[node]<-curroot
           }
           else{  # since ekatouch==0, prevroot>0
              sibling[lastsib]<-curroot
           }

           number[node]<-number[node]+number[curroot]
           center[,node]<-(center[,node]*volume[node]+center[,curroot]*volume[curroot])/(volume[node]+volume[curroot])
           volume[node]<-volume[node]+volume[curroot]


           atomlist[node,(atomnumb[node]+1):(atomnumb[node]+atomnumb[curroot])]<-atomlist[curroot,1:atomnumb[curroot]]
           atomnumb[node]<-atomnumb[curroot]+atomnumb[node]

           radius[node]<-min(distat[ord[node]],distat[ord[curroot]])

           # attach box of curroot
           rec1<-boundrec[node,]
           rec2<-boundrec[curroot,] 
           boundrec[node,]<-boundbox(rec1,rec2)
           # poista "curroot" listasta 
           highestNext[prevroot]<-highestNext[curroot]

        }     
        # if curroot was not removed, we update prevroot
        # else curroot was removed, we update lastsib
        if (istouch==0) prevroot<-curroot else lastsib<-curroot 
        curroot<-highestNext[curroot]
    }
    j<-j+1
}

root<-1 #ord[1]  #root is the barycenter

maxdis<-distat[ord[length(ord)]]
#center<-matrix(0,d,lkm)
#for (i in 1:lkm){
#    for (j in 1:d){
#       ala<-grid[down[infopointer[i],j],j]
#       yla<-grid[high[infopointer[i],j],j]
#       center[j,i]<-(ala+yla)/2
#    }
#}

lf<-list(
parent=parent,volume=volume,center=center,level=radius,
root=root,
infopointer=infopointer,
maxdis=maxdis,
dendat=pcf$dendat,
atomlist=atomlist,atomnumb=atomnumb)

return(lf)
}


